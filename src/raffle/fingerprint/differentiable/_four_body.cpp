#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <thread>
#include <vector>

template <typename F> static void parallel(int n, int workers, F fn) {
    if (workers <= 1) {
        for (int i = 0; i < n; ++i) fn(i);
        return;
    }
    std::atomic<int> next{0};
    std::vector<std::thread> threads;
    try {
        for (int t = 0; t < std::min(workers, n); ++t)
            threads.emplace_back([&] {
                for (;;) {
                    int i = next.fetch_add(1, std::memory_order_relaxed);
                    if (i >= n) break;
                    fn(i);
                }
            });
    } catch (...) {
        for (auto &t : threads) t.join();
        throw;
    }
    for (auto &t : threads) t.join();
}

// Each environment owns its histogram and gradient buffers. Metadata gives
// the offsets of its pair normals, fourth normals, near and fourth neighbors,
// followed by their counts and a single-species marker.
//
// The Gaussian is evaluated on a uniform grid. Writing the offset from the
// nearest bin as m and the residual as s, the kernel obeys
//
//     exp(-c (s - m)^2) = exp(-c s^2) * exp(2 c s)^m * exp(-c m^2),
//
// so only two transcendental calls are needed per contribution; the per-bin
// values follow from two short geometric recurrences. ``kernel_table`` holds
// exp(-c m^2) for every offset within the truncated support.
//
// Reciprocal distances and per-pair weight factors are hoisted out of the
// fourth-neighbour loop. Backward accumulators live in per-environment
// scratch buffers so the many contributions to a shared normal or distance
// do not repeatedly write the same global element.
//
// ``kernels`` is addressed by signed bin offset from ``nearest``. Both sides
// of the support use four independent geometric chains so the dependent
// multiply latency does not serialise the whole recurrence.
template <typename T>
static inline void fill_kernels(T *kernels, const T *table, int radius,
                                int lower, int upper, T first, T ratio) {
    kernels[0] = first * table[radius];
    const T ratio2 = ratio * ratio;
    const T ratio4 = ratio2 * ratio2;
    T p0 = first * ratio, p1 = first * ratio2;
    T p2 = p1 * ratio, p3 = first * ratio4;
    int offset = 1;
    for (; offset + 3 <= upper; offset += 4) {
        kernels[offset] = p0 * table[offset + radius];
        kernels[offset + 1] = p1 * table[offset + 1 + radius];
        kernels[offset + 2] = p2 * table[offset + 2 + radius];
        kernels[offset + 3] = p3 * table[offset + 3 + radius];
        p0 *= ratio4; p1 *= ratio4; p2 *= ratio4; p3 *= ratio4;
    }
    for (; offset <= upper; ++offset) {
        kernels[offset] = p0 * table[offset + radius];
        p0 *= ratio;
    }
    const T inverse = T(1) / ratio;
    const T inverse2 = inverse * inverse;
    const T inverse4 = inverse2 * inverse2;
    T q0 = first * inverse, q1 = first * inverse2;
    T q2 = q1 * inverse, q3 = first * inverse4;
    offset = -1;
    for (; offset - 3 >= lower; offset -= 4) {
        kernels[offset] = q0 * table[offset + radius];
        kernels[offset - 1] = q1 * table[offset - 1 + radius];
        kernels[offset - 2] = q2 * table[offset - 2 + radius];
        kernels[offset - 3] = q3 * table[offset - 3 + radius];
        q0 *= inverse4; q1 *= inverse4; q2 *= inverse4; q3 *= inverse4;
    }
    for (; offset >= lower; --offset) {
        kernels[offset] = q0 * table[offset + radius];
        q0 *= inverse;
    }
}
template <typename T, bool backwards>
static void run_impl(int n, int bins, int workers,
                     const int64_t *meta, const int64_t *left, const int64_t *right,
                     const T *a, const T *na, const T *b, const T *nb,
                     const T *d, const T *fd, const T *win, const T *cp,
                     const T *wf, const T *wr, const T *den, const T *grid,
                     T eta, T width, T cutoff, int radius, T epsilon,
                     T amplitude_factor, const T *kernel_table, T *out,
                     T *ga, T *gna, T *gb, T *gnb, T *gd, T *gfd, T *gwin,
                     T *gcp, T *gwf, T *gwr, T *gden, const T *upstream) {
    const int span = 2 * radius + 1;
    const T two_eta = T(2) * eta;
    const T c = eta * width * width;
    const T tiny = std::numeric_limits<T>::min();
    // ``grid`` is uniform by construction, so only its spacing (``width``) and
    // origin (``cutoff``) are needed by the recurrence.
    (void)grid;
    parallel(n, workers, [&](int e) {
        const int64_t *m = meta + 8 * e;
        const int64_t po = m[0], bo = m[1], no = m[2], fo = m[3];
        const int pairs = int(m[4]), near = int(m[5]), nf = int(m[6]);
        const int channel = int(m[7]);
        const T denominator = std::max(den[e], tiny);
        const T amplitude = amplitude_factor / denominator;
        const T center_probability = cp[e];
        const int orientations = channel == -1 ? 1 : 2;
        T *out_env = out + int64_t(e) * bins;
        const T *upstream_env = backwards ? upstream + int64_t(e) * bins : nullptr;
        std::vector<T> kernels(span);
        std::vector<T> inverse_near(near), inverse_fourth(nf);
        std::vector<T> weighted_win(nf), inverse_weighted_win(nf), weighted_cp(nf);
        // Accumulators shared by many fourth-neighbour contributions.
        std::vector<T> sum_gwin, sum_gfd, sum_gd;
        T sum_gcp = T(0), sum_gden = T(0);
        if constexpr (backwards) {
            sum_gwin.assign(nf, T(0));
            sum_gfd.assign(nf, T(0));
            sum_gd.assign(near, T(0));
        }
        for (int i = 0; i < near; ++i)
            inverse_near[i] = T(1) / d[no + i];
        for (int i = 0; i < nf; ++i)
            inverse_fourth[i] = T(1) / fd[fo + i];
        for (int p = 0; p < pairs; ++p) {
            const int j = int(left[po + p]), k = int(right[po + p]);
            const T pair_distance_sq = d[no + j] * d[no + j] * d[no + k] * d[no + k];
            // Per-fourth-neighbour factors shared by both orientations.
            for (int l = 0; l < nf; ++l) {
                const T product = pair_distance_sq * fd[fo + l] * fd[fo + l];
                const T inv = T(1) / std::max(product, T(1e-24));
                const T win_inv = win[fo + l] * inv;
                weighted_win[l] = win_inv * center_probability;
                inverse_weighted_win[l] = win_inv;
                weighted_cp[l] = center_probability * inv;
            }
            const int64_t pi = po + p;
            const T *a_p = a + 3 * pi;
            for (int rev = 0; rev < orientations; ++rev) {
                const T norm_a = na[pi];
                if (norm_a <= T(1e-12)) continue;
                const T sign = rev ? T(-1) : T(1);
                const T weight_pair = rev ? wr[pi] : wf[pi];
                const int second = rev ? j : k;
                const T *b_base = b + 3 * int64_t(bo + int64_t(second) * nf);
                const T *nb_base = nb + bo + int64_t(second) * nf;
                T pair_gw = T(0), pair_gna = T(0);
                T pair_ga[3] = {T(0), T(0), T(0)};
                for (int l = 0; l < nf; ++l) {
                    const T norm_b = nb_base[l];
                    if (norm_b <= T(1e-12)) continue;
                    const T dl = fd[fo + l];
                    const T product = pair_distance_sq * dl * dl;
                    const int64_t bi = bo + int64_t(second) * nf + l;
                    const T *b_l = b_base + 3 * l;
                    const T dot = sign * (a_p[0] * b_l[0] + a_p[1] * b_l[1]
                                          + a_p[2] * b_l[2]);
                    const T prodnorm = norm_a * norm_b;
                    const T inverse_prodnorm = T(1) / std::max(prodnorm, T(1e-24));
                    const T cosang = dot * inverse_prodnorm;
                    const T clamped = std::max(T(-1) + epsilon,
                                               std::min(T(1) - epsilon, cosang));
                    const T angle = std::acos(clamped);
                    const T weight = weight_pair * weighted_win[l];
                    const T scaled = (angle - cutoff) / width;
                    const int nearest = int(std::floor(scaled + T(0.5)));
                    const T residual = scaled - T(nearest);
                    const int lo = std::max(0, nearest - radius);
                    const int hi = std::min(bins - 1, nearest + radius);
                    const T first = std::exp(-c * residual * residual) * amplitude;
                    const T ratio = std::exp(T(2) * c * residual);
                    T *kernels_shifted = kernels.data() + radius;
                    fill_kernels(kernels_shifted, kernel_table, radius,
                                 lo - nearest, hi - nearest, first, ratio);
                    if constexpr (backwards) {
                        // The kernel-weighted sums share factors, so accumulate
                        // them in registers and update the histograms once.
                        T sum_weight = T(0), sum_weight_offset = T(0);
                        for (int q = lo; q <= hi; ++q) {
                            const int offset = q - nearest;
                            const T uk = upstream_env[q] * kernels_shifted[offset];
                            sum_weight += uk;
                            sum_weight_offset += uk * T(offset);
                        }
                        const T grad_weight = sum_weight;
                        const T grad_angle = -two_eta * width * weight
                                             * (residual * sum_weight - sum_weight_offset);
                        if (den[e] >= tiny)
                            sum_gden -= weight * sum_weight / denominator;
                        pair_gw += grad_weight * weighted_win[l];
                        sum_gwin[l] += grad_weight * weight_pair * weighted_cp[l];
                        sum_gcp += grad_weight * weight_pair * inverse_weighted_win[l];
                        if (product >= T(1e-24)) {
                            const T gr = -T(2) * grad_weight * weight;
                            sum_gd[j] += gr * inverse_near[j];
                            sum_gd[k] += gr * inverse_near[k];
                            sum_gfd[l] += gr * inverse_fourth[l];
                        }
                        if (cosang >= T(-1) + epsilon && cosang <= T(1) - epsilon) {
                            const T gc = -grad_angle
                                         / std::sqrt(T(1) - clamped * clamped);
                            for (int axis = 0; axis < 3; ++axis) {
                                pair_ga[axis] += gc * sign * b_l[axis] * inverse_prodnorm;
                                gb[3 * bi + axis] += gc * sign * a_p[axis] * inverse_prodnorm;
                            }
                            if (prodnorm >= T(1e-24)) {
                                pair_gna -= gc * cosang / norm_a;
                                gnb[bi] -= gc * cosang / norm_b;
                            }
                        }
                    } else {
                        const int count = hi - lo + 1;
                        const T *kernel_ptr = kernels_shifted + (lo - nearest);
                        T *out_ptr = out_env + lo;
                        for (int i = 0; i < count; ++i)
                            out_ptr[i] += kernel_ptr[i] * weight;
                    }
                }
                if constexpr (backwards) {
                    (rev ? gwr : gwf)[pi] += pair_gw;
                    gna[pi] += pair_gna;
                    for (int axis = 0; axis < 3; ++axis)
                        ga[3 * pi + axis] += pair_ga[axis];
                }
            }
        }
        if constexpr (backwards) {
            for (int l = 0; l < nf; ++l) {
                gwin[fo + l] += sum_gwin[l];
                gfd[fo + l] += sum_gfd[l];
            }
            for (int i = 0; i < near; ++i)
                gd[no + i] += sum_gd[i];
            gcp[e] += sum_gcp;
            gden[e] += sum_gden;
        }
    });
}

template <typename T>
static void run(int n, int bins, int workers, bool backwards,
                const int64_t *meta, const int64_t *left, const int64_t *right,
                const T *a, const T *na, const T *b, const T *nb,
                const T *d, const T *fd, const T *win, const T *cp,
                const T *wf, const T *wr, const T *den, const T *grid,
                double eta_, double width_, double cutoff_, int radius,
                double epsilon_, T *out,
                T *ga, T *gna, T *gb, T *gnb, T *gd, T *gfd, T *gwin,
                T *gcp, T *gwf, T *gwr, T *gden, const T *upstream) {
    const T eta = T(eta_), width = T(width_), cutoff = T(cutoff_);
    const T eps = T(epsilon_);
    const T factor = T(std::sqrt(eta_ / 3.14159265358979323846));
    const T c = eta * width * width;
    std::vector<T> kernel_table(2 * radius + 1);
    for (int offset = -radius; offset <= radius; ++offset)
        kernel_table[offset + radius] = std::exp(-c * T(offset) * T(offset));
    if (backwards)
        run_impl<T, true>(n, bins, workers, meta, left, right, a, na, b, nb, d, fd,
                          win, cp, wf, wr, den, grid, eta, width, cutoff, radius,
                          eps, factor, kernel_table.data(), out, ga, gna, gb, gnb,
                          gd, gfd, gwin, gcp, gwf, gwr, gden, upstream);
    else
        run_impl<T, false>(n, bins, workers, meta, left, right, a, na, b, nb, d, fd,
                           win, cp, wf, wr, den, grid, eta, width, cutoff, radius,
                           eps, factor, kernel_table.data(), out, ga, gna, gb, gnb,
                           gd, gfd, gwin, gcp, gwf, gwr, gden, upstream);
}

#if defined(_WIN32)
#define RAFFLE_EXPORT __declspec(dllexport)
#else
#define RAFFLE_EXPORT __attribute__((visibility("default")))
#endif

extern "C" RAFFLE_EXPORT int raffle_project_four(
    int dtype, int n, int bins, int workers, int backwards,
    const int64_t *meta, const int64_t *left, const int64_t *right,
    const void *a, const void *na, const void *b, const void *nb,
    const void *d, const void *fd, const void *win, const void *cp,
    const void *wf, const void *wr, const void *den, const void *grid,
    double eta, double width, double cutoff, int radius, double epsilon,
    void *out, void *ga, void *gna, void *gb, void *gnb, void *gd,
    void *gfd, void *gwin, void *gcp, void *gwf, void *gwr,
    void *gden, const void *upstream) {
#define INVOKE(T) run<T>(n,bins,workers,backwards,meta,left,right, \
    (const T*)a,(const T*)na,(const T*)b,(const T*)nb,(const T*)d, \
    (const T*)fd,(const T*)win,(const T*)cp,(const T*)wf,(const T*)wr, \
    (const T*)den,(const T*)grid,eta,width,cutoff,radius,epsilon, \
    (T*)out,(T*)ga,(T*)gna,(T*)gb,(T*)gnb,(T*)gd,(T*)gfd,(T*)gwin, \
    (T*)gcp,(T*)gwf,(T*)gwr,(T*)gden,(const T*)upstream)
    try {
        if (dtype == 4) INVOKE(float); else INVOKE(double);
        return 0;
    } catch (...) {
        return 1;
    }
#undef INVOKE
}
