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
    parallel(n, workers, [&](int e) {
        const int64_t *m = meta + 8 * e;
        const int64_t po = m[0], bo = m[1], no = m[2], fo = m[3];
        const int pairs = int(m[4]), nf = int(m[6]), channel = int(m[7]);
        const T denom = std::max(den[e], std::numeric_limits<T>::min());
        for (int p = 0; p < pairs; ++p) {
            int j = int(left[po + p]), k = int(right[po + p]);
            for (int l = 0; l < nf; ++l) {
                T dj = d[no+j], dk = d[no+k], dl = fd[fo+l];
                T product = dj*dj*dk*dk*dl*dl;
                T inv = T(1) / std::max(product, T(1e-24));
                for (int rev = 0; rev < 2; ++rev) {
                    if (rev && channel == -1) continue;
                    int64_t pi = po + p, bi = bo + (rev ? j : k)*nf+l;
                    T norm_a = na[pi], norm_b = nb[bi];
                    if (norm_a <= T(1e-12) || norm_b <= T(1e-12)) continue;
                    T sign = rev ? T(-1) : T(1);
                    T dot = T(0);
                    for (int c = 0; c < 3; ++c) dot += sign*a[3*pi+c]*b[3*bi+c];
                    T prodnorm = norm_a*norm_b;
                    T cosang = dot / std::max(prodnorm, T(1e-24));
                    T clamped = std::max(T(-1)+eps, std::min(T(1)-eps, cosang));
                    T angle = std::acos(clamped);
                    T weight_pair = rev ? wr[pi] : wf[pi];
                    T weight = weight_pair*win[fo+l]*cp[e]*inv;
                    int nearest = int(std::floor((angle-cutoff)/width+T(.5)));
                    int lo = std::max(0, nearest-radius), hi = std::min(bins-1, nearest+radius);
                    T gradweight = 0, gradangle = 0;
                    for (int q = lo; q <= hi; ++q) {
                        T delta = angle-grid[q];
                        T kernel = std::exp(-eta*delta*delta)*factor/denom;
                        if (!backwards) out[e*bins+q] += kernel*weight;
                        else {
                            T u = upstream[e*bins+q];
                            gradweight += u*kernel;
                            gradangle += u*kernel*weight*(-T(2)*eta*delta);
                            if (den[e] >= std::numeric_limits<T>::min())
                                gden[e] -= u*kernel*weight/denom;
                        }
                    }
                    if (!backwards) continue;
                    // Differentiate the projected feature weight and angle;
                    // PyTorch propagates these local derivatives through the
                    // previously calculated normals, distances and species.
                    T common = win[fo+l]*cp[e]*inv;
                    (rev ? gwr : gwf)[pi] += gradweight*common;
                    gwin[fo+l] += gradweight*weight_pair*cp[e]*inv;
                    gcp[e] += gradweight*weight_pair*win[fo+l]*inv;
                    if (product >= T(1e-24)) {
                        T gr = -T(2)*gradweight*weight;
                        gd[no+j] += gr/dj;
                        gd[no+k] += gr/dk;
                        gfd[fo+l] += gr/dl;
                    }
                    if (cosang >= T(-1)+eps && cosang <= T(1)-eps) {
                        T gc = -gradangle/std::sqrt(T(1)-clamped*clamped);
                        T div = std::max(prodnorm, T(1e-24));
                        for (int c = 0; c < 3; ++c) {
                            ga[3*pi+c] += gc*sign*b[3*bi+c]/div;
                            gb[3*bi+c] += gc*sign*a[3*pi+c]/div;
                        }
                        if (prodnorm >= T(1e-24)) {
                            gna[pi] -= gc*cosang/norm_a;
                            gnb[bi] -= gc*cosang/norm_b;
                        }
                    }
                }
            }
        }
    });
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
