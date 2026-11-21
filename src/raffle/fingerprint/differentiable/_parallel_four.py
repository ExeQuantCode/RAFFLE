"""Optional CPU four-body projection and first-order analytical backward."""

import ctypes
from functools import lru_cache
import math
from pathlib import Path
import sys

import torch


@lru_cache(maxsize=1)
def _load_native():
    import raffle

    suffix = ".dll" if sys.platform == "win32" else ".dylib" if sys.platform == "darwin" else ".so"
    name = "_raffle_four_body" + suffix
    candidates = [Path(__file__).with_name(name)] + [
        Path(root) / "fingerprint" / "differentiable" / name for root in raffle.__path__
    ]
    path = next((candidate for candidate in candidates if candidate.is_file()), None)
    if path is None:
        raise RuntimeError(
            "The parallel four-body library is unavailable. Rebuild/install RAFFLE "
            "with Python support, or omit parallel_workers."
        )
    try:
        library = ctypes.CDLL(str(path))
    except OSError as error:
        raise RuntimeError(
            "The parallel four-body library is unavailable. Rebuild/install RAFFLE "
            "with Python support, or omit parallel_workers."
        ) from error
    function = library.raffle_project_four
    pointer = ctypes.c_void_p
    function.argtypes = (
        [ctypes.c_int] * 5 + [pointer] * 15
        + [ctypes.c_double] * 3 + [ctypes.c_int, ctypes.c_double]
        + [pointer] * 13
    )
    function.restype = ctypes.c_int
    return function


def _pointer(tensor):
    return ctypes.c_void_p(tensor.data_ptr()) if tensor is not None else ctypes.c_void_p()


def _call(metadata, left, right, grid, inputs, settings, output, gradients, upstream):
    eta, width, minimum, radius, epsilon, workers = settings
    gradients = gradients if gradients is not None else [None] * 11
    status = _load_native()(
        grid.element_size(), metadata.shape[0], grid.numel(), workers, int(upstream is not None),
        _pointer(metadata), _pointer(left), _pointer(right),
        *(_pointer(value) for value in inputs), _pointer(grid),
        eta, width, minimum, radius, epsilon, _pointer(output),
        *(_pointer(value) for value in gradients), _pointer(upstream),
    )
    if status:
        raise RuntimeError("Parallel four-body calculation failed to create or run its workers")


class _ParallelFour(torch.autograd.Function):
    @staticmethod
    def forward(ctx, metadata, left, right, grid, eta, width, minimum, radius, epsilon, workers, *inputs):
        ctx.save_for_backward(metadata, left, right, grid, *inputs)
        ctx.settings = (eta, width, minimum, radius, epsilon, workers)
        output = grid.new_zeros((metadata.shape[0], grid.numel()))
        _call(metadata, left, right, grid, inputs, ctx.settings, output, None, None)
        return output

    @staticmethod
    def backward(ctx, upstream):
        if torch.is_grad_enabled():
            raise RuntimeError(
                "Parallel four-body projection does not support higher-order derivatives; "
                "omit parallel_workers to use the PyTorch implementation"
            )
        metadata, left, right, grid, *inputs = ctx.saved_tensors
        gradients = [torch.zeros_like(value) for value in inputs]
        _call(metadata, left, right, grid, inputs, ctx.settings, None, gradients, upstream.contiguous())
        return (None,) * 10 + tuple(gradients)


def project(environments, grid, sigma, width, minimum, tolerance, n_elements, workers, raw):
    if not environments:
        return raw
    offsets = [0, 0, 0, 0]
    metadata, lefts, rights, channels = [], [], [], []
    fields = [[] for _ in range(11)]
    for (alpha, pair_normals, pair_norms, fourth_normals, fourth_norms, distances,
         fourth_distances, win4, center_probability, forward_weights, reverse_weights,
         denominator, left, right) in environments:
        n_pair, n_near, n_fourth = left.numel(), distances.numel(), fourth_distances.numel()
        metadata.append((*offsets, n_pair, n_near, n_fourth, alpha if n_elements > 1 else -1))
        offsets[0] += n_pair
        offsets[1] += n_near * n_fourth
        offsets[2] += n_near
        offsets[3] += n_fourth
        lefts.append(left)
        rights.append(right)
        channels.append(alpha)
        values = (
            pair_normals.reshape(-1), pair_norms, fourth_normals.reshape(-1),
            fourth_norms.reshape(-1), distances, fourth_distances, win4,
            center_probability.reshape(1), forward_weights, reverse_weights,
            denominator.reshape(1),
        )
        for field, value in zip(fields, values):
            field.append(value)
    eta = 1.0 / (2.0 * sigma * sigma)
    radius = math.ceil(math.sqrt(-math.log(tolerance) / eta) / width)
    epsilon = max(torch.finfo(grid.dtype).eps * 10.0, 1e-12)
    outputs = _ParallelFour.apply(
        torch.tensor(metadata, dtype=torch.int64), torch.cat(lefts), torch.cat(rights),
        grid, eta, width, minimum, radius, epsilon, workers,
        *(torch.cat(field) for field in fields),
    )
    return raw.index_add(1, torch.tensor(channels), outputs.T)
