"""Split a node core/memory budget across concurrent Gaussian jobs."""

from __future__ import annotations


_MIN_GAUSSIAN_JOB_MEM_GB = 4


def split_nproc_for_items(
    nproc: int,
    n_items: int,
    *,
    prefer_depth: bool = False,
    min_inner: int | None = None,
    flatten_nested: bool = True,
) -> tuple[int, int]:
    """Split ``nproc`` into ``(n_outer_workers, n_inner_per_worker)``.

    When ``prefer_depth`` is True, keep at least ``min_inner`` cores per outer
    worker (default 1). When ``flatten_nested`` is True, never return both
    outer and inner greater than 1.
    """
    nproc = max(1, int(nproc))
    n_items = max(1, int(n_items))
    if n_items == 1:
        return 1, nproc
    if not prefer_depth:
        n_outer = min(nproc, n_items)
        n_inner = max(1, nproc // n_outer)
    else:
        if min_inner is None:
            min_inner = 1
        min_inner = max(1, int(min_inner))
        max_outer = min(n_items, max(1, nproc // min_inner))
        n_outer, n_inner = 1, nproc
        best_used, best_outer = nproc, 1
        for cand_outer in range(1, max_outer + 1):
            cand_inner = nproc // cand_outer
            if cand_inner < min_inner:
                continue
            used = cand_outer * cand_inner
            if used > best_used or (used == best_used and cand_outer > best_outer):
                n_outer, n_inner = cand_outer, cand_inner
                best_used, best_outer = used, cand_outer
    if flatten_nested and n_outer > 1 and n_inner > 1:
        if prefer_depth:
            return 1, nproc
        return min(nproc, n_items), 1
    return n_outer, n_inner


def split_gaussian_orientation_budget(
    nproc: int,
    n_jobs: int,
    mem_gb: int,
    *,
    min_mem_gb: int = _MIN_GAUSSIAN_JOB_MEM_GB,
) -> tuple[int, int, int]:
    """Split cores and GB across concurrent Gaussian orientation ESP jobs.

    Returns ``(n_workers, job_nproc, job_mem_gb)`` with
    ``n_workers * job_nproc <= nproc`` and
    ``n_workers * job_mem_gb <= mem_gb`` (each at least 1).
    """
    nproc = max(1, int(nproc))
    n_jobs = max(1, int(n_jobs))
    mem_gb = max(1, int(mem_gb))
    min_mem_gb = max(1, int(min_mem_gb))
    floor = min(min_mem_gb, mem_gb)
    max_by_mem = max(1, mem_gb // floor)
    n_cap = min(n_jobs, nproc, max_by_mem)
    n_workers, job_nproc = split_nproc_for_items(
        nproc,
        n_cap,
        prefer_depth=True,
        min_inner=1,
        flatten_nested=False,
    )
    n_workers = max(1, int(n_workers))
    job_nproc = max(1, int(job_nproc))
    job_mem = max(1, mem_gb // n_workers)
    return n_workers, job_nproc, job_mem
