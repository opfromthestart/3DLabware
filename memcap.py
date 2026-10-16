"""Cap this process's memory, so a runaway boolean fails instead of freezing.

A union that goes wrong does not raise. It allocates until the kernel starts
swapping, and the whole machine stalls for a minute before the OOM killer picks
a victim -- which may not even be the process at fault. `setrlimit(RLIMIT_AS)`
turns that into a MemoryError inside the offending call: one process reports a
traceback, the rest of the machine keeps running.

The cap comes from `AIGEN_MAX_RAM_GB`, defaulting to 4 GB, so a deliberately
large build can be given more:

    AIGEN_MAX_RAM_GB=12 python3 make_column.py

Call `memcap.cap()` first thing in a script that runs booleans.
"""

import os
import resource

DEFAULT_GB = 4.0


def cap(gb=None, quiet=False):
    """Limit this process's address space. Returns the limit in bytes."""
    if gb is None:
        gb = float(os.environ.get("AIGEN_MAX_RAM_GB", DEFAULT_GB))
    if gb <= 0.0:
        return None
    limit = int(gb * 1024 ** 3)
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    if hard != resource.RLIM_INFINITY:
        limit = min(limit, hard)
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit, hard))
    except (ValueError, OSError):
        return None
    if not quiet:
        print(f"  memory cap        : {limit / 1024 ** 3:.1f} GB "
              f"(AIGEN_MAX_RAM_GB to change)")
    return limit


def used_mb():
    """Resident size in MB, for reporting how close a run came to the cap."""
    try:
        with open("/proc/self/statm") as fh:
            return int(fh.read().split()[1]) * os.sysconf("SC_PAGE_SIZE") \
                / (1024.0 * 1024.0)
    except (OSError, MemoryError):
        return -1.0          # reporting must not be what crashes


if __name__ == "__main__":
    cap(quiet=False)
    print(f"  python, before work: {used_mb():.0f} MB resident")
