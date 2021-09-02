
from shutil import which


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    # from whichcraft import which

    return which(name) is not None
