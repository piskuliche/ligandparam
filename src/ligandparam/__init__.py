from importlib.metadata import PackageNotFoundError, version as _version

__logging_name__ = "ligandparam"

try:
    # Single source of truth: the version declared in pyproject.toml, read back
    # from the installed distribution metadata. Keeping a literal here caused it
    # to drift (it read "0.3.0" for every release up to and including 1.0.0).
    __version__ = _version("ligandparam")
except PackageNotFoundError:  # running from a source tree that was never installed
    __version__ = "0.0.0+unknown"
