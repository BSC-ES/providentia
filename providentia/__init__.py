__date__ = "2026-07-10"
__version__ = "3.1.0"

def __getattr__(name):
    if name == "Providentia":
        from .library import Providentia
        return Providentia
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")