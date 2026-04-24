"""Read ``fragmentomYcs.cfg`` for pipeline configuration.

The file is searched in order: ``FRAGMENTOMYCS_CONFIG`` env-var directory →
cwd → parent directories → module directory → project root.
"""

from __future__ import annotations

import configparser
import os
from pathlib import Path
from typing import Optional

_CONFIG_FILENAME = "fragmentomYcs.cfg"
_SECTION = "bcftools"
_KEY_BACKEND = "backend"

# Cached parsed config — loaded once per interpreter session.
_config: Optional[configparser.ConfigParser] = None


# ---------------------------------------------------------------------------
# Config file discovery
# ---------------------------------------------------------------------------

def _find_config_file() -> Optional[Path]:
    """Return the first ``fragmentomYcs.cfg`` found, or ``None``."""
    def _candidates():
        env_dir = os.environ.get("FRAGMENTOMYCS_CONFIG", "")
        if env_dir:
            yield Path(env_dir) / _CONFIG_FILENAME
        yield Path.cwd() / _CONFIG_FILENAME
        yield from (p / _CONFIG_FILENAME for p in Path.cwd().parents)
        module_dir = Path(__file__).parent
        yield module_dir / _CONFIG_FILENAME
        yield module_dir.parent / _CONFIG_FILENAME

    return next((p for p in _candidates() if p.is_file()), None)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------

def load_config(force_reload: bool = False) -> configparser.ConfigParser:
    """Return the parsed ``ConfigParser`` instance.

    The result is cached after the first successful load.  Pass
    ``force_reload=True`` to discard the cache (useful in tests).
    """
    global _config
    if _config is not None and not force_reload:
        return _config

    cfg = configparser.ConfigParser()
    cfg_path = _find_config_file()
    if cfg_path is not None:
        cfg.read(cfg_path)

    _config = cfg
    return _config


def get_bcftools_backend() -> str:
    """Return the configured normalisation backend.

    Reads the ``[bcftools] backend`` key from ``fragmentomYcs.cfg``.
    Returns ``"auto"`` when the file or key is absent.

    Valid values (case-insensitive): ``auto``, ``pybcftools``, ``bcftools``.
    """
    cfg = load_config()
    if cfg.has_option(_SECTION, _KEY_BACKEND):
        return cfg.get(_SECTION, _KEY_BACKEND).strip().lower()
    return "auto"
