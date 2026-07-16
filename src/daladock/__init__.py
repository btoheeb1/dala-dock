"""
dala-dock: a small, portable molecular-docking teaching pipeline.

Five stages, one per module, mirroring the command-line pipeline:
  1. prep_receptor  - target (PDB/PDB-id/SMILES) -> receptor.pdbqt
  2. box            - define the search box (reference ligand / fpocket / blind)
  3. prep_ligand    - candidate (SMILES/file) -> ligand.pdbqt
  4. dock           - Vina (CPU, anywhere) or AutoDock-GPU (if a GPU is available)
  5. analyze        - collect results into a ranked leaderboard

`viz` provides py3Dmol helpers to see the receptor, the box, and docked poses.
"""
import os as _os
import sys as _sys

# Sanitize the environment on import so subprocess tools (obabel, autogrid4, ...)
# resolve to THIS environment, not a leaked one (e.g. via a stray PYTHONPATH or
# another conda env earlier on PATH). Fixes the classic "wrong meeko/no gemmi" clash.
_os.environ.pop("PYTHONPATH", None)
_envbin = _os.path.dirname(_sys.executable)
if _envbin and _envbin not in _os.environ.get("PATH", "").split(_os.pathsep):
    _os.environ["PATH"] = _envbin + _os.pathsep + _os.environ.get("PATH", "")

from . import prep_receptor, box, prep_ligand, dock, analyze, viz  # noqa: F401,E402

__all__ = ["prep_receptor", "box", "prep_ligand", "dock", "analyze", "viz"]
__version__ = "0.1.0"
