"""
Stage 4 - DOCK a candidate against the receptor.

Two engines, one interface:
  * "vina"   - AutoDock Vina (CPU). Runs anywhere, incl. a laptop. No GPU, no maps.
  * "adgpu"  - AutoDock-GPU. Fast, needs an NVIDIA GPU + a built binary + grid maps.
  * "auto"   - use adgpu if a GPU + binary are found, else vina.

dock(receptor_pdbqt, ligand_pdbqt, center, size, engine="auto") -> result dict:
  {engine, best_score (kcal/mol), poses (file), scores (list)}
"""
import os
import shutil
import subprocess
import numpy as np

DEFAULT_LIG_TYPES = "C A N NA OA SA HD S F Cl Br I P".split()


def _run(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def gpu_available():
    if not shutil.which("nvidia-smi"):
        return False
    return _run("nvidia-smi -L").returncode == 0


def find_adgpu():
    """Locate an AutoDock-GPU binary via $ADGPU or PATH."""
    env = os.environ.get("ADGPU")
    if env and os.path.isfile(env):
        return env
    for name in ("autodock_gpu_128wi", "autodock_gpu_64wi", "autodock_gpu"):
        p = shutil.which(name)
        if p:
            return p
    return None


def detect_engine(verbose=True):
    if gpu_available() and find_adgpu():
        eng = "adgpu"
    else:
        eng = "vina"
    if verbose:
        print(f"  engine=auto -> {eng}"
              + ("" if eng == "adgpu" else "  (no GPU/AutoDock-GPU found; using Vina)"))
    return eng


# ---------------------------------------------------------------- Vina engine
def _dock_vina(receptor_pdbqt, ligand_pdbqt, center, size, out,
               exhaustiveness=8, n_poses=10, seed=42, verbose=True):
    from vina import Vina
    v = Vina(sf_name="vina", seed=seed, verbosity=0)
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=list(map(float, center)), box_size=list(map(float, size)))
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    poses = f"{out}_vina.pdbqt"
    v.write_poses(poses, n_poses=n_poses, overwrite=True)
    scores = [float(e[0]) for e in v.energies(n_poses=n_poses)]
    best = min(scores) if scores else None
    if verbose:
        print(f"  [vina] best affinity = {best:.2f} kcal/mol  ({len(scores)} poses)")
    return {"engine": "vina", "best_score": best, "poses": poses, "scores": scores}


# ----------------------------------------------------------- AutoDock-GPU engine
def _write_gpf(base, receptor, ligtypes, center, size, spacing):
    rectypes = sorted({l[77:79].strip() for l in open(receptor)
                       if l[:6].strip() in ("ATOM", "HETATM")})
    npts = np.clip((np.array(size) / spacing).round().astype(int), 2, 126)
    npts = (npts + (npts % 2)).tolist()
    g = [f"npts {npts[0]} {npts[1]} {npts[2]}", f"gridfld {base}.maps.fld",
         f"spacing {spacing}", "receptor_types " + " ".join(rectypes),
         "ligand_types " + " ".join(ligtypes), f"receptor {receptor}",
         f"gridcenter {center[0]:.3f} {center[1]:.3f} {center[2]:.3f}", "smooth 0.5"]
    g += [f"map {base}.{t}.map" for t in ligtypes]
    g += [f"elecmap {base}.e.map", f"dsolvmap {base}.d.map", "dielectric -0.1465"]
    open(f"{base}.gpf", "w").write("\n".join(g) + "\n")
    return f"{base}.gpf"


def build_maps(receptor_pdbqt, center, size, base="maps", spacing=0.375,
               ligand_types=None, verbose=True):
    """Run autogrid4 to build grid maps (for the AutoDock-GPU engine)."""
    ligtypes = ligand_types or DEFAULT_LIG_TYPES
    gpf = _write_gpf(base, receptor_pdbqt, ligtypes, center, size, spacing)
    r = _run(f"autogrid4 -p {gpf} -l {base}.glg")
    fld = f"{base}.maps.fld"
    if not os.path.isfile(fld):
        raise RuntimeError(f"autogrid4 failed:\n{r.stdout[-400:]}\n{r.stderr[-400:]}")
    if verbose:
        print(f"  maps built -> {fld}")
    return fld


def _dock_adgpu(receptor_pdbqt, ligand_pdbqt, center, size, out,
                maps_fld=None, nrun=50, spacing=0.375, verbose=True):
    binary = find_adgpu()
    if not binary:
        raise RuntimeError("no AutoDock-GPU binary ($ADGPU or on PATH)")
    if maps_fld is None:
        maps_fld = build_maps(receptor_pdbqt, center, size, base=out + "_maps",
                              spacing=spacing, verbose=verbose)
    _run(f"{binary} --ffile {maps_fld} --lfile {ligand_pdbqt} "
         f"--nrun {nrun} --resnam {out}_adgpu")
    dlg = f"{out}_adgpu.dlg"
    import re
    # extract docked poses into a standard multi-model PDBQT (the dlg stores them as
    # 'DOCKED: ...' lines) so viz + clustering work the same as for Vina.
    poses = f"{out}_adgpu_poses.pdbqt"
    scores = []
    with open(dlg, errors="ignore") as fh, open(poses, "w") as out_fh:
        for l in fh:
            if l.startswith("DOCKED: "):
                out_fh.write(l[8:])
            if "Estimated Free Energy of Binding" in l:
                m = re.search(r"=\s*([-\d.]+)\s*kcal/mol", l)
                if m:
                    scores.append(float(m.group(1)))
    best = min(scores) if scores else None
    if verbose:
        print(f"  [adgpu] best dG = {best:.2f} kcal/mol")
    return {"engine": "adgpu", "best_score": best, "poses": poses, "dlg": dlg,
            "scores": scores, "maps_fld": maps_fld}


# ---------------------------------------------------------------- public API
def dock(receptor_pdbqt, ligand_pdbqt, center, size, engine="auto", out="dock",
         exhaustiveness=8, n_poses=10, maps_fld=None, nrun=50, spacing=0.375,
         seed=42, verbose=True):
    if engine == "auto":
        engine = detect_engine(verbose=verbose)
    if engine == "vina":
        return _dock_vina(receptor_pdbqt, ligand_pdbqt, center, size, out,
                          exhaustiveness, n_poses, seed, verbose)
    if engine == "adgpu":
        return _dock_adgpu(receptor_pdbqt, ligand_pdbqt, center, size, out,
                           maps_fld, nrun, spacing, verbose)
    raise ValueError(f"unknown engine '{engine}'")
