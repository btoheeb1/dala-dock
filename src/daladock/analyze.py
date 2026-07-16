"""
Stage 5 - ANALYZE results -> ranked leaderboard.

leaderboard(results) takes a dict {name: dock_result} (the dicts returned by
dock()) and returns a pandas DataFrame sorted by best score (most negative first).
"""
import numpy as np
import pandas as pd


import re

_COLS = ["pose", "affinity_kcal_mol", "rmsd_to_best_A", "cluster", "cluster_size"]


def _parse_poses(poses_file):
    """
    Parse a multi-model pose file - works for Vina PDBQT (MODEL/ATOM +
    'REMARK VINA RESULT') and AutoDock poses ('USER ... Estimated Free Energy').
    Returns (energies, coords): one energy (may be None) and one heavy-atom
    coordinate array per pose.
    """
    energies, coords, cur, en = [], [], [], None
    for l in open(poses_file, errors="ignore"):
        tag = l[:6]
        if tag.startswith("MODEL"):
            cur, en = [], None
        elif tag.startswith("ENDMDL"):
            if cur:
                energies.append(en); coords.append(np.array(cur)); cur = []
        elif tag.strip() in ("ATOM", "HETATM"):
            if len(l) >= 79 and l[77:79].strip() not in ("H", "HD"):   # heavy atoms
                try:
                    cur.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                except ValueError:
                    pass
        elif "VINA RESULT" in l:
            try: en = float(l.split()[3])
            except (IndexError, ValueError): pass
        elif "Estimated Free Energy of Binding" in l:
            m = re.search(r"=\s*([-\d.]+)", l)
            if m: en = float(m.group(1))
    if cur:
        energies.append(en); coords.append(np.array(cur))
    return energies, coords


def _rmsd(a, b):
    """In-place heavy-atom RMSD (no superposition) - the docking convention."""
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def cluster_poses(poses_file, affinities=None, cutoff=2.0):
    """
    Cluster docked poses the way AutoDock does: order poses best-energy first,
    put each into the first cluster whose representative is within `cutoff` A
    (RMSD), else start a new cluster. Works for Vina and AutoDock-GPU output.
    Returns a DataFrame (one row per pose, ranked best-first).
    """
    energies, coords = _parse_poses(poses_file)
    n = len(coords)
    if n == 0:
        return pd.DataFrame(columns=_COLS)
    if affinities:                                   # backfill any missing energies
        for i in range(n):
            if energies[i] is None and i < len(affinities):
                energies[i] = affinities[i]
    order = sorted(range(n), key=lambda i: (energies[i] is None, energies[i]))
    best = order[0]
    reps, assign = [], {}
    for i in order:                                  # walk best-energy first
        for ci, rep in enumerate(reps):
            if _rmsd(coords[i], coords[rep]) <= cutoff:
                assign[i] = ci; break
        else:
            assign[i] = len(reps); reps.append(i)
    sizes = {}
    for i in order:
        sizes[assign[i]] = sizes.get(assign[i], 0) + 1
    rows = [{
        "pose": rank,
        "affinity_kcal_mol": round(energies[i], 2) if energies[i] is not None else None,
        "rmsd_to_best_A": round(_rmsd(coords[i], coords[best]), 2),
        "cluster": assign[i] + 1,
        "cluster_size": sizes[assign[i]],
    } for rank, i in enumerate(order, 1)]
    return pd.DataFrame(rows)


def leaderboard(results, csv=None):
    rows = []
    for name, r in results.items():
        scores = r.get("scores") or []
        rows.append({
            "ligand": name,
            "best_score_kcal_mol": r.get("best_score"),
            "engine": r.get("engine"),
            "n_poses": len(scores),
            "mean_top3": (round(sum(sorted(scores)[:3]) / min(3, len(scores)), 2)
                          if scores else None),
        })
    df = pd.DataFrame(rows).sort_values(
        "best_score_kcal_mol", na_position="last").reset_index(drop=True)
    df.index += 1
    df.index.name = "rank"
    if csv:
        df.to_csv(csv)
        print(f"  wrote {csv}")
    return df
