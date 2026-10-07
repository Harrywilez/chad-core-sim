"""Compute the unit-scale field grids for every Chad-core configuration (cached in cache/)."""
import sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from ccsim.runner import unit_catalogue, build_grids
t0 = time.time()
cat = unit_catalogue()
grids = build_grids(cat, Path(__file__).resolve().parent.parent / "cache", n=41, with_A=True)
print(f"{len(grids)} grids ready in {time.time()-t0:.0f}s")
