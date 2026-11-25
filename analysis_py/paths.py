"""Where this pipeline reads and writes.

Every path is resolved from this file's own location, so the repository runs
wherever it is cloned. Set AD_AUDIT_ROOT to point the pipeline at a checkout
other than the one this file lives in.
"""
import os
import pathlib

ROOT = pathlib.Path(os.environ.get("AD_AUDIT_ROOT",
                                   pathlib.Path(__file__).resolve().parent.parent))
HERE = ROOT / "analysis_py"

DATA = str(HERE / "data") + os.sep      # tables this pipeline writes
INPUT = str(HERE / "inputs") + os.sep   # donor-level tables it reads
FIG = str(HERE / "figures") + os.sep    # figures it writes

for _d in (DATA, INPUT, FIG):
    os.makedirs(_d, exist_ok=True)
