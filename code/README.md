# Calculations and validation

This directory combines the existing calculation files from the two papers. The archived source directories are not runtime dependencies and must not be edited.

## Preserved notebooks

All 11 notebooks are exact copies, including their saved outputs, metadata and manual calculations. They have not been re-executed or reformatted during the copy.

| Source | Notebook |
| --- | --- |
| FWM code | `Four Wave Mixing Case.ipynb` |
| FWM code | `The canonical coordinates of FWM.ipynb` |
| FWM code | `Numeric plots for four-wave mixing.ipynb` |
| FWM notes | `Mobius transforms and invariance.ipynb` |
| Coupler code | `The Coherent Coupler system of Jensen.ipynb` |
| Coupler code | `The general 2 mode case of the coherent coupler.ipynb` |
| Coupler code | `Numeric plots for the coherent coupler.ipynb` |
| Coupler code | `Degenerate FWM to coherent coupler transform check.ipynb` |
| Coupler code | `The projection from degenerate FWM to the coherent coupler.ipynb` |
| Coupler code | `The connection between a 2-mode coupler and a 2-mode multimode system.ipynb` |
| Coupler code | `The connection between a 2-mode coupler and a polarisation system.ipynb` |

`source-manifest.json` records original paths, destination paths and SHA-256 hashes of the copied notebooks and code. Shared `weierstrass_modified.py`, `wpk.py`, `__init__.py` and `requirements.txt` were identical between the two papers and are included once. The manifest contains both source records for shared files.

## Minimal adaptations

The copied mathematical helpers and audit scripts are unchanged. The only Python source adaptation at setup is in `asymmetric_physical_example.py`: its figure destination is now `../latex/plots/coupler/asymmetric_physical.pdf`, resolved relative to the script. No equations or numerical parameters have changed. Its original-copy hash is retained in the manifest as provenance, not as a claim that the adapted script remains identical.

Run commands from the repository root using `.venv/bin/python`; examples are in the root README. The old audit docstrings mentioning `code/.venv/bin/python` describe their historical environment location.

To regenerate figures without changing notebook files:

```sh
MPLCONFIGDIR=/tmp/fwm-coupler-matplotlib .venv/bin/python code/audit_numeric_notebook.py --figures
MPLCONFIGDIR=/tmp/fwm-coupler-matplotlib .venv/bin/python code/asymmetric_physical_example.py
```

The FWM runner writes its original filenames to `latex/plots/`. The coupler physical example writes to `latex/plots/coupler/`. Importing `audit_coupler` in the physical-example script also executes its existing checks; that behaviour is inherited unchanged.

## Historical notebook paths

For interactive work, launch notebooks with this directory as the working directory so the local helper imports resolve.

The FWM numerical notebook's active `../latex/plots/` export paths already match the new layout. The coupler numerical notebook retains its original `../paper/plots/` export paths. Those paths do not match the new layout; running its export cells unchanged will fail unless the old destination exists. Preserve the notebook record and use the audited physical-example script for the retained coupler figure. If additional coupler plots are needed later, redirect those export paths in memory in an audit runner rather than rewriting the notebook.

The exploratory `Mobius transforms and invariance.ipynb` also retains a historical `numerical_evaluation.weierstrass_modified` import. That package is not provided here; its direct counterpart is the unchanged local `weierstrass_modified.py`. No import shim or notebook edit has been introduced. This exploratory notebook is preserved as calculation history and is not part of the current audited execution path.

Copy integrity does not imply every historical notebook runs from a fresh kernel. Use the audit scripts for reproducible checks and discuss any necessary notebook changes before making them.
