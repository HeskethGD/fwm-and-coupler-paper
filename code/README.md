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

The mathematical helpers and original notebooks are unchanged. Two scripts adapt figure destinations to the current manuscript directory:

- `asymmetric_physical_example.py` exports to `../josa-b-latex/plots/coupler/asymmetric_physical.pdf`, resolved relative to the script.
- `audit_numeric_notebook.py --figures` executes the notebook's figure cells with export paths redirected in memory to `../josa-b-latex/plots/`.

These adaptations change output paths only. `source-manifest.json` retains the original-copy hashes as provenance; it is not a claim that the adapted scripts remain byte-identical to their historical copies. The notebooks themselves remain byte-identical.

Run commands from the repository root using `.venv/bin/python`; examples are in the root README. The old audit docstrings mentioning `code/.venv/bin/python` describe their historical environment location.

`projection_numeric_check.py` is a new standalone check. It independently integrates the six-component degenerate-FWM system and the four-component gauged coupler, then compares the projected source fields with the target fields at 501 points over `0 <= z <= 8`. It uses the manuscript's asymmetric physical-coupler initial data, verifies initial invariant matching at 70 digits, and repeats the double-precision integration with tighter tolerances. It prints the matching constants, initial states, errors, invariant drifts and library versions as JSON, without evaluating elliptic envelopes or importing the historical audit scripts.

```sh
.venv/bin/python code/projection_numeric_check.py
```

The tight-run maximum absolute difference is approximately `1.6e-12`, compared with `6.2e-11` at looser tolerances. This is a discrepancy between two independently integrated solutions, not a bound on either integration's error. The original mathematical calculations, helpers and notebooks are retained; the path adaptations are listed above.

To regenerate figures without changing notebook files:

```sh
MPLCONFIGDIR=/tmp/fwm-coupler-matplotlib .venv/bin/python code/audit_numeric_notebook.py --figures
MPLCONFIGDIR=/tmp/fwm-coupler-matplotlib .venv/bin/python code/asymmetric_physical_example.py
```

The FWM runner writes its original figure filenames to `josa-b-latex/plots/`. The coupler physical example writes to `josa-b-latex/plots/coupler/`. Importing `audit_coupler` in the physical-example script also executes its existing checks; that behaviour is inherited unchanged.

## Historical notebook paths

For interactive work, launch notebooks with this directory as the working directory so the local helper imports resolve.

The FWM notebook retains its historical `../latex/plots/` exports and the coupler notebook retains its historical `../paper/plots/` exports. Those folders are not working destinations. Use the commands above for figure regeneration: the runner redirects FWM export paths in memory, and the physical-example script supplies the retained coupler figure. For direct interactive use, redirect export paths in memory to `../josa-b-latex/plots/` (or its `coupler/` subdirectory) before running export cells. Do not rewrite the stored notebook records to change paths.

The exploratory `Mobius transforms and invariance.ipynb` also retains a historical `numerical_evaluation.weierstrass_modified` import. That package is not provided here; its direct counterpart is the unchanged local `weierstrass_modified.py`. No import shim or notebook edit has been introduced. This exploratory notebook is preserved as calculation history and is not part of the current audited execution path.

Copy integrity does not imply every historical notebook runs from a fresh kernel. Use the audit scripts for reproducible checks and discuss any necessary notebook changes before making them.
