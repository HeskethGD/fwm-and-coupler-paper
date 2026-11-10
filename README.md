# Four-wave mixing and coherent couplers

Working repository for a combined paper on elliptic complex-envelope solutions, normal forms and the projection from degenerate four-wave mixing to coherent couplers.

## Layout

- `code/`: copied symbolic notebooks, shared elliptic-function helpers and numerical audits.
- `josa-b-latex/`: the manuscript, supplement, bibliography databases, publisher template and figure assets; the sole working LaTeX directory.

The manuscript entry point is `josa-b-latex/fwm-and-coupler.tex`, accompanied by `josa-b-latex/supplement.tex`. The paper develops the FWM complex envelopes and their Baker–Akhiezer factors, then applies that analytic structure to the FWM normal form and the coherent-coupler projection. The supplement contains the FS determinant representations, Kronecker/theta formulas, detailed spectral calculations, parameter definitions and additional numerical material. Build both PDFs together so their cross-references resolve.

All 11 original Jupyter notebooks, including the exploratory FWM notes notebook, have been copied byte-for-byte into `code/`. See `code/README.md` for provenance, historical paths and validation commands.

The repository reference for the combined work is [HeskethGD/fwm-and-coupler-paper](https://github.com/HeskethGD/fwm-and-coupler-paper).

## Python checks

Use Python 3.13 and install the existing dependencies from `requirements.txt`. The root `.venv` is the working environment; no environment from either archived paper is needed.

```sh
python3.13 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python code/audit_initial.py
.venv/bin/python code/audit_final.py
.venv/bin/python code/audit_cubic.py
.venv/bin/python code/audit_coupler.py
.venv/bin/python code/projection_numeric_check.py
MPLCONFIGDIR=/tmp/fwm-coupler-matplotlib .venv/bin/python code/audit_numeric_notebook.py
```

The notebook audit executes selected cells in memory without rewriting the notebook. Interactive use of the notebooks additionally needs IPython and a Jupyter front end in the chosen environment.

## LaTeX

Build from `josa-b-latex/` with a TeX installation containing `latexmk`, `pdflatex`, BibTeX and biber:

```sh
cd josa-b-latex
make
```

See `josa-b-latex/README.md` for figure paths and the current manuscript status.
