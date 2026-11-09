# Four-wave mixing and coherent couplers

Working repository for a combined paper on elliptic complex-envelope solutions, normal forms and the projection from degenerate four-wave mixing to coherent couplers.

## Layout

- `latex/`: the new manuscript, combined bibliography and figure assets.
- `code/`: copied symbolic notebooks, shared elliptic-function helpers and numerical audits.
- `fwm-paper/` and `coherent-coupler-paper/`: local, Git-ignored source archives. Leave these unchanged.

The manuscript entry point is `latex/combined-paper.tex`. The first combined draft retains the FWM elliptic derivation and determinant treatment, states the asymmetric coupler solution, develops its reciprocal gauge and degenerate-FWM projection, and gives a shared Kronecker/Baker–Akhiezer discussion. Numerical evaluation follows the theory. Build it to produce `latex/combined-paper.pdf`.

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
MPLCONFIGDIR=/tmp/fwm-coupler-matplotlib .venv/bin/python code/audit_numeric_notebook.py
```

The notebook audit executes selected cells in memory without rewriting the notebook. Interactive use of the notebooks additionally needs IPython and a Jupyter front end in the chosen environment.

## LaTeX

Build from `latex/` with a TeX installation containing `latexmk`, `pdflatex` and `biber`:

```sh
cd latex
latexmk -pdf -interaction=nonstopmode -halt-on-error combined-paper.tex
```

See `latex/README.md` for figure paths and the current manuscript status.
