# Combined manuscript workspace

`combined-paper.tex` is the entry point for the first combined draft, titled **Elliptic complex-envelope solutions and normal forms for four-wave mixing and coherent couplers**.

The draft retains the FWM physical model, elliptic derivation, complete FS treatment, cubic-case discussion and normal-form construction. It adds a compact statement of the asymmetric coupler solution and the full reciprocal-gauge and degenerate-FWM projection construction. A shared section develops both Kronecker/Baker–Akhiezer identifications, Wronskians and Lamé equations. Numerical evaluation follows the theory, with three figures in this first draft.

The coupler's repeated power/envelope derivation and FS calculation are referenced to its preprint rather than reproduced. Parameter definitions and signed elliptic-point data needed to use the stated solution remain in the combined paper. Title, abstract, introduction, contribution list and conclusion cover both systems.

The original manuscripts remain untouched in their Git-ignored archive directories. Make subsequent manuscript changes here.

`references.bib` contains the union of the two original bibliographies, retaining the FWM entry when a citation key is shared, plus the new combined-repository entry `heskethFWMAndCouplerPaper`. DOI and arXiv identifier comparisons found no same-identifier/different-key duplicates; this does not constitute an exhaustive bibliographic novelty review. Historical repository references remain available in the database; the combined code-availability statement cites the new repository. Internal cross-references replace the coupler paper's references back to the FWM preprint where that material is now in the combined draft.

## Figure assets

- `plots/case1_plot_1.pdf` and `plots/case2_plot_1.pdf`: original FWM plots. Their paths also match the preserved FWM numerical notebook.
- `plots/coupler/`: all five original coupler plots, including `asymmetric_physical.pdf`. The separate directory prevents collisions with FWM figure filenames.

All seven original figures are available. The first combined draft includes both FWM examples and the asymmetric physical coupler example. The other four coupler figures are retained as assets but are not included. The coupler physical figure has been regenerated using the copied script and its existing parameters. A separate numerical projection-comparison plot has not been added; the full ratio projection is checked symbolically by the existing audit.

## Build

From this directory:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error combined-paper.tex
```

The shared `iopjournal.cls` is unchanged and was identical in the two source papers. The source defaults to the existing preprint article format. Build products are Git-ignored.

Existing mathematical notation is preserved. The coupler and lift sections explicitly reset the scope of the original, hatted and barred variables and the chi parameter. The common tilde vector field provides the bridge, without identifying unrelated trajectories or reconstruction constants. Imported LaTeX labels use a `cc:` namespace so that equation references remain unique without altering printed mathematical symbols.

The 67 displayed mathematical environments in the original FWM body through the normal-form construction are unchanged. Most imported coupler displays are copied in full; the scalar solution and Lambda coefficient definitions are selected from their original displays to avoid repeating the derivation. The source notebooks and both archived manuscripts remain unchanged.
