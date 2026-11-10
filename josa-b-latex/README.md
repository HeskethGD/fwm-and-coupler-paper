# FWM and coherent-coupler manuscript

This is the sole working manuscript directory. The entry point is `fwm-and-coupler.tex`, accompanied by `supplement.tex`.

## Build

From the repository root:

```sh
make -C josa-b-latex
```

The outputs are `fwm-and-coupler.pdf` and `supplement.pdf`. Keep them together so links between the documents work. The Makefile builds them sequentially, twice, to resolve their mutual references from a clean checkout. Requires Python 3, `latexmk`, `pdflatex`, BibTeX, biber and the packages loaded by the publisher class. Do not build the two documents concurrently while their auxiliary files are being written.

## Manuscript structure

The main paper presents three results: complex-envelope reconstruction and its Baker–Akhiezer structure; a form-preserving FWM normalisation; and the degenerate-FWM projection onto a gauged asymmetric coherent coupler. The BA identification closes the general-solution section and motivates the two transformations. Assumptions appear with the equations and maps that first use them.

The supplement contains the FS determinant calculations, expanded derivations, signed elliptic-point data, spectral identities, Kronecker/theta representation, full parameter data and additional numerical checks. The main paper's bibliography uses the publisher's native `opticajnl.bst`; the supplement uses its own biblatex database `supplement-references.bib`. The required derivations are included with the submission.

The author is identified as an independent researcher based in the United Kingdom. The Baker–Akhiezer function is denoted by `\Psi`, separate from the phase symbols `\phi`. Other mathematical symbols retain their existing system-specific scopes.

## Publisher format and page estimate

`optica-article.cls` is the official **single-column universal submission format**, obtained from the author's downloaded Optica package. The class and publisher style assets are copied unchanged; `template-provenance.json` records their hashes. The source is [Optica Publishing Group's template page](https://opg.optica.org/content/author/portal/item/templates-default).

The publisher template estimates final two-column JOSA B length by multiplying universal-template pages, measured in quarter-page increments, by **60%**. Equation displays, figure sizing and production decisions affect the result. This is a planning estimate, not a production page count.

The current build has **16 main-paper pages** (including references) and **33 supplemental pages**. Applying the 60% rule to all 16 main pages gives an approximate planning estimate of **9.6 published pages**; the final reference page is partially filled.

## Figures and calculations

- `plots/case1_plot_1.pdf`: the general complex FWM example in the supplement.
- `plots/case2_plot_1.pdf`: the physical FWM example in the main paper.
- `plots/coupler/asymmetric_physical.pdf`: the asymmetric physical-coupler example.

Additional original figure assets are retained under `plots/coupler/`. See `../code/README.md` for the generation commands. The calculation runners export directly to this directory; all 11 historical notebooks remain unchanged.

Generated manuscript PDFs and auxiliary files are Git-ignored. Template files, manuscript source, bibliography databases and figures are included so this directory builds independently of the archived original papers.

## Cross-document references in submission packages

The documents load `main-paper-references.tex` or `supplement-references.tex` to resolve references to the other document. The generated reference data must accompany a standalone build: submission platforms may discard auxiliary files. For the submission ZIPs in `dist/`, embed the contents of the relevant reference input directly in the packaged manuscript, replacing its input command. This leaves exactly one `.tex` file per ZIP, as required by the platform. The Makefile refreshes them from the latest builds using `export-references.py`. Each separate package builds without the other document or any pre-existing `.aux` files.

Cross-document references are deliberately printed without hyperlinks (`\ref*` and `\externaleqref`), because submission-platform PDF handling can redirect external links to the wrong document. Internal references retain their hyperlinks.
