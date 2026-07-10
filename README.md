# Crop Rotations, Transition Costs, and Farm Delinquencies

This repository holds the LaTeX source for the manuscript, conference poster, and presentation materials for the project **"Crop Rotations, Transition Costs, and Farm Delinquencies."** The project studies how a corn–soybean farm's *rotation history* shapes its risk of operating-loan delinquency, using a field-level simulation that couples rotation-conditioned yield distributions with a dynamic operating-credit model.

**Authors:** Lawson Connor, Victor Funes-Leal, Eunchun Park
**Affiliation:** Department of Agricultural Economics and Agribusiness, University of Arkansas

## Overview

Operating credit is the primary financial lifeline for row-crop farmers, yet standard credit-risk models treat farm revenue as a single aggregate stream and ignore the agronomic channel through which rotation history affects subsequent-year yields. Delinquency is a *path-dependent* event: a bad year rolls unpaid debt forward at compounding interest, raising the principal financed the next year and the probability of cascading default even when later harvests are only average.

The project pursues three objectives:

1. Assess whether adopting a more complex rotation pattern improves financial outcomes (delinquency probability, operating debt, present value of net cash flow), or whether the added agronomic complexity yields diminishing financial returns relative to a simpler alternating strategy.
2. Estimate the probability that a farm incurs at least one operating-loan delinquency over a six-year horizon, conditioning yield distributions on the field's full six-year rotation history.
3. Compare delinquency probabilities across rotation strategies while accounting for the negative correlation between yields and prices.

## Method summary

Yield distributions are estimated for each unique six-year rotation state, and a transition-path model tracks how the conditioning state evolves as a planned rotation unfolds from a given history. These rotation-conditioned yield generators feed a dynamic operating-credit model that simulates the full annual cycle: debt origination at planting, nine-month interest accrual, harvest-time repayment from realized revenue, and balance rollover when revenue is insufficient to retire the obligation. The negative yield–price correlation is estimated from the panel and embedded in the simulation via copula methods.

## Data sources

- **Yields:** field-level corn and soybean yields estimated with the QDANN yield-mapping model.
- **Rotation histories:** six-year crop sequences for Illinois fields, 2008–2022, derived from the USDA Cropland Data Layer (CDL).
- **Prices and costs:** grain elevator prices used to calibrate marginal distributions and the yield–price correlation.
- **Weather:** GRIDMET.
- **Vapor pressure deficit and Palmer Drought Severity Index (PDSI):** TerraClimate.
- **Soil characteristics:** gSSURGO.

## Repository structure

| Path | Description |
|------|-------------|
| `crop_rotation_delinquencies.tex` | Main manuscript source. |
| `abstract.tex` | Standalone abstract. |
| `AAEA_2026_poster.tex` | Conference poster (AAEA 2026), built with `beamerposter`. |
| `waea_presentation.tex` | WAEA presentation slides. |
| `framework_beamer.tex` | Beamer slides describing the modeling framework. |
| `references.bib` | BibTeX bibliography. |
| `beamerposter.sty` | Poster layout package. |
| `beamerthemeconfposter.sty` | Poster theme. |
| `beamerthemeconfposteralt.sty` | Alternate poster theme. |
| `figs/` | Figures used across the manuscript and presentations. |

## Building the documents

Each `.tex` file is a standalone document. A typical build for the manuscript:

```bash
pdflatex crop_rotation_delinquencies.tex
bibtex   crop_rotation_delinquencies
pdflatex crop_rotation_delinquencies.tex
pdflatex crop_rotation_delinquencies.tex
```

The poster and presentation files (`AAEA_2026_poster.tex`, `waea_presentation.tex`, `framework_beamer.tex`) compile the same way. The poster depends on the bundled `beamerposter.sty` and `beamerthemeconfposter*.sty` files, so compile it from the repository root. Any recent TeX Live or MiKTeX distribution should provide the remaining dependencies.

## Citation

If you reference this work, please cite it as:

> Connor, L., Funes-Leal, V., and Park, E. *Crop Rotations, Transition Costs, and Farm Delinquencies.* Working paper.

## License

No formal license is currently attached to this repository, which by default means all rights are reserved and the material may not be reused without permission. If you intend to make the work reusable, add a `LICENSE` file and update this section. For an academic manuscript and its figures, [Creative Commons Attribution 4.0 (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) is a common choice, allowing others to share and adapt the material with attribution. Note that the bundled poster template files (`beamerposter.sty`, `beamerthemeconfposter*.sty`) carry their own upstream licenses, which take precedence for those files.

## Contact

For questions about the project, contact the authors at the Department of Agricultural Economics and Agribusiness, University of Arkansas.
