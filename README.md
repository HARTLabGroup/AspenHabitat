---
editor_options: 
  markdown: 
    wrap: 72
---

## Description of the project

This repository contains all code for *Climate-driven decreases in
aspen's distribution and opportunities for future expansion across the
Southern Rocky Mountains*.

## Organization of the project

The project has the following structure:

-   [Code](https://github.com/HARTLabGroup/AspenHabitat/tree/base/Code):
    All code written for this project

    -   [AspenHabitatModelling.docx](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Code/AspenHabitatModelling.docx):
        A Microsoft Word version of the the manuscript generated from
        AspenHabitatModelling.Rmd

    -   [AspenHabitatModelling.Rmd](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Code/AspenHabitatModelling.Rmd):
        An R Markdown file containing code and text for generating
        AspenHabitatModelling.docx

    -   [GIDS-Downscaling-SJH.R](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Code/GIDS-Downscaling-SJH.R):
        An R file that performs a multi-step procedure that spatially
        downscales gridded climate data using GIDS (Gradient and Inverse
        Distance-Squared) of Nalder and Wein (1998) and Flint and Flint
        (2012). It was lightly modified from code written by Rodman et
        al. (2020).

    -   [Template.docx](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Code/Template.docx):
        A Microsoft Word document used to a style template for
        generating AspenHabitatModelling.docx

    -   [references.bib](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Code/references.bib):
        A BibTeX file containing all information for all references used
        in the project

-   [Documents](https://github.com/HARTLabGroup/AspenHabitat/tree/base/Documents):
    All code written for this project

    -   [ODMAP-Table.xlsx](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Documents/ODMAP-Table.xlsx):
        Documentation of the species distribution modeling methods used
        herein following standard protocols outlined by Zurell et al.
        (2020)

    -   [PredictorScreeningTable.xlsx](https://github.com/HARTLabGroup/AspenHabitat/blob/base/Documents/PredictorScreeningTable.xlsx){.uri}:
        A table containing predictors screened for inclusion in the
        models.

-   Figures: Any figures that were not generated from the provided code

    -   StudyArea.jpg: A map of the study area generated in
        [QGIS](https://www.qgis.org/).

-   .gitignore

-   AspenHabitat.Rproj

-   CITATION.cff

-   License.md

-   README.md

-   citationstyle.csl : A [Citation Style
    Language](https://citationstyles.org/) file used to format
    references in the main text of AspenHabitatModelling.docx\

## License

This project is licensed under the MIT License - see the
[LICENSE.md](https://github.com/HARTLabGroup/AspenHabitat/blob/base/License.md)
file for details

## Citation

Please cite this work following the
[Citation.cff](https://github.com/HARTLabGroup/AspenHabitat/blob/base/CITATION.cff)
file

## References

Flint, L. E., and A. L. Flint. 2012. Downscaling future climate
scenarios to fine scales for hydrologic and ecological modeling and
analysis. Ecological Processes 1:2.

Nalder, I. A., and R. W. Wein. 1998. Spatial interpolation of climatic
Normals: test of a new method in the Canadian boreal forest.
Agricultural and Forest Meteorology 92:211–225.

Rodman, K., T. Veblen, M. Battaglia, M. Chambers, P. Fornwalt, Z.
Holden, T. Kolb, J. Ouzts, and M. Rother. 2020, August 18. Data from: A
changing climate is snuffing out post-fire recovery in montane forests.
Dryad.
[*https://doi.org/10.5061/DRYAD.QZ612JMB7*](https://doi.org/10.5061/DRYAD.QZ612JMB7)

Zurell, D., Franklin, J., König, C., Bouchet, P.J., Dormann, C.F.,
Elith, J., Fandos, G., Feng, X., Guillera-Arroita, G., Guisan, A.,
Lahoz-Monfort, J.J., Leitão, P.J., Park, D.S., Peterson, A.T.,
Rapacciuolo, G., Schmatz, D.R., Schröder, B., Serra-Diaz, J.M.,
Thuiller, W., Yates, K.L., Zimmermann, N.E., Merow, C., 2020. A standard
protocol for reporting species distribution models. Ecography 43,
1261–1277.
[*https://doi.org/10.1111/ecog.04960*](https://doi.org/10.1111/ecog.04960)
