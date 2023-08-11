---
title: "README file"
author: "Sarah Hart and Asha Paudel"
contact: "Sarah.Hart@colostate.edu"
date: "June 15, 2023"
output: rmarkdown::html_vignette
---

## Description of the project

This objective of this project is to model the habitat suitability of aspen across the Southern Rocky Mountains and forecast future habitat suitability under several climate change scenarios.

## Organization of the project

The project has the following structure:

```
AspenHabitat
|-- AspenHabitat.Rproj              # RStudio project
|-- CITATION.md                     # information on how to cite
|-- Code                            # All code written for this project
    |-- AspenHabitatModelling.rmd   # The code for all analyses
    |-- Back                        # A folder containing any old code that is no longer used
    |-- BMV1_Code.R                 # Code for downscaling climate data modified from Rodman et al. 2020 (doi:10.1111/geb.13174)
    |-- GIDS-Downscaling-SJH.R      # Code for downscaling climate data modified from Rodman et al. 2020 (doi:10.1111/geb.13174)
    |-- WordTemplate                # A word template for generating word files 
|-- Documents                       # Any supporting documents
|-- ecology.csl                     # Citation style language file for generating references following *Ecology* guidlines
|-- .gitignore                      # A gitignore file that specifies intentionally untracked files, which here includes all spatial data
|-- LICENSE                         # License file
|-- README.md                       # README file
|-- Results                         # Results    
```

## License
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)

## Citation
Please cite this work as:
Hart, S.J. and A. Paudel (in prep)
