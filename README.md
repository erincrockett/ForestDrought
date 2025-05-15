# ForestDrought
This R code runs the core analyses for the paper: 

Crockett ETH, Qingfeng G, Atkins JE, Sun G, Potter KM, Coztanza J,
Ollinger S, Woodall C, McNulty  S, Trettin C, Holgerson, J, and Xiao J.
Influences of structural and species diversity on forest resistance to drought.

(c) Erin Crockett, 2025
erin.crockett@unbc.ca

Files:
- Analyses_Run.R - provides code to run the core analyses of the paper
- Functions_for_Analyses.R - provides helper functions necessary to run the main code
- ForestDroughtData.csv - provides example data useful to run the analyses. These are not the real data because as noted in the paper, the forest plot coordinates and derivative products are protected for privacy reasons. Folks interested in the exact plot coordinates may contact the US Forest Service to apply for special permissions and go through the data access process.

The code in this repository depends on several R packages that are listed in the code. It also depends on Jarret Byrnes excellent work for spatial corrections: https://github.com/jebyrnes/spatial_correction_lavaan
