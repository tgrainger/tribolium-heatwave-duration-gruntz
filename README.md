This repository contains the data and R code used to reproduce the 
analyses and figures presented in:

Gruntz and Grainger (submitted) "Persistent effects of experimental heatwave duration on population dynamics"

## Repository Structure

Data:

population counts  
body size 
sex ratios     
fecundity   

Scripts:

analysis 

Output:

Figures


## Software Requirements

Analyses were conducted in:

- R 4.5.0
- RStudio 2025.05

Required packages:

- tidyverse
- lme4
- lmerTest
- glmmTMB
- DHARMa
- car
- emmeans
- cowplot
- ggplot2
- broom
- purrr
- ggpubr
- patchwork


## Data Descriptions

#population counts.csv

-count: Census number (1, 2, 3, 4, 5)
-monthssinceheatwave: Months passed between end of heatwave and measurement (0, 1, 3, 5, 7)
-treatment: Heatwave duration treatment (control, 5 day, 10 day, 15 day)
-replicate: Population replicate (p1 to p10)
-live: Number of live adult beetles
-dead: Number of dead adult beetles
-birthswithnegs: births calculated according to Formula 1 in manuscript
-births: biths with negative values removed (counting error)

# fecundity.csv

-count: Census nunmber (1, 2, 3, 4, 5)
-monthssinceheatwave: Months passed between end of heatwave and measurement (0, 1, 3, 5, 7)
-treatment: Heatwave duration treatment (control, 5 day, 10 day, 15 day)
-replicate: Population replicate (p1 to p10)
-femalerep: Replicate female (1, 2, 3)
-eggs: Number of eggs laid in 48 hours

# body size.csv

-monthssinceheatwave: Months passed between end of heatwave and measurement (0, 1, 3, 5, 7)
-treatment: Heatwave duration treatment (control, 5 day, 10 day, 15 day)
-replicate: Population replicate (p1 to p10)
-sex: Sex of beetles (male or female)
-individualreplicate: Replicate individual (1-5)
-weighting: Weight in grams

# sex ratios.csv

-monthssinceheatwave: Months passed between end of heatwave and measurement (0, 1, 3, 5, 7)
-treatment: Heatwave duration treatment (control, 5 day, 10 day, 15 day)
-replicate: Population replicate (p1 to p10)
-males: Number of male beetles
-females: Number of female beetles

## Contact

Tess Grainger

Department of Integrative Biology  
University of Guelph

Email: tess.grainger@uoguelph.ca