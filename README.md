## CO2-Based Ventilation Rate Estimation :four_leaf_clover: :four_leaf_clover: :four_leaf_clover:
Fall2023  Biostat 615 - Group6 
### Description
This R package provides tools for estimating ventilation rates in indoor environments, crucial for understanding and mitigating airborne disease transmission risks such as COVID-19. It implements the Build-Up Method and the Transient Mass Balance Method, adapted to analyze time-series CO2 data.

### Key Features
- **Build-Up Method**: Utilizes increases in CO2 concentration to estimate air exchange rates without steady state assumptions.
- **Transient Mass Balance Method**: Employs a mass balance model for CO2 in an occupied room to estimate air exchange rates.
- **Robust Data Analysis**: Functions to calculate running averages, variances, and identify periods of CO2 concentration growth or decay.
- **Flexible Input**: Accommodates various inputs including CO2 data, room volume, and number of individuals.
- **Simulated Data Testing**: Capabilities to test the methodology using simulated data under varying conditions.

### Current Progress
- Development and testing of the Build-Up Method algorithm.
- Algorithm for identifying exponential growth and decay in CO2 concentrations.
- Simulated CO2 data generation for algorithm testing.
- Ongoing testing and optimization.

### Key Reference
Batterman S. Review and Extension of CO₂-Based Methods to Determine Ventilation Rates with Application to School Classrooms. Int J Environ Res Public Health. 2017;14(2):145. Published 2017 Feb 4. doi:10.3390/ijerph14020145

### Installation
```R
# Install the package from GitHub
devtools::install_github("SandyZ98/EVCO2")
