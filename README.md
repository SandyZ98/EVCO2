## CO2-Based Ventilation Rate Estimation :four_leaf_clover: :four_leaf_clover: :four_leaf_clover:
Fall2023  Biostat 615 - Group6 
### Description
The primary problem our project addresses is the absence of accessible, transparent tools for calculating ventilation rates in indoor environments. Accurate estimations are critical for evaluating disease transmission risks and devising effective mitigation strategies in these spaces. Existing methods, which often involve tracer gasses or sophisticated airflow measurement tools, are impractical for widespread use due to their complexity, expense, and need for specialized equipment. Further, There is a noticeable lack of R packages that can calculate ventilation rates from gaseous concentration data within a given room.

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
```
### How to use
Please check the final report(https://docs.google.com/document/d/1MH8s1ZH_H15_jIEA2GDgmWwwx_fedcp-hY88ixm1nRc/edit) and the documentation and examples for each function in the package.
