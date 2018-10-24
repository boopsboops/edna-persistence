# edna-persistence

Final data and code for the marine eDNA persistence paper published in _Communications Biology_:

Collins, R. A., Wangensteen, O. S., O'Gorman, E. J., Mariani, S., Sims, D. W., Genner, M. J. (2018) Persistence of environmental DNA in marine systems. _Communications Biology_. [https://doi.org/10.1038/s42003-018-0192-6](https://doi.org/10.1038/s42003-018-0192-6).

### Contents

* **`code/`** - R and shell scripts to run analyses.
    - `libs.R` - load required R packages
    - `map.R` - create map of study area 
    - `MFEprimer.sh` - run _in silico_ PCR
    - `model_selection.R` - fit eDNA decay model and create plots
    - `primer3.sh` - generate candidate primer pairs
    - `processing.R` - clean and format the raw data
    - `rates.R` - literature review of eDNA decay rates 
    - `species_list.R` - generate a list of UK fish/crustacean species and download COI data from GenBank  
    - `summary_stats.R` - generate descriptive statistics

* **`data/`** - Data and settings used in the analyses.
    - `crab.p3` - Primer3 settings for crabs
    - `efficiency.csv` - qPCR efficiency data 
    - `literature_rates.csv` - eDNA decay rates from the literature
    - `qPCR_data.csv` - qPCR data for main experiment
    - `shanny.p3` - Primer3 settings for shanny
    - `treatments_tanks.csv` - treatments and parameters for each experimental tank

* **`temp/`** - Temporary file directory that is not committed to the repository, but needs to be created locally first, to run the scripts and see the output. Ignored by git.