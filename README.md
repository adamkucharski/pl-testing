# pl-testing
Analysis of Premier testing data

### Quick start guide

First, set local path in R to GitHub directory, e.g.:
`
setwd("~/Documents/GitHub/pl-testing/")
`
Code to load data and generate the plots is in `scripts/run_analysis_script.R`. This loads the following datasets:

* `incidence_data.csv` - Weekly PCR tests and positives in the testing dataset
* `incidence_data_lft.csv` - Weekly LFT tests and positives in the testing dataset
* `ct_values.csv` - Individual participant Ct values (arranged by unique `User_ID`), along with time since first positive (`diff_date`)
* `first_test_wt.csv` - Individual participant Ct values from positives during the clear pre-Alpha period (i.e. before 1st Nov 2020). Note `User_ID` is unique to this dataset and not matched with above datasets.
* `first_test_alpha.csv` - Individual participant Ct values from positives during the clear Alpha period (i.e. 1st Jan 2021 to 1st Apr 2021). Note `User_ID` is unique to this dataset and not matched with above datasets.
* `proportion_data_pcr.csv` - Number tested and number positive by PCR, by day since first positive PCR test (day = 0), for participants who tested positive where repeat tests are available. Results outside 28 days before or 25 days after not shown (i.e. NA).
* `proportion_data_lft.csv` - Number tested and number positive by LFT, by day since first positive PCR test (day = 0), for participants who tested positive where repeat tests are available.

### Reference
[Kucharski AJ, Russell TW, Hellewell J et al. SARS-CoV-2 Dynamics in the Premier League Testing Program, United Kingdom. Emerg Infect Dis. 2024](https://wwwnc.cdc.gov/eid/article/30/9/24-0853_article)
