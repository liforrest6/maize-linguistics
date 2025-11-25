## Human-Maize co-evolution project - linguistics analyses

Data and scripts for the linguistics analyses of the human-maize co-evolution project.
Scripts are generally sorted into four stages:

1. Generate data
* createDistanceMatrices.R
* createDistanceMatrices_Haynie.R
* generateGlobalGrid.R

2. HPC Scripts

### For GWAS
* prepare_GWAS_DATA.SH
* makeGWASbinaryfiles.R
* transpose.awk
* One_hot_encoding.R
* language-gwas.sh

### For plotting
* getGenotype.sh
* makeManhattanPlot.R
* submit_makeManhattanPlot.sh

## For FEEMS
* prepare_FEEMS_data.sh
* run_FEEMS.py
* submit_run_FEEMS.sh

## For regression
* prepare_regression_data.sh
* regressionSNPs.sh
* submit_regressionSNPs.sh

### Other analysis
* popGenScripts.sh
* performPCA.R

3. Load data
* languageFunction.R
* loadData.R

4. Analyses
* gwasAnalysis.R
* glottologCentroidValidaation.R
* languageRegression.R

5. Plots
* plotMaps.R
* plotManhattans.R
* makeFigurePlots.R

