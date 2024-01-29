# ABDS
ABDS: Tool Suite for Analyzing Biologically Diverse Samples

## CONTENTS
### Mechanism-integrated Group-wise Pre-Imputation (MGpI):
- 'MGpI.R' houses all the MGpI functions, with a quick intro to the input and parameters.
- 'MGpI_110554.R' is used for conducting the comparison experiment on the real-world GSE110554 data.
- 'MGpI_19380.R' is used for conducting the comparison experiment on the real-world GSE19380 data.
- 'MGpI_LAD45.R' is used for conducting the comparison experiment on the real-world LAD45 proteomic data.
- 'helper functions.R' is sourced in the experiments mentioned above.
- 'imputation' folder houses all the peer methods mentioned in our paper for conducting the comparison experiment.
- 'LAD45.csv' is used in 'MGpI_LAD45.R' as the input data.

### Enumerated Cosine-based One-sample Test (eCOT):
- 'eCOT.R' houses all the eCOT functions.
- 'eCOT_Experiment.R' is used to conduct experiments mentioned in our paper related to eCOT.
- 'scatter_plot.R' is a helper function for illustrating the results in 'eCOT_Experiment.R'.

### Unified HeatMap (uniHM):
- 'uniHM.R' houses all the unified heatmap functions.
- 'example_uniHM.R' is used to conduct experiments mentioned in our paper related to uniHM.
- 'sample_data_uniHM.xlsx' is the toy data to demonstrate uniHM.


