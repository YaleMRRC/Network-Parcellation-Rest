# NetworkParcellation

<p align="center">
	<img src ="images/rest_individualzied_networks.jpg" width="480" height="347" />
</p>

These set of functions are used to generate the results in paper "An exemplar-based approach to individualized parcellation reveals the need for sex specific functional networks," Salehi et al., NeuroImage 2018. This paper uses data from two resting-state sessions (REST1 and REST2) of the Human Connectome Project (HCP) S900 release. This repository includes MATLAB functions to:
1) Import HCP data into MATLAB's workspace: **import_HCP_data.m**
2) Generate individualized resting-state parcellations: **Parcellation_rest.m**
3) Computes and compares internal clustering quality measures such as Dunn index and Davies-Bouldin index (DB): **Comparison_DB_Dunn.m**
4) Compute the reproducibility of tha parcellations between different sessions of resting-state: **REST1_REST2_reproducibility.m**
5) Develop predictive models based on node-to-network assignment (NNA) vectors using gradient boosting machines (GBM): **Exemplar_based_predictive_analysis.m**
6) Compute the importance of features (here nodes) in the predictive model's performance: **Feature_importance_analysis.m**


The main function in this set is **"Parcellation_rest.m"**, which generates the individualzied functional networks from the node-level (i.e., parcellated by grouping voxels into nodes) fMRI data. You can execute this function via the following:

1. Using the link on the top right corner of this page, clone this repository to your local computer using the following command: 
```bash
git clone https://github.com/YaleMRRC/Network-Parcellation.git
``` 
2. In MATLAB add the path to your script or in your workspace: 
```matlab
addpath('localpath/Network-Parcellation/')
```
3. Load in HCP resting-state data from your local directory (datadir) by running the following function:
```matlab
[M, HCP_subj] = import_HCP_data(datadir)
```
4. Run the following: 
```matlab
Parcellation_rest(M,HCP_subj,label_134_cort,K_max)
```
   where **label_134_cort** is a vector indicating whether each node belongs to cortex, subcortex, or cerebellum. This is used to automatically include/exclude specific parts of the brain from the analysis. **K_max** indicates the number of networks.

The rest of the functions are either auxiliary functions that are called by the meain function, or are used for other analyses in the paper. Below you can find the details about the maind functions in this set:


For further questions please raise an issue [here](https://github.com/YaleMRRC/Network-Parcellation/issues).
