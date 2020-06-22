# SCLC-NK-scoring
Computational scripts to accompany the manuscript **Harnessing natural killer immunity in metastatic small cell lung cancer** by SA Best *et al.* (2020). *Journal of Thoracic Oncology*.
- DOI: [10.1016/j.jtho.2020.05.008](https://www.sciencedirect.com/science/article/pii/S1556086420303944)

For further information about the manuscript, please contact:
  - Corresponding author, Dr Kate Sutherland :: sutherland.k (at) wehi.edu.au
  - Lead author, Dr Sarah Best :: best (at) wehi.edu.au

For further information on the computational analyses, please contact:
  - Bioinformatics Lab Head, A/Prof Melissa Davis :: davis.m (at) wehi.edu.au
  - Bioinformatics author, Dr Joe Cursons :: joseph.cursons (at) monash.edu
  
  
## Overview
Please refer to the manuscript listed above for details. Briefly, small cell lung cancer (SCLC) is an aggressive subtype of lung cancer with particularly poor patient outcomes. Immunotherapies which augment the activity of a patient's immune system have revolutionised cancer treatment in recent years, and in this study we explore the ability of natural killer (NK) cells to target small cell lung cancer. The script hosted on this repository was used for the analysis of related public data sets and generation of several figure panels included within the manuscript.


## Scripts
Computational analyses generating Figures 1C-F, Figures 2A-C, and Supplementary Figure 1C are all contained within the script *Best_2020_SCLC.py*. This script has been split into the *PathDir*, *PreProcess* and *Plot* classes as detailed below.

### PathDir
This class can be used to specify data and/or output locations, or default to local locations with some user prompts.


### PreProcess

### Plot

* figure_one()

* figure_two()

* supp_fig_one_scatter()

## Data
This project analyses several public data sources on SCLC, including a patient cohort study by George *et al* (2015), and small cell lung cancer data from the Cancer Cell Line Encyclopedia (CCLE).

### Data Sources
George, J *et al.* (2015). *Nature*. **524**, pp. 47–53.
- DOI: [10.1038/nature14664](https://www.nature.com/articles/nature14664)
- In particular, the [supplementary data tables](https://static-content.springer.com/esm/art%3A10.1038%2Fnature14664/MediaObjects/41586_2015_BFnature14664_MOESM72_ESM.xlsx).
  - Note that the transcript abundance data have been embedded within an Excel spreadsheet and accordingly a subset of gene names have been corrupted due to date autoconversion.
  
Mollaoglu, G *et al.* (2017). MYC Drives Progression of Small Cell Lung Cancer to a Variant Neuroendocrine Subtype with Vulnerability to Aurora Kinase Inhibition. *Cancer Cell*.
  - DOI: [10.1016/j.ccell.2016.12.005](https://www.cell.com/cancer-cell/fulltext/S1535-6108(16)30600-6)

Rudin, CM *et al* (2019). Molecular subtypes of small cell lung cancer: a synthesis of human and mouse model data. *Nature Reviews: Cancer*. 
  - DOI: [10.1038/s41568-019-0133-9](https://www.nature.com/articles/s41568-019-0133-9)
  - **NB**: This analysis uses sample classifications from Rudin et al where available and an author correction was published for this material (https://www.nature.com/articles/s41568-019-0164-2), which has an updated [Supplementary Table](https://static-content.springer.com/esm/art%3A10.1038%2Fs41568-019-0164-2/MediaObjects/41568_2019_164_MOESM1_ESM.xlsx)
  
Cursons, J. *et al*. (2019). A Gene Signature Predicting Natural Killer Cell Infiltration and Improved Survival in Melanoma Patients. *Cancer Immunology Research*. 
  - DOI: [10.1158/2326-6066.CIR-18-0500](https://cancerimmunolres.aacrjournals.org/content/early/2019/05/31/2326-6066.CIR-18-0500)
    - The natural killer cell signature used for gene set scoring is given within [Supplementary Table 1](https://cancerimmunolres.aacrjournals.org/highwire/filestream/37099/field_highwire_adjunct_files/0/205802_2_supp_5410854_prn3pl.xlsx) of this manuscript.


### Dependencies
**singscore**
- [Python library (pysingscore)](https://github.com/DavisLaboratory/PySingscore)
- [R/Bioconductor package](https://www.bioconductor.org/packages/release/bioc/html/singscore.html)
- [Corresponding Manuscript](https://dx.doi.org/10.1186/s12859-018-2435-4)
	
## Acknowledgements
The authors would like to acknowledge the following groups; without their contributions this work would not have been possible.

### Patient donors
All medical research requires contributions from patients who have donated tissue samples and consented to sequencing of tumor samples for scientific research.

### Open Source software developers
As demonstrated by the large number of libraries imported in the header of the associated python script, this work is heavily dependent upon the work of Open Source developers. We apologise if we have missed citations for any of the associated packages, but in particular we would like to thank contributers for the following projects:
- matplotlib
- numpy
- pandas
- scikit-learn
- scipy
- seaborn