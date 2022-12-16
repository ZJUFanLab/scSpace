# scSpace reproducibility
Code to reproduce the analysis and figures in the paper [Reconstruction of the cell pseudo-space from 
single-cell RNA sequencing data with scSpace](https://www.biorxiv.org/content/10.1101/2022.05.07.491043v1). The scSpace
python package is located [here](https://github.com/ZJUFanLab/scSpace).

## Obtaining Datasets
The data used and/or generated in this study can be found on [Google Drive](https://drive.google.com/drive/folders/1_9qS78yDi0agT4x5HQuldg3I1a2WWdYc?usp=sharing).

## Running scSpace
The analysis scripts to produced scSpace results can be accessed [here](scripts). For example, the script 
[constructSimulations.R](scripts/constructSimulations.R) is used to generate simulated scRNA-seq and spatial transcriptomics
data in benchmarking step.

## Generating Main Figures
We provide R Markdown files that were used to create the main figures, 
the results consistent with the figures can be accessed [here](output):
* [Performance of scSpace compared with other methods on simulated datasets](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure1.html) (Figure 1)
* [Validation of scSpace on human DLPFC 10x Visium data](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure2.html) (Figure 2)
* [Validation of scSpace on mouse primary visual cortex STARmap data](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure2.html) (Figure 2)
* [Validation of scSpace on human embryonic heart scRNA-seq data](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure3.html) (Figure 3)
* [Validation of scSpace on human MTG scRNA-seq data](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure4.html) (Figure 4)
* [Application of scSpace on human melanoma scRNA-seq data](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure5.html) (Figure 5)
* [Application of scSpace on human COVID-19 scRNA-seq data](https://raw.githack.com/ZJUFanLab/scSpace/master/AnalysisPaper/figures/figure6.html) (Figure 6)

## Supplementary Figures
TODO:
* [Figure S1]()
* [Figure S2]()
* [Figure S3]()
* [Figure S4]()
* [Figure S5]() 
