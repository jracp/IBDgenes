---
title: "IBDgenes"
authors: "Javad Rahimipour Anaraki and Hamid Usefi"
date: '21/01/19'
---

## Notice
The code `PFS.m` selects the most relevant and independant genes associated to Inflammatory Bowel Disease (IBD). If you need more details and explanation about the algorithm, please contact [Javad Rahimipour Anaraki](https://jranaraki.github.io/) or [Hamid Usefi](http://www.math.mun.ca/~usefi/). 

## Use case
To determine the most important features using the algorithm described in "A Comparative Study of Feature Selection Methods on Genomic Datasets" by Javad Rahimipour Anaraki and Hamid Usefi 

Here is a link to the paper: https://ieeexplore.ieee.org/document/8787392

## Compile
This code can be run using MATLAB R2006a and above

## Run
To run the code, open `PFS.m` and choose a dataset to apply the method to. The code strats reading the selected dataset using `readLargeCSV.m` written by [Cedric Wannaz](https://www.mathworks.com/matlabcentral/profile/authors/1078046-cedric-wannaz). Then it rank features and returns resulting classification accuracy by `cAcc.m` using SVM classifier. The IBD dataset is downloaded from GEO under accession number [GSE3365](https://www.ncbi.nlm.nih.gov/pubmed/16436634).

## Note
 - Dataset should have no column and/or row names, and the class values should be all numeric
