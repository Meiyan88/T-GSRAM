# TSGRAM

This is a matlab package for "Imaging Genetics Study Based on a Temporal Group Sparse Regression and Additive Model for Biomarker Detection of Alzheimer’s Disease"


Author: Xiumei Chen;  Meiyan Huang
@Southern Medical University

"TSGRAM" is a package written in Mtalab and the name stands for "Temporal Group Sparse Regression and Additive Model"

This repository contains the following files:

TSGRAM.m is the main function of the TSGRAM package.
generate_options.m is an options struct, with the following fields:  "outputFile":  name of file in outputDir in which to store results;  "lambdaVals.v": vector of lambda values to try;  "numSnps": number of SNPs to select;  "kernelType":type of kernel to use for smoother matrix;  "bandScaler": caler applied to bandwidth used for smoother matrix";  "maxIter": maximum number of backfitting iterations;  "convgTol": tolerance for backfitting stopping criterion;  


Please contact Meiyan Huang (huangmeiyan16@163.com) or Xiumei Chen (chenxiumei97@163.com) for any comments or questions.
