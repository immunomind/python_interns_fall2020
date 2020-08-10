# Python task


## Background

The gene expression data in single-cell sequencing experiment is summarised as expression matrix where each column corresponds to a gene and each row corresponds to a single cell, the matrix values contain the observed counts. The success of downstream analysis depends on proper quality control that should account for poor quality cells, which should be removed. Failure at this stage may add technical noise which has the potential to obscure the biological signals of interest.

One of the methods for single-cell data preprocessing is filtering the gene expression matrix. Following this approach, researchers filter cells that express too few genes or too many to remove poor quality or doublet cells. Similarly they filter genes that are expressed in too few cells and are thus not informative of the cellular heterogeneity.


## Problems

## Task 1

Write a python function to filter columns of the gene expression matrix, i.e. genes, so that you keep genes that have at least `min_counts` counts or are expressed in at least `min_cells` cells or have at most `max_counts` counts or are expressed in at most `max_cells` cells.

The function receives a matrix (dense or sparse) and returns filtered matrix of the same type. For speedup the function should work on GPU.

## Task 2

A software developer has written a function that filters cells and made a pull request. Review their code and write your comments and suggestions. Please see the file with code attached.


## Data

More than 50,000 cells are typically interrrogated in single-cell experiment and the largest dataset contains gene expression data for more than 1.3 million mouse brain cells. The technology allows researchers to analyse the expression of up to 20,000 genes at the same time. However, due to stochastic nature of gene expression and single-cell analysis the data appears to have many zero elements.

The use of real-world data for the task is not necessary. You can generate a test dataset using the following code:

```
from scipy.sparse import random

S = random(10**5, 2*10**3, density=6e-2)
```


## Directions

1. Develop your code in a separate python file
2. Follow the Google Python Code Style when developing
3. Use type hinting
4. Deposit your code in github repository
5. Set Google Colab notebook for demonstration purposes
6. Create a `code_review.md` file for the second task and put it in the same directory


