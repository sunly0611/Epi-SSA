# Epi-SSA Algorithm

## Introduction

The Epi-SSA algorithm is an algorithm based on the Sparrow Search Algorithm for detecting epistasis interactions in GWAS data. The article corresponding to Epi-SSA is currently being submitted for publication.

## Usage

### On Simulated Data

`java -jar epi-ssa.jar -maxG 2000 -n 60`

### On Real Data

`java -jar epi-ssa.jar -maxG 2000 -type 1 -pathIn example`

### Parameters Of DE-GWO

|parameter name|description|
|:------:|:------|
|help|display help description|
|pathIn|path of the GWAS data, default 'data.txt'.|
|pathO|path of the file recording the results, default 'result.txt'.|
|type|type of the GWAS data, 0 for simulated, 1 for PLINK tped format, default 0.|
|seed|seed of random, default 0.|
|d|the maximum order of SNP combinations, default -1.|
|maxL|the maximum length of contingency table, default -1.|
|cG|threshold, default cG=0.05.|
|maxG|the max iterations of Sparrow Search Algorithm.|
|n|number of the sparrows in the algorithm.|
|pd|producer ratio, default 0.4.|
|sd|the ratio of sparrows who perceive the danger, default 0.2.|
|st|safe threshold, default 0.8.|



