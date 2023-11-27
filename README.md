---
title: "hCNVbundles project- ELIXIR UK"
author: "Khaled Jumah, Krzysztof Poterlowicz" 
date: "11/04/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# WP3 - Exploitation of the datasets by the Galaxy Community. 


## 1.  Intergrating a Copy number variant (CNV) detecting tools into Galaxy. 

Although  a number of the CNV detectiong tools have been developed over the recent years (ref Khaled) only a few of them were intergated into the Galaxy
and only couple of them are suported and functional (Khaled put tools table here).
 
Whole exome sequning CNV detection tools according to [reference paper] 

Title|Description|Source|Tool version|CONDA|CONDA source|Wrapper source|Wrapper owner|Wrapper version|Comment|Pull requests (Linked)|Status|Assignees
------|------|------|------|------|------|------|------|------|------|------|------|------
AbsCN-seq|"AbsCN-seq is an R package for estimation of tumor purity, ploidy, and copy number variation (CNV) in next-generation sequencing data."|https://bioinformaticshome.com/tools/cnv/descriptions/AbsCN-seq.html|1|NO|||||||To add|
Accurity|"Accurity is a software tool for inferring tumor cell ploidy, tumor purity, and absolute allelic copy numbers for somatic copy number alterations in whole-genome sequencing (WGS) data."|https://www.yfish.org/display/PUB/Accurity||NO||||| (whole exome may work too)||To add|
ADTEx|"A tool for detection of copy number variation (CNV) in whole-genome exome data from paired or healthy tumor samples. The ADTex algorithm uses Hidden Markov Models (HMM) to predict CNV counts, genotypes, polyploidy, aneuploidy, cell contamination, and baseline shifts. The authors originally named ADTex CoNVEX, but they changed the name because of the conflict with another tool."|https://bioinformaticshome.com/tools/cnv/descriptions/ADTEx.html|2|NO|||||||To add|
aCNViewer|"aCNViewer is a tool for visualization of copy number variation (CNV) and loss of heterozygosity (LOH) in tumor samples. aCNViewer uses the Docker application. Requires: R, Python, Mode, Affymetrix power tools, Ascat, sequenza, and ggplot2."|https://bioinformaticshome.com/tools/cnv/descriptions/aCNViewer.html|2.2|NO|||||||To add|
Bamgineer|A tool to simulate haplotype-phased allele-specific copy number variation (CNV) into a Binary Alignment Mapping (BAM) and cancer CNVs in exome and targeted cell-free DNA sequencing data.|https://bioinformaticshome.com/tools/cnv/descriptions/Bamgineer.html|2|NO|||||||To add|
BubbleTree|"BubbleTree is an R package for the analysis of tumor samples in next-generation sequencing (NGS) data. The BubbleTree algorithm can estimate cancer-causing cell impurity, ploidy, clonality, and allele-specific copy number variation (CNV). It displays the results in a graph format. The Authors claim that BubbleTree outperforms THetA2, ABSOLUTE, AbsCN-seq, and ASCAT tools."|https://bioinformaticshome.com/tools/cnv/descriptions/BubbleTree.html|2.14.0|Yes|https://anaconda.org/search?q=BubbleTree||||||To add|
CANOES|CANOES is a tool for the detection of copy number variation (CNV) in whole-genome exome sequencing data. The CANOES algorithm uses the negative binomial distribution to model read counts and a regression-based method on user-selected reference samples.|https://bioinformaticshome.com/tools/cnv/descriptions/CANOES.html|*** the tool is not currently|NO|||||||To add|

 
 
 
Galaxy training network community [web reference] provide a copmprehensive tutotial [ web link to the tutorial / https://training.galaxyproject.org/training-material/topics/dev/tutorials/tool-from-scratch/tutorial.html]  to instruct the reader in the full process of integrating a tool into Galaxy thorugh the process of  
 - the creation of a bioconda recipe for a new tool
 - writing a Galaxy tool wrapper
 - finally the testing and deployment of this tool into both a local and public Galaxy environment. 

This document present a case study of intergrating  [CNVkit](https://cnvkit.readthedocs.io/en/stable/) into the Galaxy using the above tutorial   

## 2. Benchmarking hCNV detection tools. 
Choosing a specific tool to carry on the analysis require information about   CNVs detection accuracy, execution time and required infrastructure.  
  
Our plan includes :  

1. Do a benchmark test for the CNV detection workflow in Galaxy to measure the run time and detected CNV

2. Compare them with a reference CNV test.    

The reference CNV workflow is available in the reference article below.  
https://www.nist.gov/programs-projects/genome-bottle  
 
The workflow and the tools/code used  
https://github.com/NCBI-Hackathons/TheHumanPangenome/tree/master/MHC/e2e_notebooks 
 

# WP4 - Training and dissemination 
Usually, The preprocessing steps for CNVs data (from the quality control to the mapping step) are the similar for all CNVs detecting workflows. The changes occur in the CNVs detection step. 

Our plan is to create a tutorial that allows the user to understand the CNVs detection process and give guidance on how to choose between the available CNVs tools. 

The process to create this tutorial is by: 

1. View the input requirements for some of the available CNVs tools. 

2. Create a shared prepossessing workflow for different tools data. 

3. Make the tool that requires the longest prepossessing as the main tool in the tutorial and locate the sections where we can introduce the different CBV tools for the data. 

4. add links for additional reading to provide the user with further knowledge on how to use a specific CNV tool 

5. Locate most of the CNV detecting tools available in/outside galaxy 

Our current progress can be found here [tutorial](https://github.com/kpbioteam/training-material/blob/project34/topics/variant-analysis/tutorials/somatic-variant-discovery/tutorial.md) that uses Contol-freec to detect CNVs which also can be the backbone tutorial for this cause.  


