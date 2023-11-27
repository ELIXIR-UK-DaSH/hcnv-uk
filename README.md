---
title: "hCNVbundles project- ELIXIR UK"
author: "Khaled Jumah, Krzysztof Poterlowicz" 
date: "11/04/2022"
output: html_document
---


# WP3 - Exploitation of the datasets by the Galaxy Community. 


## 1.  Intergrating a Copy number variant (CNV) detecting tools into Galaxy. 

Although  a number of the CNV detectiong tools have been developed over the recent years (ref Khaled) only a few of them were intergated into the Galaxy
and only couple of them are suported and functional (Khaled put tools table here).
 
Whole exome sequning CNV detection tools according to [reference paper] 

Title|Description|Source|Tool version|CONDA|CONDA source|Wrapper source|Wrapper owner|Wrapper version|Comment|Pull requests (Linked)|Status|Assignees|
|------|------|------|------|------|------|------|------|------|------|------|------|------|
|AbsCN-seq|"AbsCN-seq is an R package for estimation of tumor purity, ploidy, and copy number variation (CNV) in next-generation sequencing |data."|https://bioinformaticshome.com/tools/cnv/descriptions/AbsCN-seq.html|1|NO|||||||To add|
|Accurity|"Accurity is a software tool for inferring tumor cell ploidy, tumor purity, and absolute allelic copy numbers for somatic copy number alterations in whole-genome sequencing |(WGS) data."|https://www.yfish.org/display/PUB/Accurity||NO||||| (whole exome may work too)||To add|
|ADTEx|"A tool for detection of copy number variation (CNV) in whole-genome exome data from paired or healthy tumor samples. The ADTex algorithm uses Hidden Markov Models (HMM) to predict CNV counts, genotypes, polyploidy, aneuploidy, cell contamination, and baseline shifts. The authors originally named ADTex CoNVEX, but they changed the name because of the conflict with another tool."|https://bioinformaticshome.com/tools/cnv/descriptions/ADTEx.html|2|NO|||||||To add|
|aCNViewer|"aCNViewer is a tool for visualization of copy number variation (CNV) and loss of heterozygosity (LOH) in tumor samples. aCNViewer uses the Docker application. Requires: R, Python, Mode, Affymetrix power tools, Ascat, sequenza, and ggplot2."|https://bioinformaticshome.com/tools/cnv/descriptions/aCNViewer.html|2.2|NO|||||||To add|
|Bamgineer|A tool to simulate haplotype-phased allele-specific copy number variation (CNV) into a Binary Alignment Mapping (BAM) and cancer CNVs in exome and targeted cell-free DNA sequencing data.|https://bioinformaticshome.com/tools/cnv/descriptions/Bamgineer.html|2|NO|||||||To add|
|BubbleTree|"BubbleTree is an R package for the analysis of tumor samples in next-generation sequencing (NGS) data. The BubbleTree algorithm can estimate cancer-causing cell impurity, ploidy, clonality, and allele-specific copy number variation (CNV). It displays the results in a graph format. The Authors claim that BubbleTree outperforms THetA2, ABSOLUTE, AbsCN-seq, and ASCAT tools."|https://bioinformaticshome.com/tools/cnv/descriptions/BubbleTree.html|2.14.0|Yes|https://anaconda.org/search?q=BubbleTree||||||To add|
|CANOES|CANOES is a tool for the detection of copy number variation (CNV) in whole-genome exome sequencing data. The CANOES algorithm uses the negative binomial distribution to model read counts and a regression-based method on user-selected reference samples.|https://bioinformaticshome.com/tools/cnv/descriptions/CANOES.html|*** the tool is not currently|NO|||||||To add|
| CODEX2           | The tool uses depth information from the sequence to determine the presence or absence of a CNV                                      | [GitHub](https://github.com/yuchaojiang/CODEX2)                                               |              | NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| Control-FREEC    | A tool for assessing copy number and allelic content using next-generation sequencing data.                                        | [GitHub](https://github.com/BoevaLab/FREEC)                                                  | v11.6        | Yes   | [CONDA](https://anaconda.org/search?q=Control-FREEC%09) | [GitHub](https://github.com/galaxyproject/tools-iuc/tree/master/tools/freec)  | iuc           | 11.6             | Up-to-date                               |                        |             |             |
| CoNIFER          | Uses exome sequencing data to find copy number variants (CNVs) and genotype the copy-number of duplicated genes                      | [SourceForge](https://sourceforge.net/projects/conifer/)                                      | 0.2.1        | Yes   | [CONDA](https://anaconda.org/search?q=CoNIFER)      |                                                      |               |                  | To add                                   |                        |             |             |
| CONTRA           | A tool for copy number variation (CNV) detection for targeted resequencing data such as those from whole-exome capture data.         | [CONTRA](https://contra-cnv.sourceforge.net/)                                                | v2.0.6       | NO    |                                                   | [Toolshed](https://toolshed.g2.bx.psu.edu/repository) | fcaramia      | 2.4              | To update                               |                        |             |             |
| CoNVaDING        | A tool for identification of copy number variation (CNV) specifically in next-generation sequencing (NGS) data.                      | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/CoNVaDING.html) | 1.2.0        | NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| CloneCNA         | CloneCNA is a tool to detect copy number variation (CNV) in heterogeneous tumor samples in whole-genome exome sequencing data.        | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/CloneCNA.html)   | 2            | NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| CNVkit           | is a Python library and command-line software toolkit to infer and visualize copy number from high-throughput DNA sequencing data. | [GitHub](https://github.com/etal/cnvkit)                                                      | 0.9.9        | Yes   | [CONDA](https://anaconda.org/search?q=CNVkit)    |                                                      |               |                  | Up-to-date                               | khaled196               |             |             |
| cn.mops          | cn.mops (Copy Number estimation by a Mixture Of PoissonS) is an R package pipeline for analysis of copy number variation (CNV)        | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/cn-mops.html)     | 1.30.0       | Yes   | [CONDA](https://anaconda.org/search?q=cn.mops)  |                                                      |               |                  | To add                                   |                        |             |             |
| CNspector        | CNspector is a web-based tool for clinical diagnosis and visualization of copy number variation (CNV) using next-generation sequencing (NGS) data. | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/CNspector.html) || NO    |                                                   |                                                      |               |                  | Special methods for WES                  |                        | To add      |             |
| CNVfinder        | CNVfinder is a tool to detect copy number variation in whole-genome exome sequencing data generated using amplicon-based enrichment technologies. | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/CNVfinder.html)  | 1.0.0b5      | NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| CNVnator         | A tool for discovery and characterization of copy number variation (CNV) in population genome sequencing data.                      | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/CNVnator.html)  || Yes   | [CONDA](https://anaconda.org/search?q=CNVnator) | [Toolshed](https://toolshed.g2.bx.psu.edu/repository) | fanruimeng    |                  | It was used in a paper to detect CNVs in WES data | Neglected               |             |             |
| cnvOffSeq        | cnvOffSeq is specifically designed for the detection of intergenic copy number variation (CNV) using off-target exome sequencing data. | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/cnvOffSeq.html)   | 0.1.2        | NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| CN_Learn         | CN_Learn is a tool for copy number variation (CNV) detection. The algorithm integrates multiple CNV detection algorithms and learns to identify CNVs based on validated CNVs. | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/CN_Learn.html)  || NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| DeAnnCNV         | DeAnnCNV (Detection and Annotation of Copy Number Variations) is a web-based tool to detect and annotate copy number variation (CNV) in whole genome exome sequencing. The algorithm is based on GPHMM tool (see links). | [BioinformaticsHome](https://bioinformaticshome.com/tools/cnv/descriptions/DeAnnCNV.html) || NO    |                                                   |                                                      |               |                  | To add                                   |                        |             |             |
| EXCAVATOR2       | EXCAVATOR2 is a tool for detecting copy number variation (CNV) in whole-genome exome sequencing data. The EXCAVATOR2 algorithm uses a hidden Markov model (HMM) and a method that classifies genomic regions into five CNV states. EXCAVATOR2 is an enhanced version of EXCAVATOR (see links). |

 
 
 
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


