# An End-to-End Bioinformatics Pipeline for Single-cell Allele-specific Expression Analysis with Multiple Individuals

Allele-specific expression (ASE) analysis investigates the differential expression of two alleles of the same gene in diploid organisms. ASE provides critical insights into diverse regulatory mechanisms, including cis-regulatory variation, epigenetic modifications, and nonsense-mediated decay. However, existing bioinformatics pipelines for deriving single-cell ASE counts vary significantly across studies. As a result, researchers often expend considerable time and effort on repetitive and technically complex programming tasks to implement these pipelines. To alleviate—and potentially eliminate—this burden, we propose a highly user-friendly bioinformatics pipeline that automatically processes raw sequencing data from multiple individuals and generates single-cell ASE counts for each individual. With minimal user input, the pipeline performs all necessary steps, streamlining the workflow and significantly reducing the manual workload for researchers.

Our work is based on [SALSA pipline](https://github.com/p4rkerw/SALSA). Thanks for their excellent work. 

## Part 1: Overview of the Whole Pipeline
### 1.1 Structure of Single-cell Data
Before discussing the pipeline, we first need to understand the general structure of the data. Due to the extremely large size of single-cell sequencing datasets, it is impractical to store all raw reads (typically saved in .sra files) in a single file. Therefore, researchers usually divide all individuals into smaller groups, each containing multiple individuals—typically 12 to 16. Each such group is referred to as a pool. Within each pool, the same set of individuals may be sampled multiple times. For instance, it is possible that all individuals in Pool 1 are measured across 20 separate experiments, with each experiment referred to as a run.

We use the [OneK1K dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196830) (GSE196830) as an example. (In this project, we focus exclusively on PBMCs and do not analyze whole blood samples.) This dataset contains a total of 77 pools, each consisting of multiple individuals (12, 14, or 16 per pool). For each pool, data from 20 runs are available.

This is the general structure of single-cell sequencing data collected from multiple individuals. Our pipeline is specifically designed to process such data and produce allele-specific read counts for each individual.

<span style="color:red">(figrues for the structure)</span>

### 1.2 Bioinformatics Pipeline
Our pipeline is built upon the SALSA pipeline, an excellent bioinformatics workflow for scASE analysis. However, our pipeline introduces two key improvements:

1. The original SALSA pipeline is designed for allele-specific expression analysis of a single individual and cannot be directly applied to datasets containing multiple individuals. In contrast, our pipeline is specifically designed for large-scale single-cell datasets with multiple individuals, following the structure described in Section 1.1. The SALSA pipeline is integrated as a core component within our broader framework.

2. The SALSA pipeline is not user-friendly; users often need to spend significant time and effort to run it successfully. Although our pipeline is more complex in functionality, it is significantly more user-friendly. It offers an accessible, streamlined solution for deriving ASE counts directly from raw RNA sequencing data.

Our pipeline consists of seven main steps. For each step, users only need to customize file paths according to their data and then wait for the process to complete.

In Section 2, we describe the required preparation:

- In Section 2.1, we prepare cell barcodes for each individual.

- In Section 2.2, we set up the required programming environment.

Once the preparation is complete, we proceed to the main pipeline in Section 3:

- In Section 3.1, we download the data from GEO. After this step, we obtain three FASTQ files (read1, read2, and index) for each run.

- In Section 3.2, we split the FASTQ files by individual.

- In Section 3.3, we recombine reads belonging to the same individual across different runs.

- In Section 3.4, we use Cell Ranger to align each individual’s data to the reference genome.

- Finally, in Section 3.5, we obtain single-cell ASE counts. This step largely follows the SALSA pipeline with some modifications to accommodate multi-individual data.

<span style="color:red">(A schematics figrues)</span>

## Part 2: Environment and Data Preparation 
In this document, we will use the [OneK1K dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196830) (GSE196830) as an example. We will show how to use our pipeline to process the **Pool11** and get the single-cell ASE for three individuals (847_848, 985_986, 997_998). 

<span style="color:red">(Talk about requirement)</span>

On your server, you fist need to create a folder only for our pipeline. And you need to move to this dictionary. For example: 

```
# we assume you will work on this folder. You can change to your reference. 
mkdir /home/user/mypipeline  
cd /home/user/mypipeline

# create two folder
mkdir cellbarcodes 
mkdir reference
```

In the folder `/home/user/mypipeline`, we need to folder: `cellbarcodes` and `reference`. In `cellbarcodes`, we will save cell barcodes for each individuals. We will show some example below. In `reference`, we will save all packages for this pipeline. Remember, the path is very important in our pipeline. Make sure your current path is under `/home/user/mypipeline` and make sure you can correctly change some paths below. 



### 2.1 Cell Barcodes 
Our pipeline is very user-friendly, most of work will be done automatically. People only need to do some work manully. One is work is to downlaod cell barcodes. So what is cell barcodes ? We can first take a look at one example: 

<div align="center">
  <img src="figures/cell_barcodes.png" alt="barcodes" width="200"/>
</div>

Here is an example of cell barcodes. The first column represents the **individual ID**, and the second column contains a short RNA sequence, referred to as the **cell barcode**. Each individual is associated with multiple barcodes, as a single individual typically contributes thousands of cells in high-throughput RNA sequencing platforms (such as the Illumina NovaSeq 6000).

In other words, while a single experiment can capture expression profiles from millions of cells, cell barcodes are essential for determining the individual origin of each cell. This mapping between barcodes and individuals enables downstream allele-specific expression analysis at the single-cell level.

Next question, how can we get cell barcodes ? 
In this document, you do not need to download it. Because we have prepare it for you under the folder **cellbarcodes** . But if you want to study other pools of OneK1K data set, you can directly download them online. For example, if you want to study Pool12, go to [the GEO websites for OneK1K Pool12](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5899884). Go to the **Supplementary file**, and find the file named `GSM5899884_OneK1K_scRNA_Sample12_Individual_Barcodes.csv.gz`, which is what we want. You need to decompress it by the command: 

```
gunzip GSM5899884_OneK1K_scRNA_Sample12_Individual_Barcodes.csv.gz
```

You can do it for other pools in a similar way. You can also use the same way to download cell barcodes of other dataset. But **remember**: 1. you only need to download cell barcodes for each pool. Because the set of individuals is the same in the same pool across different runs. 2. Upload the file to the folder: `/home/user/mypipeline/cellbarcodes`

If you have unzip the files and upload it to your server, all have been done. Very easy, right ? 

### 2.2 Create Environment


## Part 3: Get Multiple Individual ASE Count

### 3.1 Download .SRA Files

### 3.2 Split FASTQ Files 


### 3.3 Recombine FASTQ Files 


### 3.4 Cellranger 


### 3.5 scASE counts









































