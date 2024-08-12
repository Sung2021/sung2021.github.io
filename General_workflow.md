
This report outlines the general workflows I utilize for the analysis of sequencing data across various platforms. I generally follow these guidelines for most of my work and incorporate additional analyses based on the specific goals and discussions related to each project.          

<br>

**scRNA-seq pipeline**   

**RNA-seq pipeline**  

**ChIP-seq/Cut&Run/ATAC-seq pipeline**   

**Examples of Data Pre-processing**   

* Fastq files QC: Fastqc reports
* 10X Genomics datasets (CellRanger)   
* Cellranger vdj (TCR-seq)   
* STAR alignment for RNA-seq
   

<br>
<br> 


# scRNA-seq pipeline  
<br>

* Fastq files QC: Fastqc reports   
* 10X Genomics datasets (CellRanger) : countmatrix   
* Other formats: Drop-seq Tools, Salmon, STARsolo   
* Demultiplexing (UMI-tools)   
* Doublets-removal: scrublet, DoubletFinder   
* Filtering by custom parameter    
* Normalization and scaling : Seurat(SCTransform)   
* Cell cycle phase : Seurat   
* Integration and clustering : Seurat(CCA, Harmony)   
* Finding markers: Seurat, custom ML algorithm (XGBoost, Random Forest), DESeq2 and edgeR   
* Visualization : Seurat functions, ggplot, plotly, shiny, heatmap, custom scripts   
* Inter-cluster differential expression : MAST, DESeq2   
* Imputation of gene expressions : MAGIC, scImpute   
* Trajectories : Monocle3, destiny   
* pathway analysis: GSEA, scGSEA, clusterprofiler    
* Network analysis: WGCNA, SCENIC   
* Cell-Cell interaction : Cellchat    
* RNA velocity : Velocyto    
* scTCR/BCR-seq : VDJtools, immunarch   

<br>


# RNA-seq pipeline 
<br>

* Fastq files QC: Fastqc reports   
* STAR,Salmon: countmatrix   
* Normalization and scaling (if necessary)   
* TPM, FPKM calculation    
* PCA plots   
* Correlation analysis of selected geneset  
* DEG analysis   
* GSEA   
* Pathway analysis: GSEA, clusterprofiler    
* K-means clustering of DEGs   
* Visualization    



# ChIP-seq/Cut&Run/ATAC-seq pipeline 
<br>

* Fastq files QC: Fastqc reports    
* Mapping: bowtie2    
* Peak calling : macs2, SICER2    
* bigwig file visualization :  Deeptools, IGV, UCSC genome browser       
* Custom analysis with bed files : BEDtools, GenomicRanges, rtracklayer, ChIPseeker      
* Differentially Expressed Peaks : DESeq, edgeR       
* TSS enrichment : Deeptools (computeMatrix)    
* Motif discovery and analysis : HOMER, GREAT    


<br>
Tools for evaluating differential enrichment   
<img src="https://hbctraining.github.io/Intro-to-ChIPseq/img/diffpeaks-software.png">  
<span style="font-size:70%">https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html</span>
<br>

# Examples of Data Pre-processing   
<br>

* Fastq files QC: Fastqc reports   

```{r, eval=F, echo=TRUE}
mkdir fastqc_result
cd $path/

for input in $(ls *.fastq.gz) ;do fastqc $path/${input} ; done
```


* Trimming  

```{r, eval=F, echo=TRUE}
# gzip files only
cd $path/
mkdir trimmed
for input in files; do java -jar path/to/trimmomatic/trimmomatic-0.39.jar PE -threads 10 $path/${input}.R1.fastq.gz $path/${input}.R2.fastq.gz $path/trimmed/${input}.R1.paired.fastq.gz $path/trimmed/${input}.R1.unpaired.fastq.gz $path/trimmed/${input}.R2.paired.fastq.gz $path/trimmed/${input}.R2.unpaired.fastq.gz CROP:110; done

```

<br>

* 10X Genomics datasets (CellRanger)  


  Cellranger run example  
```{r, eval=F, echo=TRUE}
# Change to the directory where your single-cell RNA-seq data is located
cd /path
# Run cellranger to process the data
cellranger run --id=sample1 \          # Set a unique run ID
              --fastqs=/path/to/fastq \ # Specify the path to the FASTQ files
              --sample=sample1 \        # Specify the sample name
              --transcriptome=/path/to/cellranger_reference \  # Path to the reference transcriptome
              --expect-cells=3000       # Estimate of the number of cells in the dataset

# The above commend includes the following steps:
# - Demultiplexing: Demultiplex and generate sample-specific FASTQ files.
# - Cellranger count: Align reads, quantify gene expression, and generate count matrices.
# - Aggregation (if multiple samples): Combine count matrices from multiple samples.
# - Analysis: Generate analysis results for downstream exploration.

# Once the analysis is complete, you can find results in the 'sample1' directory specified by --id.
```


```{r, eval=F, echo=TRUE}
# Change to the directory 
cd ~/analysis_results

# Run cellranger to process the data
cellranger run --id=my_dataset \            # Set a unique run ID
              --fastqs=/data/cellranger-tiny \ # Specify the path to the FASTQ files
              --sample=my_sample \          # Specify the sample name
              --transcriptome=/data/refdata-cellranger-GRCh38-3.0.0 \  # Path to the reference transcriptome
              --expect-cells=10000 \           # Estimated number of cells in the dataset
              --localcores=16 \                # Number of CPU cores to use
              --localmem=64 \                 # Amount of memory to use (in gigabytes)

```


  Cellranger aggr example   
```{r, eval=F, echo=TRUE}
nohup /mnt/centers/bxdx/mvb77/tools/cellranger-6.1.2/cellranger aggr --id P30279G --csv=P30279_gex_aggr.csv

```


  Cellranger vdj example  
```{r, eval=F, echo=TRUE}
# This is an example
nohup /mnt/centers/bxdx/mvb77/tools/cellranger-6.1.2/cellranger vdj --id=TE --reference=/mnt/centers/bxdx/mvb77/tools/cellrangerdbs/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0 --fastqs=/mnt/centers/bxdx/srp68/data/DFCI_Barbie/2023_Koji/raw_data/TE --sample=TE_CKDL220031043-1A_HMLG2DSX5 --localcores=4 --chain TR > TE.23.01.14.v1.out
```


    The V(D)J pipeline outputs files
    the output contains   
  
    * web_summary.html - similar to gene expression    
    * metrics_summary.csv - similar to gene expression    
    * annotation CSV/JSONs - filtered_contig_annotations.csv, clonotypes.csv    
    * FASTQ/FASTAs - filtered_contig.fasta/filtered_contig.fastq    
    * barcoded BAMs - consensus alignment mapping files    
    * cell_barcodes.json - barcodes which are identified as targeted cells.    


<br>

* STAR alignment for RNA-seq    

```{r, eval=F, echo=TRUE}
### STAR alignment for RNA-seq data
cd $path/star_out
ulimit -n 1000
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
 sampleid="${i}"
mkdir $path/${sampleid}
gzip -d $path/trimmed/${sampleid}.R1.paired.fastq.gz
gzip -d $path/trimmed/${sampleid}.R2.paired.fastq.gz
~/Desktop/software/STAR-2.7.9a/bin/MacOSX_x86_64/STAR --runThreadN 8 \
--readFilesIn $path/trimmed/${sampleid}.R1.paired.fastq $path/trimmed/${sampleid}.R2.paired.fastq   \
--genomeDir ~/Desktop/ref/mm10/star_index \
--outFileNamePrefix $path/star_out/${sampleid} \
--outSAMtype BAM SortedByCoordinate ;
done 


cd $path/star_out
### sorting star output
mkdir $path/sorted
cd $path/star_out
for i in wt1 wt2 wt3 ko1 ko2 ko3;
do
 sampleid="${i}"
 samtools sort -@ 10 $path/star_out/${sampleid} -o $path/sorted/${sampleid}.sorted.bam
 samtools index -@ 10 $path/sorted/${sampleid}.sorted.bam 
done
```

<br>

* featureCount for RNA-seq  
```{r, eval=F, echo=T}
### featureCounts
mkdir $path/features
cd $path/sorted
for i in files;
do
 sampleid="${i}"
featureCounts -a ~/Desktop/ref/mm10_refgenes.gtf  -o $path/features/${sampleid}_featurecounts.txt $path/sorted/${sampleid}.sorted.bam 2> $path/features/${sampleid}.log
done

### cleaning up the featurecounts matrix
cd $path/features
for i in files;
do
 sampleid="${i}"
cut -f 1,7,8,9,10,11,12 $path/features/${sampleid}_featurecounts.txt > $path/features/${sampleid}_featurecounts.matrix.txt 
done
```



# Frequently used Report formats   

Markdown(html)   
Rmarkdown(html)   
PPT slides(presentation)   
Excel and excel dashboard    
Word format (for SOPs, method/guideline documents)   
Notion(tracking multiple projects)    





# Machine Learning application to scRNA-seq data    



<br>
<br>

# Useful links :

Single-cell RNA-seq:   
Integration: Seurat   
https://hbctraining.github.io/scRNA-seq_online/lessons/06_integration.html    
Integration: Harmony    
https://hbctraining.github.io/scRNA-seq_online/lessons/06a_integration_harmony.html    


RNA-seq:
RNA-Seq Analysis in R using Rsubread   
https://www.alzheimersworkbench.ucsd.edu/EndToEndAnalysis_RNASeq.html    


ChIP-seq:   
Differential Peak calling using DiffBind<br>    https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html  



