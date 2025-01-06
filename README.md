# Bulk-RNA-pipeline-UofM

## Installation

`git clone https://github.com/parmejohn/Bulk-RNA-pipeline-UofM.git`

### Dependencies
- NextFlow >= 23.10.1
- Apptainer(singularity); tested with version 1.2.5-1.el7
- Highly reccommend using the [nfcore/rnaseq](https://nf-co.re/rnaseq/3.17.0/) pipeline for QC, alignment, and quantification for raw bulk RNA-seq

## Usage
Run the pipeline in the directory of your choice

```
nextflow run */Bulk-RNA-pipeline-UofM/bulk_rnaseq.nf \
	--counts_matrix <path> \
	--normalized_counts <path> \
	--sample_table <path>
	--species <homosapiens|musmuculus> \
	--bind <path> \
	--outdir <path> \
	--majiq_config [path] \
	--gff [path] \
	--chr_filter [string] \
	--genes_gsea_filter [string]
	--run_genewalk [**true**|false]
	--run_neuroestimator [**true**|false]
	...
```

### Arguments
- <> = required
- [] = optional
- ... = nextflow options, use `nextflow run -h` for all possible options

- **counts_matrix**: Tab-delimited raw count matrix. Column 1 should be labelled `gene_id` and contain Ensembl IDs, column 2 should be labelled `gene_name` and have the corresponding gene symbol. The rest of the columns should match the sample names provided in the `sample_table` file.
- **normalized_counts**: Same as the raw, except the counts should be normalized
- **sample_table**: Comma-separated file. Header is required to be `sample,condition`. See example below for more info.
- **species**: Only `homoesapiens` or `musmusculus` is supported for now.
- **bind**: Paths where input or output data is moving. This is because apptainer(singularity) without root permission needs explicit file paths. Multiple paths can be defined with a ',' to seperate them. 
	- For instance, when in a file system where you cannot access the root, input files may be located in `/research/` but your config and reference files are in `/home/USER/`. This would look like `--bind /research/,/home/USER/`
- **outdir**: Folder where results will be outputted
- **majiq_config**: Required information for running MAJIQ. This settings config that needs to match MAJIQ's specifications as seen from the example below. More info on their [documentation](https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/builder.html). If this parameter is left empty, MAJIQ will be skipped
	- Will require BAM files and preferably knowing the strandness of your samples.
		- Strandness can be determined from the nfcore/rnaseq multiqc report or other methods
- **gff**: A GFF3 file that will be used in MAJIQ.
- **chr_filter**: Chromosomes to filter out during differential gene expression analyses. e.g. `--chr_filter X,1,2,3`
- **genes_gsea_filter**: Genes to filter out during gene-set enrichment analysis (GSEA). Can be applied when there are knocked out genes. e.g. `--genes_gsea_filter FAN1,PRDM1`
- **run_genewalk**: Perform [Genewalk](https://churchman.med.harvard.edu/genewalk) analysis. Determines relevant functions for individual genes given different conditions.
- **run_neuroestimator**: Perform [NEUROeSTIMator](https://research-git.uiowa.edu/michaelson-lab-public/neuroestimator) analysis. Quantifies neuronal activity from RNA-seq data.

#### counts_matrix/normalized_matrix example
```
gene_id gene_name       npc-1-A9        npc-1-C8        npc-1-D9        npc-1-F8	...
ENSG00000000003 TSPAN6  6000.92038      8579.833        6068.5994       8649.32162      ...
ENSG00000000005 TNMD    2       4       61      2       ...
```

#### sample_table example
```
sample,condition
npc-1-A9,RTT
npc-1-C8,RTTFAN1KO
npc-1-D9,CTRL
npc-1-F8,CTRLFAN1KO
```
#### majiq_config example
```
[info]
bamdirs=/nfcore_rnaseq_output/input/star_salmon/
genome=hg38
strandness=reverse

[experiments]
CTRL=npc-1-D9.markdup.sorted,...
CTRLFAN1KO=npc-1-F8.markdup.sorted,...
RTT=npc-1-A9.markdup.sorted,...
RTTFAN1KO=npc-1-C8.markdup.sorted,...
```
One can also use MAJIQ's online [Command line builder](https://biociphers.bitbucket.io/majiq-docs-academic/commandbuilder.html#command-builder) to generate the this file.

### Outputs
#### Summary
- [Differential gene expression](#differential-gene-expression)
- [Genewalk](#genewalk)
- [Differential splicing (MAJIQ)](#differential-splicing)
	- [Column interpretations](#column-interpretations)
- [Quantifying neuronal activity (NEUROeSTIMator)](#quantifying-neuronal-activity)

#### Differential gene expression
<details>
<summary>Click to expand</summary>
<br>

- Perform differential gene expression analysis using DESeq2
- Using the log2FC values as rank, performs GSEA using the fgsea package
	- All genes are used, as a opposed to an over-representation analyses, which would only use DEGs after a certain cutoff adjusted p-value and/or log2FC cutoff.
	- To find specific DEGs and their related pathways, one will have to filter the list by their leading edge in the resulting data files.
- Files = dge/
  - Plots
    - **deg_heatmap_X_vs_Y.pdf**: Heatmap of DEGs for a given comparison. Rows are hierarchically clustered and expression is Z-scaled. |log2FC| >= 2 & padj < 0.05
    - **deseq2_volcano_X_vs_Y.pdf**: Volcano plot of all genes for a given comparison. Please note that cutoff uses padj but y-axis uses the unadjusted p-value. |log2FC| >= 2 & padj < 0.05
    - **gsea_X_vs_Y.pdf**: GSEA for a given comparison. Positive (red) normalized enrichment score (NES)  are upregulated in X; Negative NES are upregulated in Y
    - **sample_similarity_heatmap.pdf**: Heatmap representing sample clustering using a distance matrix. Counts are transformed via variance stabilizing transformation (VST) in DESeq2. Hierarchical is done on both the rows and columns.
    	- Notebly, large portion of genes will not be differentially expressed, so it can lead to some samples having high similarity scores to one another.
    - **sample_similarity_pca.pdf**: PCA plot of all samples. Grouping is done by condition by default.
    - **deg_upset_upreg.pdf**: Upset plot showing intersection of all upregulated genes (Genes upregulated in condition X when compared to Y). log2FC >= 2 & padj < 0.05
    - **deg_upset_dnreg.pdf**: Upset plot showing intersection of all downregulated genes (Genes upregulated in condition Y when compared to X). log2FC <= -2 & padj < 0.05
  - Data
    - **deseq2_X_vs_Y_res.txt**: Tab-delimited file with DESeq2 results, unfiltered
    - **gsea_X_vs_Y.txt**: Tab-delimited file with GSEA results, unfiltered
</details>


#### Genewalk
<details>
<summary>Click to expand</summary>
<br>

- Identifies relevant functions for individual genes
	- Uses the DESeq2 results (ranked by log2FC)
	- Determines the importance of the gene and what biological pathways it impacts
- Files = genewalk/
  - **deseq2_X_vs_Y_res/:** Please refer to their thorough [documentation](https://churchman.med.harvard.edu/genewalk) and [GitHub](https://github.com/churchmanlab/genewalk) for file descriptions.
</details>

#### Differential splicing (MAJIQ)
<details>
<summary>Click to expand</summary>
<br>

- 3 modules as outlined [here](https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/quick-overview.html)
	- Builder: Define splice graphs and local splicing variations (LSVs)
		- Creates a splicegraph from the gene annotation and aligned reads
			- This will contain junctions, intron retention scores, and exon information
	- Quantifier: Quantifies the relative abundance (Ψ) of LSVs and changes in the relative abundance (ΔΨ) between conditions
	- Voila: Create visualizations and interpretable files
- Aids in discovering differential splice variants
- Files = majiq/
  - **build/**: Files to use downstream MAJIQ functions
  - **voila_tsv/**: .tsv files produced by majiq quant to be use in downstream voila visualizations. The .tsv's contain the relative LSV abundance (Ψ) and changes in the relative LSV abundance (delta PSI) between conditions.
  - **modulized/**:
  	- **summary.tsv**: Each row is a splicing module and list the total counts for each type of splicing event in the module
		- A module is defined as single entry and exit regions of the splicegraph
			- Unique single source and single target exon
	- **heatmap.tsv**: When using ΔΨ, the junctions with the max absolute ΔΨ value is chosen from each module, and one can see all of the deltapsi values for each comparison listed here.
    - ... Other files that can be parsed for additional information

- For more in-depth information, please visit their [website](https://majiq.biociphers.org/)
</details>

##### Column interpretations
<details>
<summary>Click to expand</summary>
<br>

- Below is a description of each column that is shown in MAJIQ. Please visit the MAJIQ [documentation](https://biociphers.bitbucket.io/majiq-docs-academic/getting-started-guide/quick-overview.html) or their [BitBucket](https://bitbucket.org/biociphers/majiq_academic) page
	- **lsv_id** = `s` or `t` denotes whether it is the source or target exon
		- can be visualized when looking using voila view (see below on how to run it)
		- [Example](https://majiq.biociphers.org/green_et_al_2017/examples/hogenesch/adr-cer-8v8/) from the MAJIQ documentation
	- **lsv_type** = rough graphical output for the voila view plot
		- From their [Google Groups forums](https://groups.google.com/g/majiq_voila/c/tOrbP179tuY)
			- starts by (s or t) being source or target
			- each '|' is a new junction representation and if there is intron_retention the last character is 'i'
			- each junction is represented by  XeY.ZoK where 
				- X is the ordinal splice site in the reference exon
				- Y is the ordinal exon connecting the lsv
				- Z is the ordinal splice site in exon Y
				- K is the total number of splice sites that Y has
	- **mean_dpsi_per_lsv_junction** = direction of change for a given junction
		- One can think of it as the the fold-change equivelent in a differential expression analysis
	- **probability_changing** = probability, that the dpsi is above <threshold 1 used>
		- default is 0.2
	- **probability_non_changing** = probability, that the dpsi is below <threshold 2 used>
		- default = 0.05
	- **\*mean_psi** = E |PSI| adds up to 100% of all of the LSV's junction
		- mean psi value from the different experiments
		- **negative** values correspond to **increased differential inclusion in condition1** compared with condition2
	- num_junctions
	- num_exons
	- **de_novo_junctions** =  which of these junctions are unannotated
		- junctions are annotated using the GFF3 file and RNA-seq files
	- seqid = unsure
	- strand
	- **junctions_coords** = the positions of the gene denoting the exon junction locations and how large they are
	- **exons_coords** = location and size of a given exon
	- **ir_coords** = intron retention coordinates
	- **ucsc_lsv_link** = genome browser view of the whole lsv region
</details>

#### Quantifying neuronal activity (NEUROeSTIMator)
<details>
<summary>Click to expand</summary>
<br>

- Estimates neuronal activation using gene expression
	- Uses a neural network approach along with 
	- Scores based off 22 neuronal activity markers
		- Each sample will have an activity score that ranges from 0-1
- Files = neuroestimator/
  - **neuroestimator_results.txt**: Tab-delimited table with samples and their predicted activity score
  - **neuroestimator_results.pdf**: Box plot for predicted activities across conditions. Kolmogorov–Smirnov (KS) test is performed to test for significance between the each condition.

- For more in-depth file descriptions please visit
</details>

### Using viola view

Below I outline how to run viola view for a *deltapsi.viola file on a Linux server using apptainer

1. Download the apptainer image `apptainer pull docker://parmejohn/bulk_rna:1.0.5`
2. Run an interactive terminal session using with the correct mounts `apptainer shell </path/to/file/>bulk_rna_1.0.5.sif -b </path/to/file/>`
    - If you have multiple paths, please separate them with commas (no whitespace)
4. Run `viola view` on a port of your choosing (note that only 1 deltapsi.voila can be viewed at a time)
    - `voila view -p <PORT> splicegraph.sql *.deltapsi.voila`
5. Connect to the server port on your local terminal `ssh -N -L <PORT>:localhost:<PORT> <USER@server.address>`
Open a browser and go to `http://localhost:<PORT>/`
