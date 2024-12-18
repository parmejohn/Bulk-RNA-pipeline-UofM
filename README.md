# Bulk-RNA-pipeline-UofM

## Installation
### Dependencies


### Scripts


## Usage

### Arguements

### Outputs
#### Differential gene expression
<details>
<summary>Click to expand</summary>
<br>

- Normalize, find variable genes, and scale the data
- Files = dge/
  - Plots
    - deg_heatmap_X_vs_Y.pdf
    - deseq2_volcano_RTTFAN1KO_vs_CTRL.pdf
    - gsea_heatmap_X_vs_Y.pdf
    - sample_similarity_heatmap.pdf
    - sample_similarity_pca.pdf
    - deg_upset_upreg.pdf
    - deg_upset_dnreg.pdf
  - Data
    - deseq2_X_vs_Y_res.txt
    - gsea_X_vs_Y.txt
</details>


#### Genewalk
<details>
<summary>Click to expand</summary>
<br>

- Normalize, find variable genes, and scale the data
- Files = genewalk/
  - deseq2_X_vs_Y_res/ = 
</details>

#### Differential splicing (MAJIQ)
<details>
<summary>Click to expand</summary>
<br>

- Normalize, find variable genes, and scale the data
- Files = majiq/
  - build/
    - *.log
    - *.majiq
    - *.sj
    - splicegraph.sql
  - voila_tsv/
    - *deltapsi_voila.tsv
  - modulized/ (only listing what I think is the most important to view)
    - heatmap.tsv
    - summary.tsv
    - ...

- For more in-depth file descriptions please visit
</details>

#### Quantifying neuronal activity (NEUROeSTIMator)
<details>
<summary>Click to expand</summary>
<br>

- Normalize, find variable genes, and scale the data
- Files = neuroestimator/
  - neuroestimator_results.txt
  - neuroestimator_results.pdf

- For more in-depth file descriptions please visit
</details>

### Example

### Using viola view

Below I outline how to run viola view for a *deltapsi.viola file on a Linux server using apptainer

1. Download the apptainer image `apptainer pull docker://parmejohn/bulk_rna:1.0.5`
2. Run an interactive terminal session using with the correct mounts `apptainer shell /path/to/file/bulk_rna_1.0.5.sif -b PATHTOFILEDIRECTORY`
    - If you have multiple paths, please separate them with commas (no whitespace)
4. Run `viola view` on a port of your choosing (note that only 1 deltapsi.voila can be viewed at a time)
    - `voila view -p PORT splicegraph.sql *.deltapsi.voila`
5. Connect to the server port on your local terminal `ssh -N -L PORT:localhost:PORT USER@server.address`
Open a browser and go to `http://localhost:PORT/`
