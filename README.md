# Bulk-RNA-pipeline-UofM

## Installation
### Dependencies


### Scripts


## Usage

### Arguements

### Outputs
#### Differential gene expression

#### Genewalk

#### Differential splicing (MAJIQ)

#### Quantifying neuronal activity (NEUROeSTIMator)

### Example

### Using viola view

Below I outline how to run viola view for a *deltapsi.viola file on a Linux server using apptainer

1. Download the apptainer image `apptainer pull docker://parmejohn/bulk_rna:1.0.5`
2. Run an interactive terminal session using with the correct mounts `apptainer shell /path/to/file/bulk_rna_1.0.5.sif -b PATHTOFILEDIRECTORY`
  - If you have multiple paths, please separate them with commas (no whitespace)
4. Run `viola view` on a port of your choosing (note that only 1 deltapsi.voila can be viewed at a time)
`voila view -p PORT splicegraph.sql *.deltapsi.voila`
5. Connect to the server port on your local terminal `ssh -N -L PORT:localhost:PORT USER@server.address`
Open a browser and go to `http://localhost:PORT/`
