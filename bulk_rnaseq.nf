params.indir = ''
params.outdir = ''
params.transcriptome_file = ''
params.species = ''
params.bind = ''

// code from https://training.nextflow.io/basic_training/rnaseq_pipeline/
process FASTQC {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug false
    
    publishDir (
    path: "$params.outdir/untrimmed/fastqc/",
    mode: 'copy',
    overwrite: true,
    pattern: "fastqc_*"
    )
  
    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

// code from https://training.nextflow.io/basic_training/rnaseq_pipeline/
process MULTIQC {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug false
    
    publishDir (
    path: "$params.outdir/",
    mode: 'copy',
    overwrite: true
    )
    
    input:
    path '*'

    output:
    path '*'

    script:
    """
    multiqc .
    """
}

process INDEX {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug false
    
    publishDir (
    path: "$params.outdir/",
    mode: 'copy',
    overwrite: true
    )
    
    input:
    path transcriptome

    output:
    path 'transcripts.idx'

    script:
    """
    kallisto index -i transcripts.idx $transcriptome
    """
}

process TRIM {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug true
    
    publishDir (
    path: "$params.outdir/trimmed/",
    mode: 'copy',
    overwrite: true,
    pattern:'*.gz'
    )
    
    publishDir (
    path: "$params.outdir/trimmed/",
    mode: 'copy',
    overwrite: true,
    pattern:'*.txt'
    )
    
    publishDir (
    path: "$params.outdir/trimmed/fastqc/",
    mode: 'copy',
    overwrite: true,
    pattern:'*fastqc*'
    )
  
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*{R1,R2}*val*.gz'), emit: trimmed_fq
    path '*fastqc*', emit: fastqc
    //path '*trimmed*'
    path '*.txt'

    script:
    """
    trim_galore --fastqc -q 20 --paired ${reads}
    """
}

process KALLISTO {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug false
  
    publishDir (
    path: "$params.outdir/kallisto/",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    path kallisto_idx
    tuple val(sample_id), path(reads)

    output:
    path "*.log", emit: logs_ch
    path "${sample_id}/", emit: quant
/*    path "*.tsv"
    path "*.json"*/

    script:
    """
    total_threads=\$(nproc)
    num_threads=\$((total_threads * 90 / 100))
    
	  kallisto quant -i $kallisto_idx -b 100 -t "\$num_threads" -o ${sample_id} ${reads} 2>&1 | tee ${sample_id}_kallisto.log
	  """
}

process DGEANALYSES {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug true
  
    publishDir (
    path: "$params.outdir/dge/",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    path kallisto_outputs
    val species

    output:
    path "*genewalk.txt", emit: deg_files
    path "*res.txt", optional: true
    path "gsea*.txt"
    path "*.pdf"

    script:
    """
	  ${projectDir}/src/DGEAnalyses.R --i $kallisto_outputs -species $species
	  """
}

process GENEWALK {
    containerOptions "--bind $params.bind"
    cache 'deep'
    debug true
  
    publishDir (
    path: "$params.outdir/genewalk/",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    file input_file
    val id_type

    output:
    path "*"

    script:
    """
    awk 'NR > 1 && \$7 < 0.05 && \$8 != "NA" {print \$8}' $input_file | sort -k3,3nr > genewalk_filtered.txt
    filename=$input_file
    genewalk --project results --genes genewalk_filtered.txt --id_type $id_type --base_folder "\${filename%.*}"
	  rm genewalk_filtered.txt
	  """
}

//basename $input_file | genewalk --project results --genes genewalk_filtered.txt --id_type $id_type --base_folder -
workflow {
  pairs_ch = Channel
    .fromFilePairs(params.indir.toString() + '*_{R1,R2}.fastq.gz')

  // QC is done before and after trimming
  fastqc_ch = FASTQC(pairs_ch)

  // Quantification and QC
  index_ch = INDEX(params.transcriptome_file)
  
  TRIM(pairs_ch)
  trimmed_pairs_ch = TRIM.out.trimmed_fq
  //trimmed_pairs_ch.view()
  
  KALLISTO(index_ch, trimmed_pairs_ch)
  
  untrimmed_qc_ch = fastqc_ch.collect()
  trimmed_qc_ch = TRIM.out.fastqc.collect()
  kallisto_qc_ch = KALLISTO.out.logs_ch.collect()
  
  qc_ch = untrimmed_qc_ch.concat(trimmed_qc_ch, kallisto_qc_ch).collect()
  MULTIQC(qc_ch)

  
  // Analyses
  DGEANALYSES(KALLISTO.out.quant.collect(), 'homosapiens')
  //deg_files_ch = Channel.fromPath(DGEANALYSES.out.deg_files)

	deg_files_ch = DGEANALYSES.out.deg_files.flatten()
  if (params.species == 'homosapiens'){
    GENEWALK(deg_files_ch, 'hgnc_symbol')
  } else if (params.species == 'musmuculus') {
    GENEWALK(deg_files_ch, 'mgi_id')
  }
  //DIFFERENTIALSPLICING
}
