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
    val report_name

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
    debug true
  
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
    path "${sample_id}/"
/*    path "*.tsv"
    path "*.json"*/

    script:
    """
	  kallisto quant -i $kallisto_idx -b 100 -t 16 -o ${sample_id} ${reads} 2>&1 | tee ${sample_id}_kallisto.log
	  """
}
//	  mv ${sample_id}_kallisto.log ${sample_id}/${sample_id}_kallisto.log

workflow {
  pairs_ch = Channel
    .fromFilePairs(params.indir.toString() + '*_{R1,R2}.fastq.gz')

  // QC is done before and after trimming
  fastqc_ch = FASTQC(pairs_ch)

  // Quantification and QC
  index_ch = INDEX(params.transcriptome_file)
  
  TRIM(pairs_ch)
  trimmed_pairs_ch = TRIM.out.trimmed_fq
  trimmed_pairs_ch.view()
  
  KALLISTO(index_ch, trimmed_pairs_ch)
  
  //MULTIQC(KALLISTO.out.logs_ch.mix(TRIM.out.fastqc_ch).collect())
  
  // Analyses
/*  SLEUTH
  GSEA
  GENEWALK*/
}