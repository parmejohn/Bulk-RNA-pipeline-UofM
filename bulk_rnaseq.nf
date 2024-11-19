params.counts_matrix = ''
params.normalized_counts = ''
params.sample_table = ''
params.outdir = ''
params.species = ''
params.bind = ''
params.chr_filter = 'none'
params.genes_gsea_filter = 'none'


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
    path counts_matrix
    path normalized_counts
    val species
    path table
    val chr_filter
    val genes_gsea_filter

    output:
    path "*genewalk.txt", emit: deg_files
    path "*res.txt", optional: true
    path "gsea*.txt"
    path "*.pdf"
    
    script:
    """
	  ${projectDir}/src/DGEAnalyses.R --i $counts_matrix -normalized $normalized_counts -species $species -table $table -chr_filter $chr_filter -genes_gsea_filter $genes_gsea_filter
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
    awk 'NR > 1 && \$7 < 0.05 && \$8 != "NA" {print \$NF}' $input_file | sort -k3,3nr > genewalk_filtered.txt
    filename=$input_file
    genewalk --project results --genes genewalk_filtered.txt --id_type $id_type --base_folder "\${filename%.*}"
	  rm genewalk_filtered.txt
	  """
}

workflow {
  // Analyses
  DGEANALYSES(params.counts_matrix, params.normalized_counts, 'homosapiens', params.sample_table, params.chr_filter, params.genes_gsea_filter)

	deg_files_ch = DGEANALYSES.out.deg_files.flatten()
  if (params.species == 'homosapiens'){
    GENEWALK(deg_files_ch, 'hgnc_symbol')
  } else if (params.species == 'musmuculus') {
    GENEWALK(deg_files_ch, 'mgi_id')
  }
  //DIFFERENTIALSPLICING
}
