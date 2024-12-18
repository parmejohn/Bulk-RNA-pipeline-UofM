params.counts_matrix = ''
params.normalized_counts = ''
params.sample_table = ''
params.outdir = ''
params.species = ''
params.bind = ''
params.chr_filter = 'none'
params.genes_gsea_filter = 'none'
params.majiq_config = 'none'
params.gff = ''
params.run_genewalk = true
params.run_neuroestimator = true


process DGEANALYSES {
    containerOptions "--bind $params.bind --no-home"
    cache 'deep'
    debug true
    label 'general_analyses'
  
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
    path "*res.txt", emit: deg_files
    path "gsea*.txt"
    path "*.pdf"
    
    script:
    """
	  ${projectDir}/src/DGEAnalyses.R --i $counts_matrix -normalized $normalized_counts -species $species -table $table -chr_filter $chr_filter -genes_gsea_filter $genes_gsea_filter
	  """
}

process GENEWALK {
    containerOptions "--bind $params.bind --no-home"
    cache 'deep'
    debug true
    label 'genewalk'
  
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
    #!/bin/bash
    
    awk 'NR > 1 && \$7 < 0.05 && \$8 != "NA" {print \$NF}' $input_file | sort -k3,3nr > genewalk_filtered.txt
    filename=$input_file
    genewalk --project results --genes genewalk_filtered.txt --id_type $id_type --base_folder "\${filename%.*}"
	  rm genewalk_filtered.txt
	  """
}

process MAJIQBUILD {
    containerOptions "--bind $params.bind --no-home"
    cache 'lenient'
    debug true
    label 'general_analyses'
  
    publishDir (
    path: "$params.outdir/majiq/",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    path input_file
    path gff
    path outdir

    output:
    path "build/*.majiq", emit: majiq_files
    path "build/*.log"
    path "build/*.sj"
    path "build/*graph.sql", emit: splicegraph
    path "build/*.sql-journal", optional: true
    val true, emit: report

    script:
    """
    #!/bin/bash
    
    majiq build -o build -c $input_file $gff
	  """
}

process MAJIQQUANT {
    containerOptions "--bind $params.bind --no-home"
    cache 'lenient'
    debug true
    label 'general_analyses'
  
    publishDir (
    path: "$params.outdir/majiq/",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    path input_file
    path outdir
    val magicbuild_signal

    output:
    path "deltapsi/*.log"
    path "deltapsi/*.tsv"
    path "deltapsi/*.voila", emit: voila_file
    
    script:
    """
    #!/bin/bash
    
    declare -A groups
    while IFS='=' read -r key value; do
        if [[ \$key != "["* ]]; then
            groups["\$key"]=\$value
        fi
    done < <(awk '/^\\[experiments\\]/,/^\$/' $input_file | tail -n +2)
    
    # Generate pairwise comparisons without duplicates
    group_names=("\${!groups[@]}")
    for ((i=0; i<\${#group_names[@]}; i++)); do
        for ((j=i+1; j<\${#group_names[@]}; j++)); do
            grp1=\${group_names[i]}
            grp2=\${group_names[j]}
    
            # Format all samples for grp1
            grp1_samples=\$(echo "\${groups[\$grp1]}" | tr ',' '\n' | sed 's|^|majiq/build/|' | sed 's|\$|.majiq|' | tr '\n' ' ')
            # Format all samples for grp2
            grp2_samples=\$(echo "\${groups[\$grp2]}" | tr ',' '\n' | sed 's|^|majiq/build/|' | sed 's|\$|.majiq|' | tr '\n' ' ')
    
            # Construct the command
            majiq deltapsi -o deltapsi -grp1 \$grp1_samples -grp2 \$grp2_samples -n \$grp1 \$grp2
        done
    done

	  """
}


process VOILA {
    containerOptions "--bind $params.bind --no-home"
    cache 'lenient'
    debug true
    label 'general_analyses'
  
    publishDir (
    path: "$params.outdir/majiq/voila_tsv",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    path voila_file
    path splicegraph

    output:
    path "*"

    script:
    """
    #!/bin/bash
    
    x=\$(basename $voila_file .voila)
    voila tsv $voila_file $splicegraph -f \${x}_voila.tsv
	  """
}

process VOILAMOD {
    containerOptions "--bind $params.bind --no-home"
    cache 'lenient'
    debug true
    label 'general_analyses'
  
    publishDir (
    path: "$params.outdir/majiq/",
    mode: 'copy',
    overwrite: true,
    pattern:'*'
    )
    
    input:
    path voila_file
    path splicegraph

    output:
    path "*"

    script:
    """
    #!/bin/bash
    
    x=\$(basename $voila_file .voila)
    voila modulize $voila_file $splicegraph -d modulized --overwrite
	  """
}

process NEUROESTIMATOR {
    containerOptions "--bind $params.bind --no-home"
    cache 'lenient'
    debug true
    label 'neuroestimator'
        
    publishDir (
      path: "$params.outdir/neuroestimator/",
      mode: 'copy',
      overwrite: true,
      pattern:'*'
    )
    
    input:
    path counts_table
    val species
    
    output:
    path "neuroestimator_results.txt"
    
    script:
    """
	  ${projectDir}/src/Neuroestimator.R $counts_table $species
    """
}

process NEUROESTIMATORPLOT {
    containerOptions "--bind $params.bind --no-home"
    cache 'lenient'
    debug true
    label 'neuroestimator_plot'
    
    publishDir (
      path: "$params.outdir/neuroestimator/",
      mode: 'copy',
      overwrite: true,
      pattern:'*'
    )
    
    input:
    path neuroestimator_results
    path sample_table
    
    output:
    path "neuroestimator_results.pdf"
    
    script:
    """
	  ${projectDir}/src/NeuroestimatorPlot.R $neuroestimator_results $sample_table
    """
}

workflow {
  // DIFFERENTIAL GENE EXPRESSION ANALYSES
  DGEANALYSES(params.counts_matrix, params.normalized_counts, params.species, params.sample_table, params.chr_filter, params.genes_gsea_filter)
	deg_files_ch = DGEANALYSES.out.deg_files.flatten()

	if (params.run_genewalk){
    if (params.species == 'homosapiens'){
      GENEWALK(deg_files_ch, 'hgnc_symbol')
    } else if (params.species == 'musmuculus') {
      GENEWALK(deg_files_ch, 'mgi_id')
    }
	}
  
  // DIFFERENTIALSPLICING
  if (params.majiq_config != 'none'){
    MAJIQBUILD(params.majiq_config, params.gff, params.outdir)
    MAJIQQUANT(params.majiq_config, params.outdir, MAJIQBUILD.out.report)

    VOILA(MAJIQQUANT.out.voila_file.flatten(), MAJIQBUILD.out.splicegraph)
    VOILAMOD(MAJIQQUANT.out.voila_file, MAJIQBUILD.out.splicegraph)
  }
  
  // NEUROESTIMATOR
  if (params.run_neuroestimator){
    neuroestimator_ch = NEUROESTIMATOR(params.counts_matrix, params.species)
    NEUROESTIMATORPLOT(neuroestimator_ch, params.sample_table)
  }
}
