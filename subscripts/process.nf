process fastp {
    maxForks 6
    conda 'bioconda::fastp=0.23.4'
    publishDir "${launchDir}/${params.output}/fastp/reports", pattern: "{*html,*json}", mode: "copy"
    publishDir "${launchDir}/${params.output}/fastp", pattern: "*gz", mode: "copy"
    input:
    tuple val(f_name), path(fa)
    output:
    tuple val(f_name), path(fa)
	tuple val(f_name), path("${f_name}{1,2}_trimmed.fq.gz"), emit: trimmed
	path "*{json,html}", emit: reports
    script:
    """
    fastp \
    -i ${fa[0]} \
    -o ${f_name}1_trimmed.fq.gz \
    -I ${fa[1]} \
    -O ${f_name}2_trimmed.fq.gz \
    -z $params.compress_level \
    -V \
    -g \
    --poly_g_min_len $params.poly_g_min_len \
    -x \
    --poly_x_min_len $params.poly_x_min_len \
    -5 \
    -3 \
    -M $params.cut_mean_quality \
    -n $params.n_base_limit \
    -e $params.average_qual \
    -l $params.length_required \
    -c \
    -w $params.fastp_threads \
    -j ${f_name}_fastp.json \
    -h ${f_name}_fastp.html
    """
}


process fastqc {
    conda 'bioconda::fastqc=0.12.1'
    publishDir "${launchDir}/${params.output}/${fqc_dir}", pattern: "*html", mode: "copy"
    publishDir "${launchDir}/${params.output}/${fqc_dir}/zip", pattern: "*zip", mode: "copy"
    input:
    tuple val(index), path(reads)
    val fqc_dir
    output:
    tuple val(index), path(reads)
    path "*", emit: reports
    script:
    """
    fastqc -t 6 *
    """
}


process multiqc {
    conda "${projectDir}/multiqc"
    publishDir "${launchDir}/${params.output}/${mqc_dir}", mode: "copy"
    input:
    path(reports)
    val mqc_dir
    output:
    path "*"
    script:
    """
    multiqc --config ${projectDir}/multiqc.yaml *
    """
}


process bwa {
    maxForks 4
    conda "bioconda::bwa==0.7.18 bioconda::samtools=1.18 conda-forge::ncurses"
    publishDir "${launchDir}/${params.output}/${map_dir}/", pattern: "unmapped_reads*", mode: "copy"
    input:
    tuple val(name), path('read?.fq.gz')
    val ref
    val map_dir
    output:
    tuple val(name), path("unmapped_reads_R{1,2}.fq.gz"), emit: unmapped_tuple
    path "*"
    script:
    """
    bash ${projectDir}/subscripts/bwa_indexer.sh $ref
    bwa mem -t ${params.map_cpus} $ref read1.fq.gz read2.fq.gz | samtools view -Sb > raw.bam
    samtools sort -@ ${params.map_cpus} raw.bam -o sorted.bam
    samtools view -f 12 sorted.bam -o unmapped.bam
    samtools fastq unmapped.bam -1 unmapped_reads_R1.fq -2 unmapped_reads_R2.fq
    gzip unmapped_reads*.fq
    """
}


process kraken2 {
	conda = 'bioconda::kraken2'
	maxForks 1
	publishDir "${launchDir}/${params.output}/kraken2", pattern: "*kreport", mode: "copy"
	publishDir "${launchDir}/${params.output}/kraken2/unclassified_reads", pattern: "*unclassified*", mode: "copy"
	input:
    tuple val(index), path(reads)
	output:
    tuple val(index), path("*.kreport"), emit: stats
	path "*"
	script:
	"""
	kraken2 --db ${projectDir}/kraken2_data --report ${index}.kreport --gzip-compressed --paired --unclassified-out ${index}_unclassified_#.fq ${reads[0]} ${reads[1]}
	"""
}


process bracken {
    conda = 'bioconda::bracken'
    publishDir "${launchDir}/${params.output}/bracken", pattern: "*bracken", mode: 'copy'
    input:
    tuple val(index), path(kreport)
    output:
    tuple val(index), path("*.bracken"), emit: stats
    path "*"
    script:
    """
    bracken -d $params.kraken_db -i $kreport -o ${index}.bracken -l $params.classification_lvl
    """
}


process metaphlan { // metaphlan --install --bowtie2db <database folder> needed
    maxForks 1
    conda "bioconda::metaphlan"
    publishDir "${launchDir}/${params.output}/metaphlan", mode: "copy"
    input:
    tuple val(name),path(reads)
    output:
    path "*.txt", emit: metaphlan_output
    path "*"
    script:
    """
    metaphlan \
    ${reads[0]},${reads[1]} \
    --input_type fastq \
    --bowtie2db ${projectDir}/metaphlan \
    --bowtie2out  ${name}.bz2 \
    -o ${name}.txt \
    -t $params.metaphlan_analysis_type \
    --nproc $params.metaphlan_cpus
    """

}


process metaphlan_merge {
    maxForks 1
    conda "bioconda::metaphlan=4.06"
    publishDir "${launchDir}/${params.output}/metaphlan", mode: "copy"
    input:
    path metaphlan_output
    output:
    path "merged_abundance_table.txt"
    script:
    """
    merge_metaphlan_tables.py *.txt > merged_abundance_table.txt
    """
}


process megahit {
    maxForks 1
    conda "bioconda::megahit=1.2.9"
    publishDir "${launchDir}/${params.output}", mode: "copy"
    input:
    path reads1
    path reads2
    output:
    path "megahit/final.contigs.fa", emit: assembly
    path "megahit"
    script:
    reads1 = reads1.join(',')
    reads2 = reads2.join(',')
    """
    megahit \
    -1 ${reads1} \
    -2 ${reads2} \
    -o megahit \
    """
}


process metawrapbin {
    conda "${projectDir}/mw-env" //"ursky::metawrap-mg=1.3.2"
    maxForks 1
    publishDir "${launchDir}/${params.output}", mode: "copy"
    input:
    path assembly
    path reads1, stageAs: 'read*_1.fastq'
    path reads2, stageAs: 'read*_2.fastq'
    output:
    path "metawrap_bins/*"
    path "metawrap_bins/metabat2_bins/bin.*.fa", emit: metabat2_bins, optional: true
    path "metawrap_bins/concoct_bins/bin.*.fa", emit: concoct_bins, optional: true
    path "metawrap_bins/maxbin2_bins/bin.*.fa", emit: maxbin2_bins, optional: true
    script:
    binners = ""
    if (params.concoct) {
        binners += "--concoct "
    }
    if (params.maxbin2) {
        binners += "--maxbin2 "
    }
    if (params.metabat2) {
        binners += "--metabat2 "
    }
    """
    echo ${projectDir}/checkm_data | checkm data setRoot
    metawrap binning -t ${params.mw_cpus} \
    $binners\
    --run-checkm \
    -a $assembly \
    -o metawrap_bins \
    $reads1 $reads2
    """
}


process metawrap_binrefinement {
    conda "${projectDir}/mw-env"
    maxForks 1
    errorStrategy "ignore"
    publishDir "${launchDir}/${params.output}", mode: "copy"
    input:
    path metabat2_bins, stageAs: "metabat2/*"
    path concoct_bins, stageAs: "concoct/*"
    path maxbin2_bins, stageAs: "maxbin2/*"
    output:
    path "metawrap_binref/*"
    path "metawrap_binref/metawrap*bins/*fa", emit: refined_bins
    script:
    not_nulls = []
    bin_dirs = ''
    if (!(metabat2_bins ==~ /.*null.*/)) {
        not_nulls.add('metabat2')
    }
    if (!(concoct_bins ==~ /.*null.*/)) {
        not_nulls.add('concoct')
    }
    if (!(maxbin2_bins ==~ /.*null.*/)) {
        not_nulls.add('maxbin2')
    }
    if (not_nulls.size == 3) {
        bin_dirs += "-A ${not_nulls[0]} -B ${not_nulls[1]} -C ${not_nulls[2]}"
    }
    if (not_nulls.size == 2) {
        bin_dirs += "-A ${not_nulls[0]} -B ${not_nulls[1]}"
    }
    if (not_nulls.size == 1) {
        bin_dirs += "-A ${not_nulls[0]}"
    }
    """
    metawrap bin_refinement \
    -o metawrap_binref \
    -t ${params.mw_cpus} \
    -m 100 \
    -c $params.metawrap_completion -x $params.metawrap_contamination \
    $bin_dirs
    """
}


process gunzip {
    input:
    path "R1_?.fq.gz"
    path "R2_?.fq.gz"
    output:
    path "R1*.fq", emit: reads1
    path "R2*.fq", emit: reads2
    script:
    """
    gunzip -k --force R*.gz
    """
}


process gtdbtk_classify {
    conda "bioconda::gtdbtk=2.3.2"
    maxForks 1
    publishDir "${launchDir}/${params.output}", mode: "copy"
    input:
    path bin, stageAs: "bins/*"
    output:
    stdout
    path "gtdbtk_classify/*"
    script:
    """
    export GTDBTK_DATA_PATH=${projectDir}/gtdbtk_data
    mkdir scratch
    gtdbtk classify_wf \
    --mash_db ${workflow.launchDir}/${params.output}/mash_db \
    --genome_dir bins/ \
    --extension fa \
    --out_dir gtdbtk_classify \
    --scratch_dir scratch \
    --cpus 20
    """
}


process merge_kr2_br {
    conda "anaconda::pandas conda-forge::numpy"
    publishDir "${launchDir}/${params.output}/kraken2", pattern: "*csv", mode: "copy"
    input:
    tuple val(index1), path(kraken2)
    tuple val(index2), path(bracken)
    output:
    path "*"
    script:
    """
    python3 ${projectDir}/subscripts/merge_kr2br.py -k . -b .
    """
}
