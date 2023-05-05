process FIX_NAMES {

    tag{sample}

    input:
    tuple val(meta), val(sample), path(reads)

    output:
    tuple val(sample), path("*.fastq"), emit : fqreads

    script:
    """
    cat $reads/* > ${sample}.fastq
    """
}
