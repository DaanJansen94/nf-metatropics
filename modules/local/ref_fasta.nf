process REF_FASTA {
    tag "$meta.id"

    container "$projectDir/images/samtools_minimap2.sif"

    input:
    tuple val(meta), path(report), path(emreads), path(rawfastq)

    output:
    tuple val(meta), path("*.fasta"), emit : seqref
    tuple val(meta), path("*.reads"), emit : headereads
    tuple val(meta), path("*.fastq"), emit : allreads
    //tuple val(meta), stdout, emit : virusout

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    produce_fasta.pl $report $emreads $rawfastq $args
    """
}
