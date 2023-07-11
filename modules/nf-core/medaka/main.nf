process MEDAKA {
    tag "${meta.id}.${meta.virus}"
    label 'process_single'

    conda "bioconda::medaka=1.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0' :
        'biocontainers/medaka:1.4.4--py38h130def0_0' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*.vcf"), emit: assembly
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bamfiles
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.virus}"
    """
    medaka_haploid_variant \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -r $assembly \\
        -o ./

    mv medaka.annotated.vcf ${prefix}.vcf
    mv calls_to_ref.bam ${prefix}.sorted.bam
    mv calls_to_ref.bam.bai ${prefix}.sorted.bam.bai
    rm medaka.vcf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
