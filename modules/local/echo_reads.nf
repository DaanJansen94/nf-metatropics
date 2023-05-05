process ECHO_READS {

    debug true

    input:
    tuple val(meta), path(reads)

    script:
    """
    echo ${reads}
    """
}
