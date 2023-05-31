//
// Check input samplesheet and get read channels
//

include { MINIMAP2_ALIGN              } from '../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT               } from '../../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX              } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FASTQ              } from '../../modules/nf-core/samtools/fastq/main'

workflow HUMAN_MAPPING {
    take:
    readsONT // file: /path/to/samplesheet.csv

    main:
    MINIMAP2_ALIGN(
        readsONT,
        params.fasta,
        true,
        false,
        false
    )
    SAMTOOLS_SORT(
        MINIMAP2_ALIGN.out.bam
    )
    //SAMTOOLS_SORT.out.bam.view()

    SAMTOOLS_FASTQ(
        SAMTOOLS_SORT.out.bam,
        false
    )

    emit:
    nohumanreads = SAMTOOLS_FASTQ.out.other                                     // channel: [ val(meta), [ reads ] ]
    versionsmini = MINIMAP2_ALIGN.out.versions // channel: [ versions.yml ]
    versionssamsort = SAMTOOLS_SORT.out.versions // channel: [ versions.yml ]
    versionssamfastq = SAMTOOLS_FASTQ.out.versions // channel: [ versions.yml ]
}
