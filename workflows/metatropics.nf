/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatropics.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { INPUT_CHECK_METATROPICS } from '../subworkflows/local/input_check_metatropics'
include { FIX } from '../subworkflows/local/subfix_names'
include { HUMAN_MAPPING } from '../subworkflows/local/human_mapping'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { ECHO_READS                  } from '../modules/local/echo_reads'
include { GUPPY_ONT                   } from '../modules/local/guppy/ont'
include { GUPPYDEMULTI_DEMULTIPLEXING } from '../modules/local/guppydemulti/demultiplexing'
include { FASTP                       } from '../modules/nf-core/fastp/main'
include { NANOPLOT                    } from '../modules/nf-core/nanoplot/main'
include { METAMAPS_MAP                } from '../modules/local/metamaps/map'
include { METAMAPS_CLASSIFY           } from '../modules/local/metamaps/classify'
include { R_METAPLOT                  } from '../modules/local/r/metaplot'
include { KRONA_KRONADB               } from '../modules/nf-core/krona/kronadb/main'
include { KRONA_KTIMPORTTAXONOMY      } from '../modules/nf-core/krona/ktimporttaxonomy/main'
include { REF_FASTA                   } from '../modules/local/ref_fasta'
include { SEQTK_SUBSEQ                } from '../modules/nf-core/seqtk/subseq/main'
include { REFFIX_FASTA                } from '../modules/local/reffix_fasta'
include { MEDAKA                      } from '../modules/nf-core/medaka/main'
include { SAMTOOLS_COVERAGE           } from '../modules/nf-core/samtools/coverage/main'
include { IVAR_CONSENSUS              } from '../modules/nf-core/ivar/consensus/main'
include { HOMOPOLISH_POLISHING        } from '../modules/local/homopolish/polishing'
include { ADDING_DEPTH                } from '../modules/local/adding_depth'
include { FINAL_REPORT                } from '../modules/local/final_report'
include { BAM_READCOUNT               } from '../modules/local/bam/readcount'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Info required for completion email and summary
def multiqc_report = []

workflow METATROPICS {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    //INPUT_CHECK.out.reads.map{it[1]}.view()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Echo reads
    //
    //ECHO_READS (
    //    INPUT_CHECK.out.reads
    //)

    //
    //MODULES_DEVELOPED_BY_ANTONIO
    //
    ch_input2 = params.input_fastq
    INPUT_CHECK_METATROPICS{
        ch_input2
    }
    //INPUT_CHECK_METATROPICS.out.reads.view()
    //ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1],it[0])}
    //ch_sample.view()

    if(params.basecall==true){
        ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1],it[0])}

        inFast5 = channel.fromPath(params.input_dir)
        //inFast5.view()
        GUPPY_ONT(
            inFast5
        )
        //ch_versions = ch_versions.mix(GUPPY_ONT.out.versions)

        GUPPYDEMULTI_DEMULTIPLEXING(
            GUPPY_ONT.out.basecalling_ch
        )
        //ch_versions = ch_versions.mix(GUPPYDEMULTI_DEMULTIPLEXING.out.versions)

        ch_barcode = GUPPYDEMULTI_DEMULTIPLEXING.out.barcodeReads.flatten().map{file -> tuple(file.simpleName, file)}
        ch_sample_barcode = ch_sample.join(ch_barcode)
        //ch_sample_barcode.view()

        //FIX_NAMES(
        FIX(
            ch_sample_barcode
        )
        //FIX.out.reads.view()
        //FIX.out.reads.map{it[1]}.view()
    }
    else if(params.basecall==false){
        ch_sample = INPUT_CHECK_METATROPICS.out.reads.map{tuple(it[1].replaceFirst(/\/.+\//,""),it[0],it[1])}
        //ch_sample.map{tuple(it[1].replaceFirst(/.fastq/,""),it[0],it[1])}
        //ch_sample.view()
        FIX(
            ch_sample
        )
        //FIX.out.reads.view()
    }

    fastp_save_trimmed_fail = false
    FASTP(
        FIX.out.reads,
        [],
        fastp_save_trimmed_fail,
        []
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    //FASTP.out.reads.view()

    NANOPLOT(
        FASTP.out.reads
    )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
    //NANOPLOT.out.txt.view()

    HUMAN_MAPPING(
        FASTP.out.reads
    )
    //HUMAN_MAPPING.out.nohumanreads.view()

    METAMAPS_MAP(
        HUMAN_MAPPING.out.nohumanreads
    )

    meta_with_othermeta = METAMAPS_MAP.out.metaclass.join(METAMAPS_MAP.out.otherclassmeta)
    meta_with_othermeta_with_metalength = meta_with_othermeta.join(METAMAPS_MAP.out.metalength)
    meta_with_othermeta_with_metalength_with_parameter = meta_with_othermeta_with_metalength.join(METAMAPS_MAP.out.metaparameters)
    //meta_with_othermeta_with_metalength_with_parameter.view()

    METAMAPS_CLASSIFY(
        meta_with_othermeta_with_metalength_with_parameter
    )


    //METAMAPS_MAP.out.metaclass.view()
    //NANOPLOT.out.totalreads.view()
    //METAMAPS_CLASSIFY.out.classlength.view()
    //METAMAPS_CLASSIFY.out.classcov.view()

    rmetaplot_ch=((METAMAPS_MAP.out.metaclass.join(METAMAPS_CLASSIFY.out.classlength)).join(METAMAPS_CLASSIFY.out.classcov)).join(NANOPLOT.out.totalreads)

    R_METAPLOT(
        rmetaplot_ch
    )

    KRONA_KRONADB();
    KRONA_KTIMPORTTAXONOMY(
        METAMAPS_CLASSIFY.out.classkrona,
        KRONA_KRONADB.out.db
    )

    reffasta_ch=(R_METAPLOT.out.reporttsv.join(METAMAPS_CLASSIFY.out.classem)).join(HUMAN_MAPPING.out.nohumanreads)
    //reffasta_ch.view()

    REF_FASTA(
        reffasta_ch
    )

    //REF_FASTA.out.headereads
    //REF_FASTA.out.allreads.view()
    //REF_FASTA.out.virusout.view()

    //SEQTK_SUBSEQ
    headers_ch = REF_FASTA.out.headereads.flatMap { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        entry[1].collect { virus ->
            [[id: id, single_end: singleEnd, virus: (virus.getBaseName()).replaceFirst(/.+\./,"")], "${virus}"]
        }
    }//.view()

    fasta_ch = REF_FASTA.out.seqref.flatMap { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        entry[1].collect { virus ->
            [[id: id, single_end: singleEnd, virus: ((virus.getBaseName()).replaceFirst(/\.REF+/,"")).replaceFirst(/.+\./,"")], "${virus}"]
        }
    }//.view()

    fastq_ch = REF_FASTA.out.allreads.flatMap { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        entry[1].collect { virus ->
            [[id: id, single_end: singleEnd, virus: (virus.getBaseName()).replaceFirst(/.+\./,"")], "${virus}"]
        }
    }//.view()

    REFFIX_FASTA(
        fasta_ch
    )
    //REFFIX_FASTA.out.fixedseqref.view()


    //fastq_ch.join(headers_ch).view()
    SEQTK_SUBSEQ(
        fastq_ch.join(headers_ch)
    )
    //SEQTK_SUBSEQ.out.sequences.view()

    //SEQTK_SUBSEQ.out.sequences.join(REFFIX_FASTA.out.fixedseqref).view()
    MEDAKA(
        SEQTK_SUBSEQ.out.sequences.join(REFFIX_FASTA.out.fixedseqref)
    )
    //MEDAKA.out.bamfiles.view()

    SAMTOOLS_COVERAGE(
        MEDAKA.out.bamfiles
    )
    //SAMTOOLS_COVERAGE.out.coverage.view()

    //MEDAKA.out.bamfiles.join(REFFIX_FASTA.out.fixedseqref).view()

    savempileup = false
    IVAR_CONSENSUS(
        MEDAKA.out.bamfiles.join(REFFIX_FASTA.out.fixedseqref),
        savempileup
    )
    //IVAR_CONSENSUS.out.fasta.view()

    //IVAR_CONSENSUS.out.fasta.join(REFFIX_FASTA.out.fixedseqref).view()
    HOMOPOLISH_POLISHING(
        IVAR_CONSENSUS.out.fasta.join(REFFIX_FASTA.out.fixedseqref)
    )
    //HOMOPOLISH_POLISHING.out.polishconsensus.view()



    //covcon_ch = SAMTOOLS_COVERAGE.out.coverage.join(HOMOPOLISH_POLISHING.out.polishconsensus)
    covcon_ch = (SAMTOOLS_COVERAGE.out.coverage.join(HOMOPOLISH_POLISHING.out.polishconsensus)).map { entry ->
    [[id: entry[0].id, single_end: entry[0].single_end], entry[1], entry[2]]
    }//.view()

    //covcon_ch.combine(R_METAPLOT.out.reporttsv, by: 0).view()
    
    addingdepthin_ch = (covcon_ch.combine(R_METAPLOT.out.reporttsv, by: 0)).map { entry ->
        def id = entry[0].id
        def singleEnd = entry[0].single_end
        def virus = entry[1].getBaseName().replaceFirst(/.+\./,"")
        [[id: id, single_end: singleEnd, virus: virus], entry[1], entry[2], entry[3]]
    }//.view()

    ADDING_DEPTH(
        addingdepthin_ch
    )
    
    //(ADDING_DEPTH.out.repdepth.map{it[1]}).collect().view()
    FINAL_REPORT(
        (ADDING_DEPTH.out.repdepth.map{it[1]}).collect()
    )
    //FINAL_REPORT.out.finalReport.view()

    BAM_READCOUNT(
        MEDAKA.out.bamfiles.join(REFFIX_FASTA.out.fixedseqref)
    )
    //BAM_READCOUNT.out.bamcount.view()






















    ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionsmini)
    ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamsort)
    ch_versions = ch_versions.mix(HUMAN_MAPPING.out.versionssamfastq)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatropics.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatropics.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
