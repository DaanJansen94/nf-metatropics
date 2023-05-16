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

//include { MINIMAP2_ALIGN              } from '../modules/nf-core/minimap2/align/main'
//include { SAMTOOLS_SORT               } from '../modules/nf-core/samtools/sort/main'
//include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
//include { SAMTOOLS_FASTQ              } from '../modules/nf-core/samtools/fastq/main'
//include { FIX_NAMES                   } from '../modules/local/fix_names'
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
