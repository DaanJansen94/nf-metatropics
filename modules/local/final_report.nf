process FINAL_REPORT {
    //tag "${meta.id}"

    input:
    path(report)

    output:
    path("all.final_report.tsv"), emit: finalReport

    script:
    """
    echo -e "Sample\tAccession\tTaxID\tVirusName\tMappedReads\tFractionMappedReads\tAbundance\tCoverage\tDepthAverage\tConsensusCov\tN_content\tMedianReadIdentities\tMeanReadLength\tMeanBaseQuality" > all.final_report.tsv
    cat *.sdepth.tsv |  grep -v VirusName >> all.final_report.tsv
    """
}
