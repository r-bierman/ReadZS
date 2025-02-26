process PICARD {
    tag "${bamName}"
    label 'process_medium'

    conda 'bioconda::picard=2.26.2'

    input:
    tuple val(inputChannel), val(bamFileID), path(bam)

    output:
    tuple val(inputChannel), val(bamFileID), path("*dedup*"), emit: bam_tuple

    script:
    outputFile = "${bamFileID}.dedup"
    metrics = "${bamFileID}.metrics"
    """
    picard MarkDuplicates -I ${bam} -O ${outputFile} -M ${metrics} --QUIET true
    """
}
