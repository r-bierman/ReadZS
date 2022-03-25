process FILTER_BAM_10X {
  tag "${chr}, ${bamFileID}"
  label 'process_medium'

  //conda 'python=3.7 bioconda::pysam pandas'

  input:
  tuple val(inputChannel), val(bamFileID), path(bam), path(bai)
  val isSICILIAN
  val isCellranger
  val libType
  each chr

  output:
  path "*.filter", emit: filter

  script:
  """
  /share/software/user/open/python/3.6.1/bin/python3 /oak/stanford/groups/horence/rob/readzs_fork/ReadZS/bin/filter.py \\
    --input_bam <(samtools view -b ${bam} ${chr} ) \\
    --isSICILIAN ${isSICILIAN} \\
    --isCellranger ${isCellranger} \\
    --libType ${libType} \\
    --bamName ${bamFileID} \\
    --inputChannel ${inputChannel} \\
    --chr ${chr}
  """
}
