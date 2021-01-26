nextflow.enable.dsl=2

params.vcfFile = null
params.buildCheckFile = null
params.outDir = null

process lift_vcfs_to_GCRHh37 {
  publishDir "$params.outDir"

  input:
    file vcfFile
    file buildCheckFile
  // output:
  //   file 'lift/*'

  script:
    """
      /app/1_lift_vcfs_to_GRCh37.slurm.sh $vcfFile $buildCheckFile .
    """

  stub:
    """
      touch dummy.BuildChecked
    """
}

workflow test{

  vcfFileChannel = Channel.fromPath(params.vcfFile)
  buildCheckFileChannel = Channel.fromPath(params.buildCheckFile)
  publishDir = params.outDir


  lift_vcfs_to_GCRHh37(vcfFileChannel, buildCheckFileChannel)

}
