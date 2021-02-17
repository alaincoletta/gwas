nextflow.enable.dsl=2

/*
 * prepare VCF files ONLY happens if input is split by chromosome and/or input is not compressed
 */


process PREPARE_VCF_FILES{

    publishDir "${params.outputDir}"
    container "gwas/imputation2"

    input:
    path(vcfFile)

    output:
        path("*output.gz")

    shell:
    '''
echo !{vcfFile}
myinput=!{vcfFile}
myoutput=!{vcfFile}.output.gz
# test if input came as split chromosome
if [[ ${myinput} == *.vcf ]]; then
      bgzip -c ${myinput} > ${myoutput}
elif [[ ${myinput} == *.vcf.gz ]]; then
    cp ${myinput} ${myoutput}
    echo "input file is merged; split by chromsome."
else
    echo "invalid file format, please check the input."
    exit
fi
    '''


    stub:
    """
        touch myvcf.vcf.output.gz
    """
}

workflow test{

    params.outputDir = "$launchDir/testData/"
    params.vcfFile = "$launchDir/testData/A4420_2020_1.vcf"


   PREPARE_VCF_FILES(Channel.of([params.vcfFile]))

}
