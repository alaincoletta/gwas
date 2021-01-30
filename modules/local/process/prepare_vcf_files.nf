nextflow.enable.dsl=2

/*
 * prepare VCF files ONLY happens if input is split by chromosome and/or input is not compressed 
 */


process PREPARE_VCF_FILES{

    publishDir "${params.outputDir}", mode: 'copy'   
    container "gwas/imputation2"
    
    input:
    tuple path(vcfFile),
    path(myOutput)
    
    output:
    path(myOutput)

    shell:
    '''
echo "!{params.vcfFile}"
echo "!{params.myOutput}"
myinput=!{params.vcfFile}
myoutput=!{params.myOutput}
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
        touch myvcf.vcf.gz
    """
}

workflow test{
    
    params.outputDir = "$launchDir/testData/"
    params.vcfFile = "$launchDir/testData/myvcf.vcf"
    params.myOutput = "$launchDir/testData/myvcf.vcf.gz"

    prepare_vcf_params_ch = Channel.of([
            params.vcfFile,
            params.myOutput,
    ])

   PREPARE_VCF_FILES(prepare_vcf_params_ch)

}