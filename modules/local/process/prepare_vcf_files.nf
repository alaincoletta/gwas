nextflow.enable.dsl=2

/*
 * prepare VCF files split by chromosome and GZIP
 */


process PREPARE_VCF_FILES{

    publishDir "${launchDir}/testData"
    //container "gwas/imputation2"
    
    input:
    path vcfFile
    
    output:
    path  "myvcf.vcf.gz" 

    shell:
    '''
myinput = !{vcfFile}
# test if input came as split chromosome
if [[ ${myinput} == *.vcf ]]; then
    bgzip -c ${myinput} > myvcf.vcf.gz
elif [[ ${myinput} == *.vcf.gz ]]; then
    cp ${myinput} myvcf.vcf.gz
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
    params.vcfFile = "${launchDir}/testData/myvcf.vcf"
    ch_vcfFile = channel.fromPath(params.vcfFile)

   PREPARE_VCF_FILES(ch_vcfFile)

}