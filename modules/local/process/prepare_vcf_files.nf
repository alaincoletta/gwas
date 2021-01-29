nextflow.enable.dsl=2

/*
 * prepare VCF files split by chromosome and GZIP
 */


process PREPARE_VCF_FILES{

    publishDir "${params.outdir}"
    container "geneplaza/imputation2"
    
    input:
    path vcf_file
    
    output:
    path  "myvcf.vcf.gz" 

    shell:
    '''
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
    params.outdir = "${launchDir}/testData"
    params.vcf_file = "${launchDir}/testData/myvcf.vcf"
    ch_vcf_file = channel.fromPath(params.vcf_file)

   PREPARE_VCF_FILES(ch_vcf_file)
//    PREPARE_VCF_FILES.out.dummyOuput 
    
//  CHECK_SORT_VCF(
//          PREPARE_VCF_FILES.out.dummyOuput  
//    )


}