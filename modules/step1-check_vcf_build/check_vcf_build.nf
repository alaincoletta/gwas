nextflow.enable.dsl = 2

//params.vcfFile

process CHECK_VCF_BUILD {
    publishDir "${params.outputDir}"
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
if [[ ${myinput} == *.vcf ]]; then
    echo "Input autosome list detected"
    my_chr1='my_chr1.txt'
    cat ${myinput} | awk '$1==1 {print$0}' > $my_chr1
    mydirname=$PWD
    filename=$(basename ${my_chr1})
    echo $PWD/$filename
elif [[ ${myinput} == *.vcf.gz ]]; then
    echo "Merged vcf.gz file detected, chr1 will be extract in the Rscript."
    mydirname=$(dirname $myinput)
    filename=$(basename $myinput)
else
    echo "Invalid file format, please check the input."
    exit
fi

Rscript /app/required_tools/check_vcf_build/check_vcf_build.R $mydirname/$filename > !{myOutput}

        '''

    stub:
    """
    touch output.BuildChecked
    """
}

workflow test {
    params.outputDir = "$launchDir/testData/"
    params.vcfFile = "$launchDir/testData/myvcf.vcf"
    params.myOutput = "$launchDir/testData/output.BuildChecked"
    
    check_vcf_params_ch = Channel.of([
            params.vcfFile,
            params.myOutput,
    ])

    CHECK_VCF_BUILD(check_vcf_params_ch)

}
