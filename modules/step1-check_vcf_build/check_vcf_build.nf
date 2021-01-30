nextflow.enable.dsl = 2

//params.vcfFile

process CHECK_VCF_BUILD {
    input:
    tuple path(vcfFile),
    path(myOutput)

    output:
    path(myOutput)

    shell:
        '''
myinput=!{vcfFile}
if [[ ${myinput} == *.vcf ]]; then
    echo "Input autosome list detected, run parallel pipeline."
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

Rscript /app/required_tools/check_vcf_build/check_vcf_build.R $mydirname/$filename > !{params.myoutput}

        '''

    stub:
    """
    touch output.BuildChecked
    """
}

workflow test {

    params.vcfFile = "$baseDir/testData/myvcf.vcf"
    params.myOutput = "$baseDir/testData/output.BuildChecked"
    
    check_vcf_params_ch = Channel.of([
            file(params.vcfFile),
            file(params.myOutput),
    ])

    CHECK_VCF_BUILD(check_vcf_params_ch)

}
