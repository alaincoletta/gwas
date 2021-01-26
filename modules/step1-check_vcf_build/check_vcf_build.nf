nextflow.enable.dsl = 2

//params.vcfFile

process CHECK_VCF_BUILD {
    input:
    path(vcfFile)

    output:
    path("*.BuildChecked")

    shell:
        '''
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

    params.myoutput = "output.BuildChecked"
    rawFileChannel = Channel.fromPath(params.vcfFile)
    check_vcf_build(rawFileChannel)

}
