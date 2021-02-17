nextflow.enable.dsl=2

/*
 * check sort VCF files using SamTools  bcftools and tabix commands
 */

 process CHECK_SORT_VCF{
        publishDir "${params.outputDir}", mode:'copy'
        container "gwas/imputation2"

        input:
        path mergedVcfFile

        output:
            path("*sorted.vcf.gz")

        shell:
        '''
    myinput=!{mergedVcfFile}
    myOutput=!{mergedVcfFile}.sorted
    mv ${myinput} input.vcf.gz
    bcftools sort input.vcf.gz -T temp -Oz -o sorted.vcf.gz
    tabix -p vcf sorted.vcf.gz
    for i in {1..22}; do
      bcftools view sorted.vcf.gz -r ${i} -Oz -o chr${i}.sorted.vcf.gz;
    done
        '''

        stub:
        """
            touch chr1.sorted.vcf.gz
        """
 }

 workflow test{
    include { PREPARE_VCF_FILES } from "$launchDir/modules/local/process/prepare_vcf_files.nf"

    //params.outdir = "$launchDir/testData/lift"
    // params.mergedVcfFile = "myvcf.vcf.gz"   --> use output from PREPARE_VCF_FILES
    //ch_mergedVcfFile = channel.fromPath(params.mergedVcfFile)
    params.outputDir = "$launchDir/testData/"
    params.vcfFile = "$launchDir/testData/A4420_2020_1.vcf"


    PREPARE_VCF_FILES(Channel.of([params.vcfFile]))

 CHECK_SORT_VCF(
         PREPARE_VCF_FILES.out
   )
 }
