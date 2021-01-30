nextflow.enable.dsl=2
// // Import generic module functions
// include { initOptions; saveFiles } from './functions'

// params.options = [:]
// def options    = initOptions(params.options)

/*
 * check sort VCF files using SamTools  bcftools and tabix commands
 */

 process CHECK_SORT_VCF{
        publishDir "${params.outdir}"
        //container "gwas/imputation2"

        input:
        path mergedVcfFile 
    
        shell:
        '''
    myinput=!{mergedVcfFile}

    bcftools sort myvcf.vcf.gz -t ${temp} -oz -o ./myvcf.sorted.vcf.gz
    tabix -p vcf ./myvcf.sorted.vcf.gz
    for i in {1..22}; do
      bcftools view ./myvcf.sorted.vcf.gz -r ${i} -oz -o chr${i}.sorted.vcf.gz;
    done
        '''

        stub:
        """
            touch chr1.sorted.vcf.gz
        """
 }

 workflow test{
    include { PREPARE_VCF_FILES } from "$baseDir/module/process/prepare_vcf_files.nf"
     
    params.outdir = "lift"
    params.mergedVcfFile = "myvcf.vcf.gz"
    ch_mergedVcfFile = channel.fromPath(params.mergedVcfFile)

    
 CHECK_SORT_VCF(
         PREPARE_VCF_FILES.out.dummyOuput 
   )


 }