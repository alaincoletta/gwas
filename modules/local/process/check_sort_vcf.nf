nextflow.enable.dsl=2
// // Import generic module functions
// include { initOptions; saveFiles } from './functions'

// params.options = [:]
// def options    = initOptions(params.options)

/*
 * check sort VCF files using SamTools  bcftools and tabix commands
 */

 process CHECK_SORT_VCF{
        publishDir "${params.outdir}", mode:'copy'
        container "gwas/imputation2"

        input:
        path mergedVcfFile 
    
        shell:
        '''
    myinput=!{mergedVcfFile}
    myoutput=!{params.myOutput}

    bcftools sort ${myinput} -t ${temp} -oz -o ${myOutput}
    tabix -p vcf ${myOutput} 
    for i in {1..22}; do
      bcftools view ${myOutput} -r ${i} -oz -o chr${i}.sorted.vcf.gz;
    done
        '''

        stub:
        """
            touch chr1.sorted.vcf.gz
        """
 }

 workflow test{
    include { PREPARE_VCF_FILES } from "$launchDir/modules/local/process/prepare_vcf_files.nf"
     
    params.outdir = "$launchDir/testData/lift"
    // params.mergedVcfFile = "myvcf.vcf.gz"   --> use output from PREPARE_VCF_FILES
    //ch_mergedVcfFile = channel.fromPath(params.mergedVcfFile)

    
 CHECK_SORT_VCF(
         PREPARE_VCF_FILES.out
   )


 }