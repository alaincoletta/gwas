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
        container "geneplaza/imputation2"

        input:
        path merged_vcf_file 
    
        shell:
        '''
    myinput=!{merged_vcf_file}

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
     params.outdir = 'outdir'
     params.merged_vcf_file = "myvcf.vcf.gz"
     ch_merged_vcf_file = channel.fromPath(params.merged_vcf_file)

    
 CHECK_SORT_VCF(
         PREPARE_VCF_FILES.out.dummyOuput  
   )


 }