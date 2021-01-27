nextflow.enable.dsl=2

// params.vcfFile = null
// params.buildCheckFile = null
// params.outDir = null

process LIFT_VCFS_TO_GCRHH37 {
  publishDir "$params.outDir"

  input:
    path vcfFile
    path buildCheckFile
  // output:
  //   file 'lift/*'

  shell:
    '''
myinput=!{vcfFile}
buildcheck=!{buildCheckFile}
myoutdir=!{params.outDir}


#Examples of the variables needed (-v)
#myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf
#buildcheck=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked
#myoutdir=/stsi/raqueld/1_lift

#How to run Example
#qsub 1_lift_vcfs_to_GRCh37.job -v myinput=/stsi/raqueld/vcf/6800_JHS_all_chr_sampleID_c2.vcf,buildcheck=/stsi/raqueld/0_check_vcf_build/6800_JHS_all_chr_sampleID_c2.BuildChecked,myoutdir=/stsi/raqueld/1_lift,custom_temp=/mnt/stsi/stsi0/raqueld/tmp -N 1_6800_JHS_all_chr_sampleID_c2


date
echo "Running on node:"
hostname
pwd

# newgrp tlabdbgap

# module purge
# module load ucsc_tools/373
# module load samtools


export filename=$(basename $buildcheck)

export plink="plink"
export plink2="plink2"



LIFT_OVER () {
    lift="/app/required_tools/lift/LiftMap.py"
    cpath="/app/required_tools/chainfiles"
    cfilename=$(grep "Use chain file" $buildcheck | tr -d ' ' | tr ':' '\t' | tr -d '"' | cut -f 2 | sed -e 's/->/ /g')
    nchains=$(echo $cfilename | tr ' ' '\n' | wc -l | awk '{print $1}')
    checknone=$(grep "Use chain file" $buildcheck | grep "none" | wc -l)

    name="chr$1"

    echo ${name}

    echo "Remove multi-allelic variants"
    bcftools view $name.sorted.vcf.gz -M 2 -m 2 | bcftools norm /dev/stdin -d both -Oz -o $name.sorted.bi.vcf.gz
    tabix -p vcf $name.sorted.bi.vcf.gz

    echo "Converting vcf to ped/map.."
    bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' $name.sorted.bi.vcf.gz > $name.sorted.bi.pos
    # plink --vcf $name.sorted.bi.vcf --make-bed --a1-allele $name.sorted.bi.pos 5 3 --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing --out $name.sorted.bi
    # plink --bfile $name.sorted.bi --recode --a1-allele $name.sorted.bi.bim 5 2 --double-id --set-missing-var-ids @:#:\$1:\$2 --out $name
    $plink --vcf $name.sorted.bi.vcf.gz --make-bed --a1-allele $name.sorted.bi.pos 5 3 --biallelic-only strict --set-missing-var-ids @:#:\$1:\$2 --vcf-half-call missing --double-id --recode ped --id-delim '_' --out $name
    # NOTE: '_' is FID_IID deliminator, keep watching if exception ID appeared.
    # TODO update to plink2 when it supported ped files.

    if [ $checknone -eq 1 ]; then
        echo "The data set is already based on the correct reference build (Grch37). Just converting and copying it."
        # plink --bfile $name.sorted.bi --make-bed --a1-allele $name.sorted.bi.bim 5 2 --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --out $name.lifted_already_GRCh37.sorted.with_dup
        $plink2 --vcf $name.sorted.bi.vcf.gz --make-bed --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --out $name.lifted_already_GRCh37.sorted.with_dup

        cut -f 2 $name.lifted_already_GRCh37.sorted.with_dup.bim | sort | uniq -d > $name.list_multi_a_markers.txt
        ndup=$(wc -l $name.list_multi_a_markers.txt | awk '{print $1}')

        if [ "$ndup" -ge 1 ]; then
            echo "Still found duplicate variant ids or multiallelic markers after filtering. Performing additional filtering."
            dupflag=$(echo -e "--exclude $name.list_multi_a_markers.txt")
        else
            dupflag=""
        fi
        # plink --bfile $name.lifted_already_GRCh37.sorted.with_dup $dupflag --a1-allele $name.sorted.bi.bim 5 2 --make-bed --out $name.lifted_already_GRCh37
        $plink2 --bfile $name.lifted_already_GRCh37.sorted.with_dup $dupflag --make-bed --out $name.lifted_already_GRCh37

        for bfile in bed bim fam; do
            mv $name.lifted_already_GRCh37.${bfile} lifted_already_GRCh37.chr${1}.${bfile}
        done

        return 1
    fi

    echo "Number of chain files: $nchains"
    echo "Lifting ped/map using chain file(s): $cfilename..."
    if [ "$nchains" -eq 1 ]; then
         liftname=$(echo $cfilename |  sed -e 's/\.chain.*//g')
         $lift -p $name.ped -m $name.map -c $cpath/$cfilename -o $name.lifted_$liftname

    fi
    if [ "$nchains" -eq 2 ]; then
         first=$(echo $cfilename | tr ' ' '\n' | head -n 1)
         second=$(echo $cfilename | tr ' ' '\n' | tail -n 1)
         liftname=$(echo $first |  sed -e 's/\.chain.*//g')
         $lift -m $name.map -c $cpath/$first -o $name.lifted_$liftname
         liftname=$(echo $second |  sed -e 's/\.chain.*//g')
         $lift -m $name.lifted_$liftname.map -c $cpath/$second -o $name.lifted_$liftname
    fi

    echo "Converting lifted output $name.lifted_$liftname to bed/bim/fam.."
    $plink --file $name.lifted_$liftname --a1-allele $name.sorted.bi.pos 5 3 --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --make-bed --out $name.lifted_$liftname.sorted.with_dup
    # TODO: update to plink2 when it supported ped files.

    cut -f 2 $name.lifted_$liftname.sorted.with_dup.bim | sort | uniq -d > $name.list_multi_a_markers.txt
        ndup=$(wc -l $name.list_multi_a_markers.txt | awk '{print $1}')
        if [ "$ndup" -ge 1 ]; then
            echo "Still found duplicate variant ids or multiallelic markers after filtering. Performing additional filtering."
            dupflag=$(echo -e "--exclude $name.list_multi_a_markers.txt")
            # plink --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --a1-allele $name.sorted.bi.bim 5 2 --allow-extra-chr --out $name.lifted_$liftname
            $plink2 --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --allow-extra-chr --out $name.lifted_$liftname
        else
            echo "No duplicate. Skip additional plink, move the filename directly."
            dupflag=""
            # plink --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --a1-allele $name.sorted.bi.bim 5 2 --allow-extra-chr --out $name.lifted_$liftname
            $plink2 --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --allow-extra-chr --out $name.lifted_$liftname
            for bfile in bed bim fam; do
                mv $name.lifted_$liftname.sorted.with_dup.${bfile} lifted_$liftname.chr$1.${bfile}
            done
        fi


}
export -f LIFT_OVER


if [ ! -z $custom_temp ]; then
    mkdir -p $custom_temp/$SLURM_JOBID
    export TEMP=$custom_temp/$SLURM_JOBID
else
    mkdir -p $TMP/$SLURM_JOBID
    export TEMP=$TMP/$SLURM_JOBID
fi

if [ ! -d $myoutdir ]; then
    mkdir -p $myoutdir
fi

cd $myoutdir

vcfstarttime=$(date +%s)
echo $myinput
# Test if input came as split chromosome
if [[ ${myinput} == *.vcf ]]; then
    bgzip -c ${myinput} > myvcf.vcf.gz
elif [[ ${myinput} == *.vcf.gz ]]; then
    cp ${myinput} myvcf.vcf.gz
    echo "Input file is merged; split by chromsome."
else
    echo "Invalid file format, please check the input."
    exit
fi

ls -l

# CHECK_SORT ALL myvcf.vcf.gz
bcftools sort myvcf.vcf.gz -T ${TEMP} -Oz -o ./myvcf.sorted.vcf.gz
tabix -p vcf ./myvcf.sorted.vcf.gz
for i in {1..22}; do
  bcftools view ./myvcf.sorted.vcf.gz -r ${i} -Oz -o chr${i}.sorted.vcf.gz;
done
# parallel SPLIT_CHR ./${inprefix}.chrALL.sorted.vcf.gz {1} ::: {1..22}


vcfendtime=$(date +%s)


liftstarttime=$(date +%s)
for i in {1..22}; do
  LIFT_OVER $i;
done
# parallel LIFT_OVER {1} ::: {1..22}
liftendtime=$(date +%s)


rm $name.map
rm *with_dup*
rm $name.lifted_$liftname.sorted.*
rm $name.list_multi_a_markers.txt
rm $name.sorted.bi.*


vcfruntime=$((vcfendtime-vcfstarttime))
liftruntime=$((liftendtime-liftstarttime))
echo "vcf convert time: $vcfruntime"
echo "liftover run time: $liftruntime"

#output example
#$inprefix.lifted_$liftname.chr${chrom}.${bfile}
    '''
    

  stub:
    """
      touch dummy.BuildChecked
    """
}

workflow test{

  params.buildCheckFile = "output.BuildChecked"
  params.vcfFile = "/root/A4420_2020_1.vcf"
  //params.outDir = "/root/lift"
  vcfFileChannel = Channel.fromPath(params.vcfFile)
  buildCheckFileChannel = Channel.fromPath(params.buildCheckFile)
  // publishDir = params.outDir


  LIFT_VCFS_TO_GCRHH37(vcfFileChannel, buildCheckFileChannel)

}
