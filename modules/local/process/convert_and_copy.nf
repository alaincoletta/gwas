shell:
'''
#--------------------------------process CONVERT_AND_COPY (lifted files uses plink)------------------------------#
  
    echo "converting lifted output $name.lifted_$liftname to bed/bim/fam.."
    $plink --file $name.lifted_$liftname --a1-allele $name.sorted.bi.pos 5 3 --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --make-bed --out $name.lifted_$liftname.sorted.with_dup
    # todo: update to plink2 when it supported ped files.

    cut -f 2 $name.lifted_$liftname.sorted.with_dup.bim | sort | uniq -d > $name.list_multi_a_markers.txt
        ndup=$(wc -l $name.list_multi_a_markers.txt | awk '{print $1}')
        if [ "$ndup" -ge 1 ]; then
            echo "still found duplicate variant ids or multiallelic markers after filtering. performing additional filtering."
            dupflag=$(echo -e "--exclude $name.list_multi_a_markers.txt")
            $plink2 --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --allow-extra-chr --out $name.lifted_$liftname
        else
            echo "no duplicate. skip additional plink, move the filename directly."
            dupflag=""
            $plink2 --bfile $name.lifted_$liftname.sorted.with_dup $dupflag --make-bed --allow-extra-chr --out $name.lifted_$liftname
            
            for bfile in bed bim fam; do
                mv $name.lifted_$liftname.sorted.with_dup.${bfile} lifted_$liftname.chr$1.${bfile}
            done
        fi
    #-------------------------------------------------------------------------------------------------------------#
#--------------------------------process CONVERT_AND_COPY (non-lifted files uses plink2)-----------------------#
    if [ $checknone -eq 1 ]; then
        echo "the data set is already based on the correct reference build (grch37). just converting and copying it."
        $plink2 --vcf $name.sorted.bi.vcf.gz --make-bed --double-id --set-missing-var-ids @:#:\$1:\$2 --allow-extra-chr --out $name.lifted_already_grch37.sorted.with_dup

        cut -f 2 $name.lifted_already_grch37.sorted.with_dup.bim | sort | uniq -d > $name.list_multi_a_markers.txt
        ndup=$(wc -l $name.list_multi_a_markers.txt | awk '{print $1}')

        if [ "$ndup" -ge 1 ]; then
            echo "still found duplicate variant ids or multiallelic markers after filtering. performing additional filtering."
            dupflag=$(echo -e "--exclude $name.list_multi_a_markers.txt")
        else
            dupflag=""
        fi
        $plink2 --bfile $name.lifted_already_grch37.sorted.with_dup $dupflag --make-bed --out $name.lifted_already_grch37

        for bfile in bed bim fam; do
            mv $name.lifted_already_grch37.${bfile} lifted_already_grch37.chr${1}.${bfile}
        done

        return 1
    fi

    echo "number of chain files: $nchains"
    echo "lifting ped/map using chain file(s): $cfilename..."
    if [ "$nchains" -eq 1 ]; then
         liftname=$(echo $cfilename |  sed -e 's/\.chain.*//g')
         $lift -p $name.ped -m $name.map -c $cpath/$cfilename -o $name.lifted_$liftname

    fi
    #----------------------------------------------------------------------------------------------------------------#
'''