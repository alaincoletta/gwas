shell:
'''
first=$(echo $cfilename | tr ' ' '\n' | head -n 1)
         second=$(echo $cfilename | tr ' ' '\n' | tail -n 1)
         liftname=$(echo $first |  sed -e 's/\.chain.*//g')
         $lift -m $name.map -c $cpath/$first -o $name.lifted_$liftname
         liftname=$(echo $second |  sed -e 's/\.chain.*//g')
         $lift -m $name.lifted_$liftname.map -c $cpath/$second -o $name.lifted_$liftname
'''