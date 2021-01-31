shell:
'''
cpath="/app/required_tools/chainfiles"
    cfilename=$(grep "use chain file" $buildcheck | tr -d ' ' | tr ':' '\t' | tr -d '"' | cut -f 2 | sed -e 's/->/ /g')
    nchains=$(echo $cfilename | tr ' ' '\n' | wc -l | awk '{print $1}')
    checknone=$(grep "use chain file" $buildcheck | grep "none" | wc -l)
'''