hell:
'''
bcftools view $name.sorted.vcf.gz -m 2 -m 2 | bcftools norm /dev/stdin -d both -oz -o $name.sorted.bi.vcf.gz
    tabix -p vcf $name.sorted.bi.vcf.gz
'''