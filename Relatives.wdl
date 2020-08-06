#Relatives WDL

    # This line would merge multiple VCF files 
    # merge /mnt/EXOMES/PX5400/annotation/PX5400.annotated.vcf.gz /mnt/EXOMES/PX5401/annotation/PX5401.annotated.vcf.gz -Oz -o  ~{input_vcf} ~{sep=' ' relative_vcfs}

    # input_vcf = select_first([merged_vcf, input_vcf])

    # This line will create SnpSift genotype extractor for single and multi sample VCFs
    # bcftools query -l ~{input_vcf} | awk '{print "GEN["$1"].GT GEN["$1"].AD GEN["$1"].DP GEN["$1"].GQ"}' | tr '\n' ' '