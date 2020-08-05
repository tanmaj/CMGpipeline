## # QualiMap
##
## This WDL tool wraps the [QualiMap](http://qualimap.bioinfo.cipf.es/) tool.
## QualiMap computes metrics to facilitate evaluation of sequencing data. 
version 1.0


## An additional option for calculating coverage using GATK
task CreateMitoFasta {
    input {
        File input_vcf
        String sample_basename

        File reference_fa
        File reference_fai
        File reference_dict
        
        File? enrichment_bed
        File refSeqFile
    }

    command <<<
    set -e
    java -Xmx5g -jar /usr/GenomeAnalysisTK.jar \
      -T FastaAlternateReferenceMaker \
      -R ~{reference_fa} \
      --variant ~{input_vcf} \
      -o mtDNA.fasta \
      -L chrM
     >>>
    
    output {
    File mtDNA_fasta = "mtDNA.fasta"
    }

    runtime {
        docker: "~{docker}"
        maxRetries: 3
        requested_memory_mb_per_core: 6000
        cpu: 1
        runtime_minutes: 60
    }
}

task MitoMap {
    input {
        File mtDNA_fasta
    }

    command {
    set -e
	cp /usr/src/app/mitomap.py ./
	cp ~{mtDNA_fasta} ./
	python mitomap.py
    }

    output {
    File mtDNA_fasta = "mtDNA.fasta"
    }

    runtime {
        docker: "alesmaver/mitomap"
        maxRetries: 3
        requested_memory_mb_per_core: 3000
        cpu: 1
        runtime_minutes: 60
    }
}