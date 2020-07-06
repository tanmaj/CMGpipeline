## # QualiMap
##
## This WDL tool wraps the [QualiMap](http://qualimap.bioinfo.cipf.es/) tool.
## QualiMap computes metrics to facilitate evaluation of sequencing data. 
version 1.0

task bamqc {
    input {
        File bam
        File enrichment_bed
        String sample_basename
        Int ncpu = 1
        Int max_retries = 1
        Int memory_gb = 8
        Int? disk_size_gb
    }

    String out_directory = basename(bam, ".bam") + '_qualimap_bamqc_results'
    String out_tar_gz_file = sample_basename + ".Qualimap.tar.gz"

    Int java_heap_size = ceil(memory_gb * 0.9)
    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command {
        set -euo pipefail
        
        qualimap bamqc -bam ~{bam} \
            -gff ~{enrichment_bed} \
            -outdir ~{out_directory} \
            -nt ~{ncpu} \
            -nw 400 \
            --java-mem-size=~{java_heap_size}g
        
        # Check if qualimap succeeded
        if [ ! -d "~{out_directory}/raw_data_qualimapReport/" ]; then
            exit 1
        fi

        tar -czf ~{out_tar_gz_file} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/qualimap:1.0.0'
        maxRetries: max_retries
    }

    output {
        File results = out_tar_gz_file
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs QualiMap's bamqc tool on the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

## An additional option for calculating coverage using GATK
task DepthOfCoverage {
    input {
        File input_bam
        File input_bam_index
        String sample_basename

        File reference_fa
        File reference_fai
        File reference_dict
        
        File? enrichment_bed
        File refSeqFile
                
        # Runtime params
        String gatk_path
        Int threads
        String docker
    }

    command {
    set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      DepthOfCoverage \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      -O ~{sample_basename}.coverage \
      ~{"-L " + enrichment_bed} \
       --omit-depth-output-at-each-base \
       -ip 2 \
       --summary-coverage-threshold 5 \
       --summary-coverage-threshold 10 \
       --summary-coverage-threshold 15 \
       --summary-coverage-threshold 20 \
       --summary-coverage-threshold 30 \
       --summary-coverage-threshold 40 \
       --summary-coverage-threshold 50 \
       --summary-coverage-threshold 60 \
       --summary-coverage-threshold 70 \
       --summary-coverage-threshold 80 \
       --summary-coverage-threshold 90 \
       --summary-coverage-threshold 100 \
       --gene-list ~{refSeqFile}

       tar -czf ~{sample_basename}.DepthOfCoverage.tar.gz ~{sample_basename}.coverage*
    }
    
    output {
    File DepthOfCoverage_output = "~{sample_basename}.DepthOfCoverage.tar.gz"
    }

    runtime {
        docker: "~{docker}"
        maxRetries: 3
    }
}