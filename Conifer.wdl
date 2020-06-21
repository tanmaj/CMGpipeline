version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION 
workflow Conifer {
  input {
    File input_bam
    File input_bam_index

    String? enrichment
    File? enrichment_bed
  }  

  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )

  if (defined(enrichment_bed)) {
    call MakeRPKM {
      input:
        input_bam=input_bam,
        input_bam_index=input_bam_index,
        sample_basename=sample_basename,
        enrichment=enrichment,
        enrichment_bed=enrichment_bed
    }
  }
}


##################
# TASK DEFINITIONS
##################

task MakeRPKM {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    String? enrichment
    File? enrichment_bed

    # Runtime parameters
    String docker = "ielis/conifer:latest"
  }
  
  command {
  set -e
  python /root/conifer/conifer.py rpkm --probes ~{enrichment_bed} --input ~{input_bam} --output ~{sample_basename}.~{enrichment}.rpkm 
  }
  runtime {
    docker: docker
  }
  output {
    File output_rpkm = "sample_basename.~{enrichment}.rpkm"
  }
}

