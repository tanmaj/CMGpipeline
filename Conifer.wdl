version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION 
workflow Conifer {
  input {
    File input_bam
    File input_bam_index

    Array[File] input_reference_rpkms 
    Int CONIFER_svd
    Float CONIFER_threshold

    String enrichment
    File enrichment_bed
  }  

  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )

  call MakeRPKM {
      input:
        input_bam=input_bam,
        input_bam_index=input_bam_index,
        sample_basename=sample_basename,
        enrichment=enrichment,
        enrichment_bed=enrichment_bed
  }

  call CONIFER_Analyze {
      input:
        input_rpkm=MakeRPKM.output_rpkm,
        input_reference_rpkms=input_reference_rpkms,
        CONIFER_svd=CONIFER_svd,
        sample_basename=sample_basename,
        enrichment=enrichment,
        enrichment_bed=enrichment_bed

  }

  call CONIFER_Call {
      input:
        input_hdf5=CONIFER_Analyze.output_hdf5,
        CONIFER_threshold=CONIFER_threshold,
        sample_basename=sample_basename
  }

    call CONIFER_Plotcalls {
      input:
        input_hdf5=CONIFER_Analyze.output_hdf5,
        input_conifer_calls=CONIFER_Call.output_conifer_calls,
        sample_basename=sample_basename
  }

  output {
    File output_conifer_calls = CONIFER_Call.output_conifer_calls
    Array[File] output_plotcalls = CONIFER_Plotcalls.output_plotcalls
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

    String enrichment
    File enrichment_bed

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command {
  set -e
  python /home/bio/conifer_v0.2.2/conifer.py rpkm --probes ~{enrichment_bed} --input ~{input_bam} --output ~{enrichment}_~{sample_basename}.txt 
  }
  runtime {
    docker: docker
  }
  output {
    File output_rpkm = "~{enrichment}_~{sample_basename}.txt"
  }
}

task CONIFER_Analyze {
  input {
    # Command parameters
    File input_rpkm
    Array[File] input_reference_rpkms 
    String CONIFER_svd
    String sample_basename

    String enrichment
    File enrichment_bed

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command <<<
  set -e
  RPKM_DIR=$(dirname ~{input_reference_rpkms[0]})
  cp ~{input_rpkm} $RPKM_DIR

  python /home/bio/conifer_v0.2.2/conifer.py analyze --probes ~{enrichment_bed} --rpkm_dir $RPKM_DIR --output ~{sample_basename}.analysis.hdf5 --svd ~{CONIFER_svd} --write_svals ~{sample_basename}.singular_values.txt --plot_scree ~{sample_basename}.screeplot.png --write_sd ~{sample_basename}.sd_values.txt
  >>>

  runtime {
    docker: docker
  }
  output {
    File output_hdf5 = "~{sample_basename}.analysis.hdf5"
  }
}

task CONIFER_Call {
  input {
    # Command parameters
    File input_hdf5
    Float CONIFER_threshold
    String sample_basename

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command {
  set -e
  python /home/bio/conifer_v0.2.2/conifer.py call --threshold ~{CONIFER_threshold} --input ~{input_hdf5} --output ~{sample_basename}.CONIFER_CALLS_POPULATION.txt

  head -n 1 ~{sample_basename}.CONIFER_CALLS_POPULATION.txt > ~{sample_basename}.CONIFER_CALLS.txt
  cat ~{sample_basename}.CONIFER_CALLS_POPULATION.txt | grep ~{sample_basename} >> ~{sample_basename}.CONIFER_CALLS.txt
  }

  runtime {
    docker: docker
  }
  output {
    File output_conifer_calls = "~{sample_basename}.CONIFER_CALLS.txt"
  }
}

task CONIFER_Plotcalls {
  input {
    # Command parameters
    File input_hdf5
    File input_conifer_calls
    String sample_basename

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command {
  set -e
  export LC_ALL="en_US.UTF-8"
  export LC_CTYPE="en_US.UTF-8"
  python /home/bio/conifer_v0.2.2/conifer.py plotcalls --input ~{input_hdf5} --calls ~{input_conifer_calls} --output ./
  }

  runtime {
    docker: docker
  }
  output {
    Array[File] output_plotcalls = glob("*.png")
  }
}
