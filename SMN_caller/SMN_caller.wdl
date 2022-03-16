version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION
workflow SMN_caller {
  input {
    File input_bam
    File input_bam_index
  }

  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )

  call SMN_caller_task {
      input:
        input_bam=input_bam,
        input_bam_index=input_bam_index,
        sample_basename=sample_basename
  }

  output {
    File output_tsv = SMN_caller_task.output_tsv
    File output_json = SMN_caller_task.output_json
  }
}



task SMN_caller_task {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    # Runtime parameters
    String docker = "alesmaver/smn_caller"
  }

  command {
    echo ~{input_bam} > manifest.txt
    python3 /SMNCopyNumberCaller/smn_caller.py --manifest manifest.txt --genome 19 --prefix ~{sample_basename}.smn_caller --outDir . --threads 6
    python3 /SMNCopyNumberCaller/smn_charts.py --sample-data ~{sample_basename}.smn_caller.json --out-dir .
    mv   smn_~{sample_basename}.marked.pdf   ~{sample_basename}.smn_caller.pdf
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 2000
    cpu: 6
    runtime_minutes: 60
  }
  output {
    File output_tsv = "~{sample_basename}.smn_caller.tsv"
    File output_json = "~{sample_basename}.smn_caller.json"
    File output_json = "~{sample_basename}.smn_caller.pdf"
  }
}
