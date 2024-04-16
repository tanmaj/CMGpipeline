version 1.0

workflow ExomeDepth {
  input {
    File target_bed
    Array[File] baseline_samples
    Array[File] baseline_samples_bais
    Array[File] test_samples
    Array[File] test_samples_bais
    String chr
    File ref
    Float transition_probability
  }

  call EXOME_DEPTH {
    input:
        target_bed = target_bed,
        test_samples = test_samples,
        test_samples_bais = test_samples_bais,
        baseline_samples = baseline_samples,
        baseline_samples_bais = baseline_samples_bais,
        ref = ref,
        transition_probability = transition_probability
  }

  output {
    File exome_depth_tsv = EXOME_DEPTH.exome_depth_tsv
  }
}


task EXOME_DEPTH {
  input {
    File target_bed
    Array[File] test_samples
    Array[File] test_samples_bais
    Array[File] baseline_samples
    Array[File] baseline_samples_bais
    File ref
    Float transition_probability
  }

  command <<<
    sed -i 's/^chr//' ~{target_bed}

    wget https://raw.githubusercontent.com/egustavsson/runExomeDepth/master/ExomeDepth.R
    Rscript ExomeDepth.R \
        --targets ~{target_bed} \
        --test-samples ~{write_lines(test_samples)} \
        --baseline-samples ~{write_lines(baseline_samples)} \
        --output-directory ./
  >>>

    runtime {
        docker: "alesmaver/exome_depth"
        maxRetries: 3
        requested_memory_mb_per_core: 6000
        cpu: 1
        runtime_minutes: 60
    }

    output {
        File exome_depth_tsv = "exome_depth.tsv"
  }
}