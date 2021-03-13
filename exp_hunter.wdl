version 1.0

workflow ExpansionHunter {
  input {
    String sample_id
    File bam_file
    File bai_file
    File reference_fasta
    File variant_catalog_file
    String expansion_hunter_docker

    File repeats_file
  }

  parameter_meta {
    sample_id: "sample name"
    bam_file: ".bam file to search for repeat expansions"
    reference_fasta: ".fasta file with reference used to align bam file"
    variant_catalog_file: "JSON array whose entries specify individual loci that the program will analyze"
    expansion_hunter_docker: "expansion hunter docker including annotation software"

    repeats_file: "Repats file used for annotation with 'stranger'"
  }

  meta {
      author: "Gaber Bergant and Aleš Maver"
      email: "cmg.kimg@kclj.si"
  }
  
  call RunExpansionHunter {
      input:
        sample_id = sample_id,
        bam_file = bam_file,
        bai_file = bai_file,
        reference_fasta = reference_fasta,
        variant_catalog_file = variant_catalog_file,
        expansion_hunter_docker = expansion_hunter_docker
    }

  call AnnotateExpansionHunter {
      input:
        sample_id = sample_id,
        repeats_file = repeats_file,
        expansion_hunter_docker = expansion_hunter_docker,
        expansion_hunter_vcf = RunExpansionHunter.expansion_hunter_vcf
    }

  output {
    File? expansion_hunter_vcf_annotated = AnnotateExpansionHunter.expansion_hunter_vcf_annotated
  }

}

task RunExpansionHunter {
  input {
    String sample_id
    File bam_file
    File bai_file
    File reference_fasta
    File variant_catalog_file
    String expansion_hunter_docker
  }

  output {
    File expansion_hunter_vcf = "~{sample_id}.vcf"
  }

  command <<<

    echo "[ RUNNING ] expansion hunter on sample ~{sample_id}"
    ExpansionHunter \
      --reads ~{bam_file} \
      --reference ~{reference_fasta} \
      --variant-catalog ~{variant_catalog_file} \
      --output-prefix ~{sample_id}

  >>>
  
  runtime {
    docker: expansion_hunter_docker
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 30
  }

}

task AnnotateExpansionHunter {
  input {
    String sample_id
    File repeats_file
    String expansion_hunter_docker
    File expansion_hunter_vcf
  }
  
  output {
    File expansion_hunter_vcf_annotated = "~{sample_id}.annotated.vcf"
  }
  
  command <<<

    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    # annotated_vcf = "${sample_id}.annotated.vcf"

    echo "[ RUNNING ] expansion hunter vcf annotation on sample ~{sample_id}"
    stranger \
      --repeats-file ~{repeats_file} \
      ~{expansion_hunter_vcf} > ~{sample_id}.annotated.vcf

  >>>

  runtime {
    docker: expansion_hunter_docker
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 30
  }

}