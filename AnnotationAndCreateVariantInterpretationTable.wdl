version 1.0

# Subworkflows
import "./AnnotationPipeline.wdl" as Annotation
import "./CreateInterpretationTable.wdl" as CreateInterpretationTable

workflow AnnotateAndTable {
  input {
    File input_variant_string

    String sample_basename

    File chromosome_list
    
    File gnomAD_vcf
    File gnomAD_vcf_index

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

    File SLOpopulation_vcf
    File SLOpopulation_vcf_index

    File ClinVar_vcf
    File ClinVar_vcf_index

    File SpliceAI
    File SpliceAI_index

    File dbscSNV
    File dbscSNV_index

    File HPO
    File HPO_index
    File OMIM
    File OMIM_index
    File gnomadConstraints
    File gnomadConstraints_index
    File CGD
    File CGD_index
    File bcftools_annotation_header

    File reference_fa
    File reference_fai
    File reference_dict

    File dbNSFP
    File dbNSFP_index

    String? targetRegions

    #String bgzip_docker = "dockerbiotools/bcftools:latest"
    String bcftools_docker = "dceoy/bcftools:latest"
    #String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
    String SnpEff_docker = "alesmaver/snpeff_v50:latest"
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    String vcfanno_docker = "clinicalgenomics/vcfanno:0.3.2"
    String R_docker = "alesmaver/r-base"
  }

  call CreateVCFfromString {
      input:
        input_variant_string = input_variant_string,
        sample_basename = sample_basename,

        docker = R_docker
    }

  call GetVariantRegions {
      input:
        input_variant_string = input_variant_string,
        separator = ":"
    }

  call Annotation.AnnotateVCF as AnnotateVCF {
    input:
      input_vcf = CreateVCFfromString.output_vcf,
      chromosome_list = chromosome_list,
      
      gnomAD_vcf = gnomAD_vcf,
      gnomAD_vcf_index = gnomAD_vcf_index,

      gnomADexomes_vcf = gnomADexomes_vcf,
      gnomADexomes_vcf_index = gnomADexomes_vcf_index,

      SLOpopulation_vcf = SLOpopulation_vcf,
      SLOpopulation_vcf_index = SLOpopulation_vcf_index,

      ClinVar_vcf = ClinVar_vcf,
      ClinVar_vcf_index = ClinVar_vcf_index,

      SpliceAI = SpliceAI,
      SpliceAI_index = SpliceAI_index,

      dbscSNV = dbscSNV,
      dbscSNV_index = dbscSNV_index,

      HPO = HPO,
      HPO_index = HPO_index,
      OMIM = OMIM,
      OMIM_index = OMIM_index,
      gnomadConstraints = gnomadConstraints,
      gnomadConstraints_index = gnomadConstraints_index,
      CGD = CGD,
      CGD_index = CGD_index,
      bcftools_annotation_header = bcftools_annotation_header,

      fasta_reference = reference_fa,
      fasta_reference_index = reference_fai,
      fasta_reference_dict = reference_dict,

      dbNSFP = dbNSFP,
      dbNSFP_index = dbNSFP_index,

      targetRegions = select_first([GetVariantRegions.regions, targetRegions]),

      #bcftools_docker = bcftools_docker,
      #SnpEff_docker = SnpEff_docker,
      gatk_docker = gatk_docker,
      gatk_path = gatk_path,
      vcfanno_docker = vcfanno_docker
  }

  # Get snpEff and dbNSFP annotations
  call GetGeneListFromVCF {
    input:
      input_vcf = AnnotateVCF.output_vcf,
      docker = SnpEff_docker
    }

  call CreateInterpretationTable.CreateInterpretationTable as CreateInterpretationTable {
    input:
      input_vcf = AnnotateVCF.output_vcf,
      input_vcf_index = AnnotateVCF.output_vcf_index,
      panel_gene_list = GetGeneListFromVCF.panel_gene_list
  }

  # Outputs that will be retained when execution is complete
  output {
    File XLSX_OUTPUT = CreateInterpretationTable.XLSX_OUTPUT
  }
}

task GetGeneListFromVCF {
	input { 
		File input_vcf

		String docker
	}

	command {
		SNPSIFT_EXTRACTFIELDS='java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields -s "," - "ANN[*].GENE" | uniq '
		zcat ~{input_vcf} | eval $SNPSIFT_EXTRACTFIELDS | tail -n +2 | paste -sd "," - 
	}

	runtime {
	    docker: docker
	    requested_memory_mb_per_core: 3000
	    cpu: 1
	    runtime_minutes: 180
	}
	output {
		String panel_gene_list = read_string(stdout())
	}
}

# Output VCF from string
task CreateVCFfromString {
  input {
    # Command parameters
    String input_variant_string
    String CreateVCFfromString_Rscript = "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/R_scripts/SCRIPTS_createVCFfromVariantString.R"
    String sample_basename

    # Runtime parameters
    String docker
  }

  command {
  wget -t 1 -T 20 ~{CreateVCFfromString_Rscript}
  # Repeat in case the proxy defined in the docker image would case problems accessing the GitHub repo
  unset https_proxy
  wget -t 1 -T 20 ~{CreateVCFfromString_Rscript}

  Rscript SCRIPTS_createVCFfromVariantString.R --variant="~{input_variant_string}" --sample="~{sample_basename}"    
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 30
  }
  output {
    File output_vcf = "~{sample_basename}.IMPORTVARIANT.vcf"
  }
}

task GetVariantRegions {
  input {
    String input_variant_string
    String separator
  }
  command <<<
    echo "~{input_variant_string}" | tr "~{separator}" \\n | awk -F' ' '{print $1":"$2-1"-"$2+1}' OFS='\t' | tr \\n ";" | sed 's/\;$//'
  >>>
  runtime {
    docker:"biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 500
    cpu: 1
    runtime_minutes: 5
  }
  output {
    Array[String] value = read_lines(stdout())
    String regions = value[0]
  }
}
