version 1.0
## Copyright CMG@KIGM, Ales Maver

## Usage
# This workflow accepts three types of inputs: Illumina FASTQ files, a BAM file, or CRAM files
# Currently, the input CRAM files should be aligned to the hg19 reference genome assembly, we will implement support for other genome formats in the future
# The CRAM output is optional and disabled by default at the moment, until production switches to CRAM

# Subworkflows
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/MitoMap.wdl" as MitoMap

# WORKFLOW DEFINITION 
workflow MitoMapWorkflow {
  input {
    File input_vcf
    String sample_basename

    File reference_fa
    File reference_fai
    File reference_dict
  }  

  call MitoMap.CreateMitoFasta as CreateMitoFasta {
    input:
    input_vcf = SelectFinalVariants.output_vcf,
    sample_basename = sample_basename,

    reference_fa = reference_fa,
    reference_fai = reference_fai,
    reference_dict = reference_dict,

    docker = "broadinstitute/gatk3:3.8-1"
  }

  call MitoMap.MitoMap as MitoMap {
    input:
    mtDNA_fasta = CreateMitoFasta.mtDNA_fasta,
    sample_basename = sample_basename
  }

  output {
    File? mitoResults_xls = MitoMap.mitoResults_xls
    File? mitoResults_txt = MitoMap.mitoResults_txt

  }
}

