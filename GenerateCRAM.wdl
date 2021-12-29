version 1.0

import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/FastqToVCFPipeline_3.wdl" as FastqToVcf

# WORKFLOW DEFINITION 
workflow GenerateCRAM {
  input {
    File input_bam
    File reference_fa
    File reference_fai
    String sample_basename
  }
  
call FastqToVcf.ConvertToCram as ConvertToCram {
    input:
      input_bam = input_bam,
      ref_fasta = reference_fa,
      ref_fasta_index = reference_fai,
      sample_basename = sample_basename
  }
  
  
  output {
    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
  }
  
  # end-workflow
}
