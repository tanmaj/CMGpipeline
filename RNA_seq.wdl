version 1.0

import "https://raw.githubusercontent.com/AlesMaver/bioWDL_RNAseq/develop/align-hisat2.wdl" as AlignHisat2
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/GenerateCRAM.wdl" as GenerateCRAM

workflow RNA_seq {
  input {
        String sample_basename

        # data for AlignHisat2:
        Array[FastqPair]+ inputReads 
        String outputDir = "."
        String library
        Array[String] readgroups
        Array[File]+ indexFiles

        # references for CRAM generation:
        File reference_fa
        File reference_fai
  }

  call AlignHisat2.AlignHisat2 as AlignHisat2 {
      input:
        inputReads = inputReads, 
        outputDir = outputDir,
        sample = sample_basename,
        library = library,
        readgroups = readgroups,
        ### platform = "illumina",
        indexFiles = indexFiles
  }

  call GenerateCRAM.GenerateCRAM as GenerateCRAM {
      input:       
        input_bam = AlignHisat2.bamFile["file"],
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        sample_basename = sample_basename
  }

  output {
    File output_bam = AlignHisat2.bamFile["file"] 
    File output_bam_index = AlignHisat2.bamFile["index"]
    
    File output_cram = GenerateCRAM.output_cram 
    File output_cram_index = GenerateCRAM.output_cram_index 
  }
}
