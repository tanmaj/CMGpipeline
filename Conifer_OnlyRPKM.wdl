version 1.0
## Copyright CMG@KIGM, Ales Maver
import "./Conifer.wdl" as Conifer

# WORKFLOW DEFINITION 
workflow Conifer.Conifer_OnlyRPKM as Conifer_OnlyRPKM {
  input {
    File input_bam
    File input_bam_index

    String? enrichment
    File? enrichment_bed
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

  output {
    File output_rpkm = MakeRPKM.output_rpkm
  }
}
