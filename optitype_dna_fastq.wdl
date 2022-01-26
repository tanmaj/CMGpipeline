version 1.0

task OptitypeDnafromFastq {
  input {
    String optitype_name = "optitype"
    File input_fq1
    File input_fq2
  }

  Int space_needed_gb = 10 + round(5*size([input_fq1, input_fq2], "GB"))
  runtime {
    memory: "64GB"
    docker: "mgibio/immuno_tools-cwl:1.0.1"
    disks: "local-disk ~{space_needed_gb} SSD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    /usr/bin/python /usr/local/bin/OptiType/OptiTypePipeline.py -i ~{input_fq1} ~{input_fq2} --dna -v -p ~{optitype_name} -o . 
    mv ~{optitype_name}_result.tsv ~{optitype_name}.optitype_result.tsv
    mv ~{optitype_name}_coverage_plot.pdf ~{optitype_name}.optitype_coverage_plot.pdf
  >>>

  output {
    File optitype_tsv = optitype_name + "_result.tsv"
    File optitype_plot = optitype_name + "_coverage_plot.pdf"
  }
}

workflow OptitypeDnafromFastq {
  input {
    String? optitype_name
    File input_fq1
    File input_fq2
  }
  call OptitypeDnafromFastq {
    input:
    optitype_name=optitype_name,
    input_fq1 = input_fq1,
    input_fq2 = input_fq2,
  }
}
