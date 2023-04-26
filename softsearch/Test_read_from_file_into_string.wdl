version 1.0

# WORKFLOW DEFINITION 
workflow TestWF {
  input {
    String sample_basename
  }
  
  call FileToArray
  
  scatter (chromosome in FileToArray.scatter_regions ) {
    call DoSomethingWithChromosome {
      input:
        chromosome=chromosome,
        sample_basename=sample_basename
    }

  }

}

task FileToArray {
    command <<<
        # Get chromosome interval list file
        wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/wgs_calling_regions.v1_mod.list.txt
        # Re-attempt the wget without the https_proxy set
        unset https_proxy 
        wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/wgs_calling_regions.v1_mod.list.txt
        
        var=$(<file)
        echo $var
    >>>
      
    runtime {
        docker: "alesmaver/softsearch"
        runtime_minutes: 5
    }
    output {
        Array[String] scatter_regions = read_lines(stdout())
    }
}


task DoSomethingWithChromosome {
  input {
    String chromosome
    String sample_basename
  }

  command <<<
    echo sample_basename: ~{sample_basename}
    echo chromosome: ~{chromosome}
  >>>

  runtime {
    docker: "alesmaver/softsearch"
    runtime_minutes: 5
  }
  output {
    File outfile = stdout()
  }
}

