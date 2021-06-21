version 1.0

# MIT License
#
# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/bcftools.wdl" as bcftools
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/bwa.wdl" as bwa
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/clever.wdl" as clever
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/common.wdl" as common
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/delly.wdl" as delly
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/duphold.wdl" as duphold
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/gridss.wdl" as gridss
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/manta/manta.wdl" as manta
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/picard.wdl" as picard
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/samtools.wdl" as samtools
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/smoove.wdl" as smoove
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/survivor.wdl" as survivor

workflow SVcalling {
    input {
        File bamFile
        File bamIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        BwaIndex bwaIndex
        String sample
        String newId = "\'%CHROM\\_%POS\'"
        Boolean excludeMisHomRef = false
        Boolean excludeFpDupDel = false
        String outputDir = "."

        Array[File] input_manta_reference_vcfs

        Map[String, String] dockerImages = {
            "bcftools": "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2",
            "clever": "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6",
            "delly": "quay.io/biocontainers/delly:0.8.1--h4037b6b_1",
            "gridss": "quay.io/biocontainers/gridss:2.9.4--0",
            "manta": "quay.io/biocontainers/manta:1.4.0--py27_1",
            "picard":"quay.io/biocontainers/picard:2.23.2--0",
            "samtools": "quay.io/biocontainers/samtools:1.10--h9402c20_2",
            "survivor": "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0",
            "smoove": "quay.io/biocontainers/smoove:0.2.5--0",
            "duphold": "quay.io/biocontainers/duphold:0.2.1--h516909a_1"
        }
    }

    meta {allowNestedInputs: true}

    String SVdir = outputDir + '/structural-variants/'

    call manta.Germline as manta {
        input:
            dockerImage = dockerImages["manta"],
            bamFile = bamFile,
            bamIndex = bamIndex,
            sample = sample,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            runDir = SVdir + 'manta/'
    }

    call MergeMantaFiles {
        input:
            input_manta_vcf = manta.mantaVCF,
            input_manta_reference_vcfs = input_manta_reference_vcfs,
            sample_basename = sample
    }

    call AnnotateMantaVCF {
        input:
            input_vcf = MergeMantaFiles.merged_vcf,
            sample_basename = sample
    }

    output {
        File mantaVcf = manta.mantaVCF
        File mantaVcfindex = manta.mantaVCFindex
        File output_sv_table = AnnotateMantaVCF.output_sv_table
    }

    parameter_meta {
        outputDir: {description: "The directory the output should be written to.", category: "common"}
        referenceFasta: { description: "The reference fasta file", category: "required" }
        referenceFastaFai: { description: "Fasta index (.fai) file of the reference", category: "required" }
        referenceFastaDict: { description: "Sequence dictionary (.dict) file of the reference", category: "required" }
        bamFile: {description: "sorted BAM file", category: "required"}
        bamIndex: {description: "BAM index(.bai) file", category: "required"}
        bwaIndex: {description: "Struct containing the BWA reference files", category: "required"}
        sample: {description: "The name of the sample", category: "required"}
        newId: {description: "Assign ID on the fly (e.g. --set-id +'%CHROM\_%POS').", category: "advanced"}
        excludeMisHomRef: {description: "Option to exclude missing and homozygous reference genotypes.", category: "advanced"}
        excludeFpDupDel: {description: "Option to exclude false positive duplications and deletions according to DUPHOLD.", category: "advanced"}
        dockerImages: {description: "A map describing the docker image used for the tasks.",
                           category: "advanced"}
    }
}

###
###
###
task MergeMantaFiles {
  input {
    # Command parameters
    File input_manta_vcf
    Array[File] input_manta_reference_vcfs 
    String sample_basename
  }

  command <<<
  MANTA_VCFS_DIR=$(dirname ~{input_manta_reference_vcfs[0]})
  cp ~{input_manta_vcf} $MANTA_VCFS_DIR
  cp $MANTA_VCFS_DIR/*.vcf ./

  gunzip ./*manta.vcf.gz
  for file in *manta.vcf; do
    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > "${file/vcf/filtered.vcf}"
  done

  ls -d ./*manta.filtered.vcf > fileList

  SURVIVOR merge fileList 500 1 1 0 0 20 ./merged.vcf
  >>>

  runtime {
    docker: "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
    requested_memory_mb_per_core: 3000
    cpu: 2
    runtime_minutes: 120
    docker_user: "root"
  }
  output {
    File merged_vcf = "merged.vcf"
  }
}

# SNPeff task
task AnnotateMantaVCF {
  input {
    # Command parameters
    File input_vcf
    String sample_basename
  }

  command { 
  java -jar /home/biodocker/bin/snpEff/snpEff.jar -noInteraction -noDownload hg19 ~{input_vcf} > ~{sample_basename}.merged.annotated.vcf
  java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'MODERATE' | ANN[*].IMPACT has 'HIGH') & SUPP < 2 & isVariant(GEN[PX7084].GT)" ~{sample_basename}.merged.annotated.vcf > ~{sample_basename}.merged.annotated.filtered.vcf
  java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields -s "," -e "." ~{sample_basename}.merged.annotated.filtered.vcf CHROM POS REF ALT SVLEN QUAL ANN[*].GENE ANN[1].IMPACT ANN[1].EFFECT GEN[*].GT > ~{sample_basename}.mantaSVs.txt

  }
  runtime {
    docker: "alesmaver/snpeff_v50"
    maxRetries: 3
    requested_memory_mb_per_core: 2000
    cpu: 8
    runtime_minutes: 120
  }
  output {
    File output_sv_table = "~{sample_basename}.mantaSVs.txt"
  }
}
