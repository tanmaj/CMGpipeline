version 1.0

workflow CreateGenesBed {
    call DownloadAndPrepareBed
}

task DownloadAndPrepareBed {
    command <<<
        curl -s ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz \
            | gunzip --stdout - \
            | awk '$3 == "gene"' \
            | convert2bed -i gff --attribute-key="gene_name" - \
            > genes.bed
    >>>

    output {
        File genes_bed = "genes.bed"
    }

    runtime {
        docker: "dceoy/bedops"
        cpu: 1
        memory: "4G"
    }

}
