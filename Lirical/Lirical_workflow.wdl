version 1.0

workflow Lirical_workflow {
    input {
        String sample_basename
        String hpo_ids
        File hg19_variants_mv_db
        File input_vcf

    }

    call Lirical {
        input:
            sample_basename = sample_basename,
            hpo_ids = hpo_ids,
            hg19_variants_mv_db = hg19_variants_mv_db,
            input_vcf = input_vcf
    }

    output {
        File lirical_tsv = Lirical.lirical_tsv
        File lirical_json = Lirical.lirical_json
        File lirical_html = Lirical.lirical_html

    }

task Lirical {
    input {
        String sample_basename
        String hpo_ids
        File hg19_variants_mv_db
        File input_vcf
    }

    command {
        set -e

        java -jar lirical-cli-2.0.0-RC2/lirical-cli-2.0.0-RC2.jar \
            prioritize \
            -d ./data/ \
            -e19 ~{hg19_variants_mv_db} \
            --assembly hg19 \
            -p ~{hpo_ids} \
            --vcf ~{input_vcf} \
            --sample-id ~{sample_basename} \
            -o ./ \
            -f html -f json -f tsv

        mv lirical.tsv  ~{sample_basename}.lirical.tsv
        mv lirical.json  ~{sample_basename}.lirical.json
        mv lirical.html  ~{sample_basename}.lirical.html
    }

    runtime {
        docker: "alesmaver/lirical"

    }

    output {
        File lirical_tsv = sample_basename + ".lirical.tsv"
        File lirical_json = sample_basename + ".lirical.json"
        File lirical_html = sample_basename + ".lirical.html"
    }

}
