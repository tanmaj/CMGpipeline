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
            -e19 /cmg1scratch/cromwell/VariantAnnotationPipeline_databases/Exomizer/2302_hg19/2302_hg19_variants.mv.db \
            --assembly hg19 \
            -p "HP:0004322, HP:0000316, HP:0000369, HP:0000508,  HP:0000637,  HP:0006313, HP:0000470" \
            --vcf /cmg1scratch/cromwell/PX12139.vcf \
            --sample-id PX12139 \
            -o /cmg1scratch/cromwell/ \
            -f html -f json -f tsv
    }

    runtime {
        docker: "alesmaver/lirical"

    }

    output {
        File lirical_tsv = optitype_name + ".optitype_result.tsv"
        File lirical_json = optitype_name + ".optitype_result.tsv"
        File lirical_html = sampleName + ".DeepVariant.visual_report.html"
    }

}
