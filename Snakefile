rule anno_upstream:
    input: "{prefix}-{ref_id}.tsv.gz"
    output: "{prefix}-{ref_id,[\w\.]+}.upstream.tsv"
    shell: """
        zcat {input} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{\
            if ($5 == "+") {{print $1, $2, $3-10000, $3, $5}} \
            else {{print $1, $2, $4, $4+10000, $5}} \
        }}' > {output}
    """

rule anno_downstream:
    input: "{prefix}-{ref_id}.tsv.gz"
    output: "{prefix}-{ref_id,[\w\.]+}.downstream.tsv"
    shell: """
        zcat {input} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{\
            if ($5 == "+") {{print $1, $2, $4, $4+10000, $5}} \
            else {{print $1, $2, $3-10000, $3, $5}} \
        }}' > {output}
    """

rule anno_body:
    input: "{prefix}-{ref_id}.tsv.gz"
    output: "{prefix}-{ref_id,[\w\.]+}.body.tsv"
    shell: """
        zcat {input} \
        | awk 'BEGIN {{FS="\t"; OFS="\t"}} {{\
            if ($5 == "+") {{print $1, $2, $3-1, $3-1+10000, $5}} \
            else {{print $1, $2, $4-1-10000, $4-1, $5}} \
        }}' > {output}
    """

rule bam_index:
    input: "{prefix}.bam"
    output: "{prefix}.bam.bai"
    shell: "samtools index {input} {output}"

rule wps:
    input:
        bam="{path}/alignments/{name}.bam",
        bai="{path}/alignments/{name}.bam.bai",
        anno="expression/transcriptAnno-GRCh37.75.body.tsv"
    params: out_dir="{path}/wps/{name}"
    output: dynamic("{path}/wps/{name}/block_{id}.tsv.gz")
    shell: """
        ./expression/extractReadStartsFromBAM_Region_WPS.py \
            --minInsert=120 --maxInsert=180 \
            -i {input.anno} \
            -o '{params.out_dir}/block_%s.tsv.gz' \
            {input.bam}
    """

rule all:
    input: rules.wps.output
    output: "{path}/{name}.plc"
    shell: "touch {output}"
