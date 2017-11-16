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
