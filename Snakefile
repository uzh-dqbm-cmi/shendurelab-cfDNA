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
        anno="expression/transcriptAnno-GRCh37.75.mini.tsv"
    output: "{path}/wps/{name}/blocks.tar"
    params: out_dir="{path}/wps/{name}"
    shell: """
        ./expression/extractReadStartsFromBAM_Region_WPS.py \
            --minInsert=120 --maxInsert=180 \
            -i {input.anno} \
            -o '{params.out_dir}/block_%s.tsv.gz' \
            {input.bam};
        cd {params.out_dir} && tar cf blocks.tar block_*.tsv.gz && rm -f block_*.tsv.gz;
    """

rule fft:
    input: "{path}/wps/{name}/blocks.tar"
    output: "{path}/fft/{name}/blocks.tar"
    params: in_dir="{path}/wps/{name}", out_dir="{path}/fft/{name}"
    shell: """
        cwd=$(pwd);
        cd {params.in_dir} && tar xf blocks.tar;
        cd $cwd && (cd {params.in_dir}; \ls block_*.tsv.gz) | xargs -n 500 \
            Rscript ./expression/fft_path.R {params.in_dir} {params.out_dir};
        cd $cwd/{params.out_dir} && tar cf blocks.tar block_*.tsv.gz && rm -f block_*.tsv.gz;
        cd $cwd/{params.in_dir} && rm -f block_*.tsv.gz;
    """

rule fft_summaries:
    input: rules.fft.output
    output:
        "{path}/fft_summaries/fft_{name}_WPS.tsv.gz",
        "{path}/fft_summaries/fft_{name}_cov.tsv.gz",
        "{path}/fft_summaries/fft_{name}_starts.tsv.gz"
    params:
        in_dir="{path}/fft/{name}",
        out_dir="{path}/fft_summaries",
        anno="expression/transcriptAnno-GRCh37.75.mini.tsv"
    shell: """
        cwd=$(pwd);
        cd {params.in_dir} && tar xf blocks.tar;
        cd $cwd && ./expression/convert_files.py \
            -a {params.anno} \
            -t {wildcards.path} \
            -r {wildcards.path} \
            -p . \
            -i {wildcards.name};
        cd $cwd/{params.in_dir} && rm -f block_*.tsv.gz;
    """
