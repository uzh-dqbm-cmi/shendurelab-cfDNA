include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "feature/fft.snake"

rule test:
    input: "data/fft_summaries/GRCh37.75-mini/fft_IH03_WPS.tsv.gz"

rule test_reset:
    shell: """
        rm -f data/fft_summaries/GRCh37.75-mini/*
        rm -f data/individuals/IH03/fft/GRCh37.75-mini.tar
        rm -f data/individuals/IH03/wps/GRCh37.75-mini.tar
    """
