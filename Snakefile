include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "feature/fft.snake"

rule test:
    input: "data/fft_summaries/GRCh37.75-mini/fft_IH03_WPS.tsv.gz"
