from os.path import exists

include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "feature/fft.snake"

if exists("sandbox.snake"):
    include: "sandbox.snake"
