from os.path import exists
from snakemake import WorkflowError
from warnings import warn

include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "feature/fft.snake"

try:
    if exists("sandbox.snake"):
        include: "sandbox.snake"
except WorkflowError:
    warn("Couldn't include sandbox.snake, but trying to proceed")
