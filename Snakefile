from os.path import exists, abspath, dirname
from snakemake import WorkflowError
from warnings import warn

def addwd(target, rule):
    """Getcwd within a rule. Will be moved to snakecharmer in the future."""
    wd = abspath(dirname(getattr(rules, rule).snakefile))
    setattr(target, "wd", wd)

include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "feature/fft.snake"

try:
    if exists("sandbox.snake"):
        include: "sandbox.snake"
except WorkflowError:
    warn("Couldn't include sandbox.snake, but trying to proceed")
