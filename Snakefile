from os.path import exists, abspath, dirname
from snakemake import WorkflowError
from warnings import warn

def addwd(target, rule, root=False, parent=False):
    """Getcwd within a rule. Will be moved to snakecharmer in the future."""
    if root:
        wd = abspath(dirname(workflow.snakefile))
    else:
        wd = abspath(dirname(getattr(rules, rule).snakefile))
    if parent:
        wd = dirname(wd)
    setattr(target, "wd", wd)

include: "feature/anno.snake"
include: "feature/alignments.snake"
include: "feature/fft.snake"

try:
    if exists("sandbox.snake"):
        include: "sandbox.snake"
except WorkflowError:
    warn("Couldn't include sandbox.snake, but trying to proceed")
