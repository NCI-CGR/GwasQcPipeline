from typing import MutableSequence, Optional

import snakemake
import typer

from cgr_gwas_qc.config import ConfigMgr


def fmt_doc_string(func):
    func.__doc__ = func.__doc__.format(ConfigMgr.SNAKEFILE)
    return func


@fmt_doc_string
def main(ctx: typer.Context):
    """Run the Gwas Qc Pipeline using Snakemake.

    This is a wrapper around snakemake. It will pre-pend the `-s` option:

        snakemake -s {} OPTIONS TARGETS

    This is how the user should run snakemake locally.

    You can run specific sub-workflows by using the `--subworkflow WORKFLOW` option.

    Where possible values of WORKFLOW are:

    - entry_points

    """
    args = ctx.args
    if is_arg(args, ["-h", "--help"]):
        print(main.__doc__)
        print("\nSnakemake help is below:\n")
        snakemake.main(args)

    # If the user provided `--subworkflow` option, then run just the
    # subworkflow. Otherwise run the main workflow.
    workflow_name = pop_arg(args, "--subworkflow")
    if workflow_name:
        snakefile = ConfigMgr.subworkflow(workflow_name)
    else:
        snakefile = ConfigMgr.SNAKEFILE.as_posix()

    check_and_prepend_arg(args, ["-s", "--snakefile"], snakefile)
    if not is_arg(args, ["--profile"]):
        check_and_prepend_arg(args, ["-j", "--cores", "--jobs"], "1")

    snakemake.main(args)


def is_arg(args, options):
    """Check if an option is already in the argument list."""
    return any(option in args for option in options)


def check_and_prepend_arg(args, options, value):
    """Prepend an option to the front of the arg list."""
    if not is_arg(args, options):
        [args.insert(i, v) for i, v in enumerate([options[0], value])]


def pop_arg(args: MutableSequence, name: str) -> Optional[str]:
    if name not in args:
        return None

    arg_location = args.index(name)
    _, value = args.pop(arg_location), args.pop(arg_location)

    return value
