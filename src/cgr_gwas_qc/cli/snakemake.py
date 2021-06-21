from typing import MutableSequence, Optional

import snakemake
import typer

from cgr_gwas_qc.config import ConfigMgr

context_settings = {
    "allow_extra_args": True,
    "ignore_unknown_options": True,
    "help_option_names": [],
}
app = typer.Typer(add_completion=False, context_settings=context_settings)


@app.command()
def main(ctx: typer.Context):
    """A light weight wrapper around snakemake.

    The purpose of this wrapper is to allow interacting with snakemake directly,
    while hiding some of the implementation details. The key feature of this
    wrapper is to tell snakemake where the workflow's snakefile exists.  We do
    this by adding the `-s` option to your snakemake command::

        snakemake -s <path to workflow install location> OPTIONS TARGETS

    In addition, we also add some snakemake options if you did not specify them.
    This is mostly for convenience.  Specifically, the we require the use of
    ``conda``, so we always will add ``--use-conda``.  Second, in recent
    versions of snakemake they no longer default to ``--cores 1`` and will throw
    an error if you did not specify a value.  We find this annoying, so we will
    always add ``--cores 1`` unless you provide a value.

    So for example, to run the full workflow locally you will typically run::

        cgr snakemake --cores 8 -k

    This will be translated into::

        snakemake -s <path to workflow install location> --cores 8 -k --use-conda

    Instead of running the entire workflow, you can also run specific
    sub-workflow by running::

        cgr snakemake --cores 8 -k --subworkflow sample_qc

    Args:
        subworkflow (str):
            Specify which sub-workflow to run [default: None].
            Must be one of [``entry_points`` | ``contamination`` | ``sample_qc`` | ``subject_qc`` | ``delivery``].
        **kwargs:
            All snakemake arguments will be passed to snakemake.
            See their documentation for possible options.

    References:
        - https://snakemake.readthedocs.io/en/stable/

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

    if not is_arg(args, ["--use-conda"]):
        args.append("--use-conda")

    if is_arg(args, ["--notemp", "--nt"]):
        # Adds config["notemp"] to help custom scripts know if they should
        # delete temp files.
        check_and_prepend_arg(args, ["--config"], "notemp=True")

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


if __name__ == "__main__":
    app()

typer_click_object = typer.main.get_command(app)  # only needed for building documentation
