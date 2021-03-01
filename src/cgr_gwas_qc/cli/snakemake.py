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
    """
    args = ctx.args
    if is_arg(args, ["-h", "--help"]):
        print(main.__doc__)
        print("\nSnakemake help is below:\n")
        snakemake.main(args)

    check_and_prepend_arg(args, ["-s", "--snakefile"], ConfigMgr.SNAKEFILE.as_posix())

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
