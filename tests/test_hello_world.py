import io
from contextlib import redirect_stdout

from cgr_gwas_qc.cli import hello_cli

def test_hello_cgr():
    expected = "Hello there pytest.\n"

    # Capture print statement
    observed = io.StringIO()
    with redirect_stdout(observed):
        hello_cli.print_hello("pytest")

    assert expected == observed.getvalue()
