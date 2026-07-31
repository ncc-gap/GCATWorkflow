"""
Compatibility shim for the ``snakemake.snakemake()`` Python API.

snakemake >= 8 (installed on Python >= 3.11) removed the top-level
``snakemake.snakemake()`` function that these tests use for dry-run validation.
On Python 3.10 snakemake 7.x is installed and still provides it.

This shim keeps the original behaviour on snakemake 7 and falls back to the
command line interface on snakemake >= 8, so the tests work on Python
3.10/3.11/3.12.
"""
import sys
import subprocess
import importlib

_sm = importlib.import_module("snakemake")


def snakemake(snakefile, workdir=None, dryrun=False, printdag=False, **kwargs):
    # snakemake 7.x: use the original Python API
    if hasattr(_sm, "snakemake"):
        return _sm.snakemake(
            snakefile, workdir=workdir, dryrun=dryrun, printdag=printdag, **kwargs
        )
    # snakemake >= 8: the Python entry point was removed, use the CLI
    cmd = [sys.executable, "-m", "snakemake", "--snakefile", str(snakefile), "--cores", "1"]
    if workdir is not None:
        cmd += ["--directory", str(workdir)]
    if printdag:
        cmd += ["--dag"]
    elif dryrun:
        cmd += ["-n"]
    return subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
