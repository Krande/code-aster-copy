# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------

import os
from pathlib import Path
from subprocess import CalledProcessError, check_call

import click

CTXT = dict(allow_extra_args=True, ignore_unknown_options=True)

# build context
ROOT = Path(__file__).parent.parent
BUILD = ROOT / "build"
ASTER_BUILD = os.environ.get("ASTER_BUILD", "mpi")


C4CHE = BUILD / "mpi" / "c4che" / "release_cache.py"
WAF = ROOT / "waf_mpi"
if ASTER_BUILD == "debug":
    C4CHE = BUILD / "mpidebug" / "c4che" / "debug_cache.py"
    WAF = ROOT / "waf_debug"


# executed under commands to have a beautiful error message
def _check():
    if ASTER_BUILD not in ("mpi", "debug", "std"):
        raise click.UsageError(
            f"ASTER_BUILD: expecting 'mpi', 'debug' or 'std', not: {ASTER_BUILD}"
        )


def call_command(cmd: list[str]):
    """Wrap `subprocess.call`."""
    print("COMMAND:", *cmd)
    try:
        check_call(cmd)
    except CalledProcessError as exc:
        raise click.ClickException(str(exc))


@click.command(context_settings=CTXT)
@click.option("--prefix", metavar="PATH")
@click.pass_context
def configure(ctx, prefix):
    """
    Configure the project
    """
    _check()
    args = list(ctx.args)
    if prefix:
        prefix = Path(prefix).absolute()
        args = ["--prefix", str(prefix)] + args
    call_command([WAF, "configure"] + list(args))


@click.command(context_settings=CTXT)
@click.pass_context
def build(ctx):
    """
    Build the project
    """
    _check()
    if not C4CHE.exists():
        ctx.invoke(configure)
    args = list(ctx.args)
    call_command([WAF, "build"] + args)


@click.command(context_settings=CTXT)
@click.pass_context
def install(ctx):
    """
    Build and install the project
    """
    _check()
    if not C4CHE.exists():
        ctx.invoke(configure)
    args = list(ctx.args)
    call_command([WAF, "install"] + args)


@click.command(context_settings=CTXT)
@click.pass_context
def doc(ctx):
    """
    Build the embedded documentation
    """
    _check()
    if ASTER_BUILD == "std":
        raise click.UsageError("doc skipped, only available in parallel")
    args = list(ctx.args)
    call_command([WAF, "doc"] + args)


@click.command(context_settings=CTXT)
@click.pass_context
def test(ctx):
    """
    Run one or more testcases

    All arguments will be treated as testcases names.
    """
    if not ctx.args:
        raise click.UsageError("no testcase name provided")
    args = []
    for testname in ctx.args:
        args.extend(["-n", testname])
    call_command([WAF, "test"] + args)


@click.command()
@click.pass_context
def bootstrap(ctx):
    """
    Run all steps: configure, install and doc
    """
    ctx.invoke(configure)
    ctx.invoke(install)
    ctx.invoke(doc)


@click.command()
def distclean():
    """
    Run distclean.
    """
    call_command([WAF, "distclean"])
