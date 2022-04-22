"""
This modules defines some utilities shared by several *wscript* files.
"""

import os
import os.path as osp

from waflib import Context, Errors, Logs, Utils


def exec_pyaster(self, pyfile, args, **kwargs):
    """Execute aster depending on the configuration"""
    cwd = kwargs["cwd"]
    logbase = osp.join(cwd, osp.basename(pyfile))
    env = self.all_envs[self.variant]
    if "env" in kwargs:
        environ = kwargs["env"]
        del kwargs["env"]
    else:
        environ = os.environ.copy()
    python = list(env.PYTHON)[0]

    python_ld_path = env.ASTERLIBDIR
    add_to_env_paths(environ, "PYTHONPATH", python_ld_path)
    add_to_env_paths(environ, "LD_LIBRARY_PATH", python_ld_path)

    cmdexe = [python, pyfile] + args
    # this position allows CATALO_CMD to define an environment variable
    # or a prefix to wrap the executable
    cmdprefix = Utils.to_list(env["CATALO_CMD"])
    if self.env.BUILD_MPI and not cmdprefix and env["CONFIG_PARAMETERS"]["require_mpiexec"]:
        cmdprefix = env["base_mpiexec"] + ["-n", "1"]
    cmds = " ".join(cmdprefix + cmdexe)
    Logs.debug("os environ: %r" % environ)
    # do not confuse with installed elementsdir
    environ["ASTER_ELEMENTSDIR"] = ""
    kwargs["output"] = Context.BOTH
    try:
        stdout, stderr = self.cmd_and_log(cmds, env=environ, shell=True, quiet=0, **kwargs)
        with open(logbase + ".out", "w") as flog:
            flog.write(stdout)
        with open(logbase + ".err", "w") as flog:
            flog.write(stderr)
    except Errors.WafError as err:
        Logs.warn("stdout: %s" % err.stdout)
        Logs.warn("stderr: %s" % err.stderr)
        Logs.info("To run manually, use:")
        Logs.info('export LD_LIBRARY_PATH="%s"' % environ["LD_LIBRARY_PATH"])
        Logs.info('export PYTHONPATH="%s"' % environ["PYTHONPATH"])
        Logs.info(" ".join(cmdprefix + cmdexe))
        raise


def add_to_env_paths(environ, name, path):
    if not path:
        return
    paths = [path] if isinstance(path, str) else path
    raw = environ.get(name, None)
    if raw is not None:
        paths += raw.split(os.pathsep)
    environ[name] = os.pathsep.join(p for p in paths)


def remove_previous(install_node, patterns):
    """Remove previously installed files (that are just copied to dest)."""
    if not install_node:
        return

    for pattern in patterns:
        for i in install_node.ant_glob(pattern):
            os.remove(i.abspath())
