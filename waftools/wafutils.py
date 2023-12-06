"""
This modules defines some utilities shared by several *wscript* files.
"""

import os
import os.path as osp
from subprocess import check_output

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
    if env.BUILD_MPI and not cmdprefix and env["CONFIG_PARAMETERS"]["require_mpiexec"]:
        cmdprefix = env["base_mpiexec"] + ["-n", "1"]
    cmds = " ".join(cmdprefix + cmdexe)
    Logs.debug("os environ: %r" % environ)
    if kwargs.get("for_catalo"):
        del kwargs["for_catalo"]
        # do not confuse with installed elementsdir
        environ["ASTER_ELEMENTSDIR"] = ""
    try:
        if Logs.verbose:
            kwargs["env"] = environ
            kwargs["shell"] = True
            self.log_command(cmds, kwargs)
            ret, stdout, stderr = Utils.run_process(cmds, kwargs)
            if ret:
                error = Errors.WafError("Command %r returned %r" % (cmds, ret))
                error.returncode = ret
                error.stderr = stderr
                error.stdout = stdout
                raise error
        else:
            kwargs["output"] = Context.BOTH
            stdout, stderr = self.cmd_and_log(cmds, env=environ, shell=True, quiet=0, **kwargs)
        with open(logbase + ".out", "w") as flog:
            flog.write(stdout or "")
        with open(logbase + ".err", "w") as flog:
            flog.write(stderr or "")
    except Errors.WafError as err:
        Logs.warn("stdout: %s" % err.stdout)
        Logs.warn("stderr: %s" % err.stderr)
        Logs.info("To run manually, use:")
        Logs.info('. "%s/profile.sh"' % env["ASTERDATADIR"])
        Logs.info(" ".join(cmdprefix + cmdexe))
        raise


def add_to_env_paths(environ, name, path):
    if not hasattr(add_to_env_paths, "pathsep"):
        add_to_env_paths.pathsep = os.pathsep
        # pathsep should be the one returned by the defined python exe interpreter
        # that may differ from the python used to run waf in case of cross compiling
        if "PYTHON" in environ.keys():
            add_to_env_paths.pathsep = check_output(
                [environ["PYTHON"], "-c", "import os; print(os.pathsep, end='')"], encoding="utf-8"
            )

    if not path:
        return
    paths = [path] if isinstance(path, str) else path
    raw = environ.get(name, None)
    if raw is not None:
        paths += raw.split(add_to_env_paths.pathsep)
    environ[name] = add_to_env_paths.pathsep.join(paths)


def remove_previous(install_node, patterns):
    """Remove previously installed files (that are just copied to dest)."""
    if not install_node:
        return

    for pattern in patterns:
        for i in install_node.ant_glob(pattern):
            os.remove(i.abspath())
