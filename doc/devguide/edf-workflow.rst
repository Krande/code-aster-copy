.. _devguide-worflow_edf:

########################
EDF development workflow
########################


This section describes the adopted workflow for the continuous integration
in the EDF repository.


Rules for EDF continuous integration
====================================

- The main branch is ``main``.
  Developers must not commit revisions in this branch.

- Developers must use their own branch to commit changes: ``XY-feature`` (use
  your initials for XY).

- Commit messages must refer to the an existing issue.


Configure your repository
=========================

.. note:: In the examples, the repositories paths are
    :file:`$HOME/dev/codeaster/src` for *src* and
    :file:`$HOME/dev/codeaster/devtools` for *devtools*.

- Update the devtools repository:

  .. code-block:: sh

        cd $HOME/dev/codeaster/devtools && hg pull && hg update default

- Configure the URL of the remote repositories:

  .. code-block:: sh

        cd $HOME/dev/codeaster/src
        install_env
        # or
        $HOME/dev/codeaster/devtools/bin/install_env

  It should print something like:

  .. code-block:: none

        Configuring 'src' repository...
        INFO     checking repository $HOME/data/dev/codeaster/src
        INFO     settings added into $HOME/data/dev/codeaster/src/.hg/hgrc
        INFO     checking repository $HOME/data/dev/codeaster/data
        INFO     settings added into $HOME/data/dev/codeaster/data/.hg/hgrc
        INFO     checking repository $HOME/data/dev/codeaster/validation
        INFO     settings added into $HOME/data/dev/codeaster/validation/.hg/hgrc
        INFO     checking repository $HOME/data/dev/codeaster/devtools
        INFO     checking asrun preferences...
        Do you want to automatically configure and build code_aster (y/n)? n
        INFO     Instructions to build code_aster:
        cd $HOME/data/dev/codeaster/src
        ./waf configure
        ./waf install


Development workflow
====================

.. note:: Despite this documentation, ``rebase`` should be prefered.


See the example below (time is botton-up):

.. code-block:: none

    @    9:1ef6a397a16f main: merge 'mc-fix-cmd'
    |\
    | o  8:fbb270649f54 mc-fix-cmd: [#45678] Other fixes.
    | |
    | o  7:424270e8c082 mc-fix-cmd: [#45678] Add changes for the same feature
    | |
    | o  6:2fedc3797079 mc-fix-cmd: update to branch 'main'
    |/|
    o |  5:06487835c81e main: merge 'mc-fix-cmd'
    |\|
    | o  4:832106d39564 mc-fix-cmd: [#45678] Reopen same developer branch for new feature
    |/
    o    3:08f9c539e56e main: merge 'mc-fix-cmd'
    |\
    | o  2:cc027a3c53fb mc-fix-cmd: [#12345] Continue development for feature 12345
    | |
    | o  1:960aca07afbb mc-fix-cmd: [#12345] Start feature for 12345
    |/
    o  0:53c1e379e218 main: Main branch: main

.. Commands to create this sample tree
.. hg init
.. echo 1 > hello
.. hg add
.. hg branch main
.. hg ci -m "Main branch: main"
.. hg branch mc-fix-cmd
.. echo 1 >> hello
.. hg ci -m '[#12345] Start feature for 12345'
.. echo 1 >> hello
.. hg ci -m '[#12345] Continue development for feature 12345'
.. hg update main
.. hg merge mc-fix-cmd
.. hg ci -m "merge 'mc-fix-cmd'"
.. hg branch -f mc-fix-cmd
.. echo 1 >> hello
.. hg ci -m "[#45678] Reopen same developer branch for new feature"
.. hg update main
.. hg merge mc-fix-cmd
.. hg ci -m "merge 'mc-fix-cmd'"
.. hg up mc-fix-cmd
.. hg merge main
.. hg ci -m "update to branch 'main'"
.. echo 1 >> hello
.. hg ci -m "[#45678] Add changes for the same feature"
.. echo 1 >> hello
.. hg ci -m "[#45678] Other fixes."
.. hg up main
.. hg merge mc-fix-cmd
.. hg ci -m "merge 'mc-fix-cmd'"
.. hg log -G --template="{rev}:{node|short} {branch}: {desc|firstline}\n"

1. Start point should always be the ``main`` branch.

#. Start branch (``hg branch mc-fix-cmd``) and hack code.

#. Continue hacking and submit your work (``hg submit``).

#. **If the checkings pass, the robot automatically merges in the main branch.**

#. Reopen the branch from ``main`` for a new feature
   (``hg update main && hg branch -f mc-fix-cmd``).
   Code and submit changes (``hg submit``).

#. **If the checkings pass, the robot automatically merges in the main branch.**

#. Update working branch with last changes from ``main``
   (``hg update mc-fix-cmd && hg merge main && hg commit -m "update to branch 'main'"``).

#. Additional developments are required, continue from the same branch.

#. More changes and submission (``hg submit``).

#. **If the checkings pass, the robot automatically merges in the main branch.**


Memo:

- Reuse the same branch name.

- Continue on the same branch if you continue on the same feature.

- Reopen the branch from ``main`` for new feature (``hg branch -f ...``).


List of checkings
=================

To be accepted, the developments must pass the following checkings.

- Merge with main branch (``main``) should be trivial
  (checked by ``check_automerge.sh``).

  *In case of conflicts you have to merge the main branch with yours first.
  If your branch has several heads you have to merge them first.*

- Check sequential and parallel builds.

- Check build of the embedded documentation (checked by ``check_docs.sh``).

- Check that sequential testcases are passed (``submit`` testlist).

- Check that parallel testcases are passed (``submit`` testlist).


``hg submit`` checks the same steps except the parallel build and the parallel
testcases.

Source files are checked by *aslint*.
Issues must be validated and changed documents must be submitted.
