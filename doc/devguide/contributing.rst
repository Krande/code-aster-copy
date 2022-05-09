.. _contributing:

############
Contributing
############


*In the todo list...*

- Add a testcase for all new features.

- Check the testcases: see :ref:`testing`.

- Coding standards: http://llvm.org/docs/CodingStandards.html

  The source code may (*should*) be formatted using
  :file:`$HOME/dev/codeaster/devtools/bin/beautify` (it uses the LLVM style,
  and the parameters from :file:`.clang-format`).

- Document all new objects and methods.

  Use :file:`./check_docs.sh` script to check that new objects have been
  included in the documentation.

- Check that the documentation can be built without warnings/errors.

- The ``submit`` command from ``devtools`` checks that the fastest testcases
  pass without error (``submit`` testlist).
  The continous integration procedure checks more testcases
  (``verification`` testlist).
