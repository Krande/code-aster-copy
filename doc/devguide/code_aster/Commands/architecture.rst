######################
Structure of a command
######################

.. automodule:: code_aster.CodeCommands
   :show-inheritance:
   :members:
   :special-members: __init__


*******************
Catalog description
*******************

Each command must define its user syntax through its catalog.
The checking of the user keywords is performed by the Catalog objects
(see :py:mod:`~code_aster.Cata.Language.SyntaxChecker`).

See :py:mod:`~code_aster.Cata.Language.SyntaxObjects` for the type of the values
of the user's keywords after the syntax checkings.
For example, the value is always a *list* if ``max > 1``. It is not necessary to
check again in the implementation of the
:py:meth:`~code_aster.Supervis.ExecuteCommand.ExecuteCommand.exec_` function of the command.


*******************************
Utility functions for executors
*******************************

.. automodule:: code_aster.CodeCommands.operator
   :show-inheritance:
   :members:
   :special-members: __init__
