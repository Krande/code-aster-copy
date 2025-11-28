#!/usr/bin/env python3
# coding: utf-8

"""
usage: python3 show_keywords.py

Show a list of the code_aster keywords.
It needs the environment to import the catalog from the code_aster installation.
"""

import os
import sys


class ExtractKeywords:
    """Loop on catalog to extract the available keywords."""

    def __init__(self):
        try:
            from code_aster.Cata.Commands import commandStore
            import code_aster.Cata.Language.SyntaxObjects as SO
        except ImportError:
            print(
                "ERROR: Environment of the developement version of code_aster "
                "not found.\n"
                "Can not import Commands descriptions."
            )
            raise

        self.commandStore = commandStore
        self.Command = SO.Command
        self.Bloc = SO.Bloc
        self.FactorKeyword = SO.FactorKeyword
        self.SimpleKeyword = SO.SimpleKeyword

    @staticmethod
    def impr_simp(comm, fact, simp, mcle):
        lines = []
        # on cherche les enumerations :
        into = mcle.definition.get("into")
        if into:
            if type(into) is str:
                lines.append("%s/%s/%s %s" % (comm, fact, simp, into))
            else:
                for valinto in into:
                    lines.append("%s/%s/%s %s" % (comm, fact, simp, valinto))
        else:
            lines.append("%s/%s/%s" % (comm, fact, simp))
        return lines


    def impr_fact(self, comm, fact, mfac):
        lines = []
        for noent in mfac.entities.keys():
            if isinstance(mfac.entities[noent], self.SimpleKeyword):
                lines.extend(self.impr_simp(comm, fact, noent, mfac.entities[noent]))
            elif isinstance(mfac.entities[noent], self.Bloc):
                lines.extend(self.impr_bloc(comm, fact, mfac.entities[noent]))
            else:
                # Certaines commandes utilisent des mots-clés facteurs emboités.
                # On imprime un message sans arrêter
                sys.stderr.write("Erreur: commande {0}/{1}/{2}\n"
                                .format(comm, fact, noent))
        return lines


    def impr_bloc(self, comm, fact, bloc):
        lines = []
        for noent in bloc.entities.keys():
            if isinstance(bloc.entities[noent], self.SimpleKeyword):
                lines.extend(self.impr_simp(comm, fact, noent, bloc.entities[noent]))
            elif isinstance(bloc.entities[noent], self.FactorKeyword):
                lines.extend(self.impr_fact(comm, noent, bloc.entities[noent]))
            elif isinstance(bloc.entities[noent], self.Bloc):
                lines.extend(self.impr_bloc(comm, fact, bloc.entities[noent]))
            else:
                raise TypeError("Erreur 2", comm, fact, noent)
        return lines

    def get(self):
        """Get all keywords."""
        lines = []
        for comm in self.commandStore.values():
            if not isinstance(comm, self.Command):
                continue
            for noent, keyword in comm.entities.items():
                if isinstance(keyword, self.SimpleKeyword):
                    lines.extend(self.impr_simp(comm.name, "--", noent, keyword))
                elif isinstance(keyword, self.FactorKeyword):
                    lines.extend(self.impr_fact(comm.name, noent, keyword))
                elif isinstance(keyword, self.Bloc):
                    lines.extend(self.impr_bloc(comm.name, "--", keyword))
                else:
                    raise TypeError("Erreur 3", comm, noent)
        return sorted(lines)


if __name__ == "__main__":
    lines = ExtractKeywords().get()
    print(os.linesep.join(lines))
