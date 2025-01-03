# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

from . import *
from .sd_carte import sd_carte
from .sd_champ import sd_champ
from .sd_compor import sd_compor
from .sd_mater import sd_mater
from .sd_util import *

# Remarque :
# ------------
# la sd_cham_mater est produite par la seule commande AFFE_MATERIAU.
# C'est pourquoi, on fera appel parfois pour la décrire au vocabulaire de cette
# commande.


class sd_cham_mater_varc(AsBase):
    # ----------------------------------
    # la sd_cham_mater_varc est la partie de la sd_cham_mater correspondant aux variables de commande
    # (mot clé AFFE_VARC)

    nomj = SDNom(fin=8)
    CVRCNOM = AsVK8()
    CVRCGD = AsVK8()
    CVRCVARC = AsVK8()
    CVRCCMP = AsVK8()

    def exists(self):
        return self.CVRCVARC.exists

    # indirection via CVRCVARC:
    def check_cham_mater_i_CVRCVARC(self, checker):
        if not self.exists():
            return
        lnom = self.CVRCVARC.get()
        for nom in lnom:
            nom2 = self.nomj()[:8] + "." + nom + ".1"
            sd2 = sd_carte(nom2)
            sd2.check(checker)

            nom2 = self.nomj()[:8] + "." + nom + ".2"
            sd2 = sd_carte(nom2)
            sd2.check(checker)

            # dans le cas de la temperature, on cree parfois des cartes sous-terraines
            # (lorsque l'on indique VALE_REF sans donner de CHAM_GD/EVOL) :
            if nom.strip() == "TEMP":
                desc = sd2.DESC.get()
                ngdmax = desc[1]
                ngdedi = desc[2]
                vale = sd2.VALE.get_stripped()
                ncmp = len(vale) // ngdmax
                assert len(vale) == ncmp * ngdmax, (ngdmax, ncmp, vale)
                for kedit in range(ngdedi):
                    assert vale[ncmp * kedit + 0] == "TEMP", (vale, kedit, ncmp)
                    if vale[ncmp * kedit + 1] == "CHAMP":
                        sd3 = sd_champ(vale[ncmp * kedit + 2])
                        sd3.check(checker)

    # vérification des objets .CVRC* :
    def check_CVRC(self, checker):
        if not self.exists():
            return
        cvrccmp = self.CVRCCMP.get()
        cvrcnom = self.CVRCNOM.get_stripped()
        cvrcgd = self.CVRCGD.get_stripped()
        cvrcvarc = self.CVRCVARC.get_stripped()

        # Les 5 objets ont la meme longueur > 0 :
        nbcvrc = len(cvrcnom)
        assert nbcvrc > 0, self
        assert len(cvrccmp) == nbcvrc, (cvrccmp, cvrcnom, self)
        assert len(cvrcgd) == nbcvrc, (cvrcgd, cvrcnom, self)
        assert len(cvrcvarc) == nbcvrc, (cvrcvarc, cvrcnom, self)

        # Les 4 objets sont "non blancs" :
        sdu_tous_non_blancs(self.CVRCCMP, checker)
        sdu_tous_non_blancs(self.CVRCNOM, checker)
        sdu_tous_non_blancs(self.CVRCGD, checker)
        sdu_tous_non_blancs(self.CVRCVARC, checker)

        # les noms des CRVC doivent etre differents:
        sdu_tous_differents(self.CVRCNOM, checker)


class sd_cham_mater(AsBase):
    # =============================
    nomj = SDNom(fin=8)

    # CHAMP_MAT est une carte contenant la liste des noms de matériaux
    # affectées sur les mailles du maillage.
    CHAMP_MAT = sd_carte()

    # si AFFE_VARC :
    varc = Facultatif(sd_cham_mater_varc(SDNom(nomj="")))

    # si AFFE_COMPOR :
    # COMPOR est une carte définissant les sd_compor affectés
    # sur les mailles du maillage
    COMPOR = Facultatif(sd_carte())

    def check_CHAMP_MAT(self, checker):
        # ----------------------------------
        # on vérifie que la carte .CHAMP_MAT contient bien des noms de
        # matériau.
        vale = self.CHAMP_MAT.VALE.get_stripped()
        desc = self.CHAMP_MAT.DESC.get()
        numgd = desc[0]
        n_gd_edit = desc[2]
        assert sdu_nom_gd(numgd) == "NOMMATER", (desc, sdu_nom_gd(numgd))
        ncmp_max = len(sdu_licmp_gd(numgd))
        assert ncmp_max == 30, ncmp_max

        for kedit in range(n_gd_edit):
            v1 = vale[kedit * ncmp_max : (kedit + 1) * ncmp_max]
            ktref = None
            for k1 in range(len(v1)):
                x1 = v1[k1]
                if x1 == "":
                    continue

                # cas particulier : 'TREF=>', suivi de la temperature de reference ecrite sur 3 k8
                if x1 == "TREF=>":
                    ktref = k1
                    continue
                if ktref:
                    assert k1 == ktref + 1, (k1, ktref)
                    # on doit retrouver la valeur de TREF :
                    x1 = v1[k1] + v1[k1 + 1] + v1[k1 + 2]
                    try:
                        tref = float(x1)
                    except:
                        assert 0, " On doit trouver la valeur de TREF: " + x1
                    break

                # cas général : x1 est un nom de sd_mater
                sd2 = sd_mater(x1)
                sd2.check(checker)

    def check_COMPOR(self, checker):
        # ----------------------------------
        # on vérifie (un peu)  la carte .COMPOR (si elle existe)
        desc = self.COMPOR.DESC.get()
        if not desc:
            return
        vale = self.COMPOR.VALE.get_stripped()
        numgd = desc[0]
        n_gd_edit = desc[2]
        assert sdu_nom_gd(numgd) == "COMPOR", (desc, sdu_nom_gd(numgd))
        ncmp_max = len(sdu_licmp_gd(numgd))

        for kedit in range(n_gd_edit):
            v1 = vale[kedit * ncmp_max : (kedit + 1) * ncmp_max]
            assert v1[3] == "COMP_INCR", v1
            sd2 = sd_compor(v1[6].split(".")[0].strip())
            sd2.check(checker)

            for x1 in v1[7:]:
                assert x1 in ("", "VIDE", "0"), v1
