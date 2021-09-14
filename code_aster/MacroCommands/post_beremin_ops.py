# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

import sys    
    
"""
WARNING: (TEMPORARY) FILE TO PUT IN
                     ~/dev/codeaster/src/code_aster/MacroCommands/

Compute the Beremin probability of failure

- Aboubakr Amzil filters are not yet implemented
- Optimization is probably possible if we loop through the
  times (could be useful for Charpy tests for example)

 - Strain model -> Stress to consider (6 components: XX, YY, ZZ, XY, XZ, YZ)
   - "PETIT"      -> "SIEF"
   - "PETIT_REAC" -> "SIEF"
   - "GDEF_LOG"   -> "VARI_ELGA_NOMME"
   - "SIMO_MIEHE" -> To retrieve Kirchhof stresses ? Not obvious ...
                     Normally, we can do without
                     Use of HHO models instead ?

 - Integral must be computed on the initial configuration

_biblio : CR-T66-2017-132, Lorentz
"""
import os
import pathlib
import math
import medcoupling as mc
import numpy as np
import aster

from ..Cata.Syntax import _F
from ..Commands import (DEFI_FICHIER, IMPR_RESU,CREA_TABLE)
from ..Helpers import LogicalUnitFile, ReservedUnitUsed

from libaster import EntityType
from ..Messages import UTMESS

# projection with medcoupling (set _FB = True) or with code_aster
# (set _FB = False)

_FB = False

def post_beremin_ops(self, RESULTAT, GROUP_MA,
                     DEFORMATION, FILTRE_SIGM, 
                     COEF_MULT=None
                     ,UNITE=None
                    ):
    """
    Beremin post-processor

    Inputs :
    ========
    - RESULTAT: aster concept of mechanical result. Must be the same as
                filename_in
    - MED_FILENAME_IN: MED filename of code_aster mechanics result
    - GROUP_MA: mesh cells group on which Beremin post-treatment is carried out
                Only 1 mesh cells group authorized
    - DEFORMATION: "PETIT", "PETIT_REAC", "GDEF_LOG" or "SIMO_MIEHE"
    - FILTRE_SIGM: "SIGM_ELGA" or "SIGM_ELMOY".
                 SIGM: FIELD WITH 4 OR 6 STRESSES COMPONENTS
    - COEF_MULT: as explained in doc aster u4.81.22 (weibull) p.15

    Outputs :
    =========
    - TABLE : Through IMPR_TABLE is available 
                         information (Weibull stress and probability of failure)
    - MED file : MED filename containing maximum of major principal
                         stress by where the plasticity is active (possible if UNITE=number is given in command file)
                         Result of type ELGA (apply ElgaFieldToSurface filter
                         in Paraview to see it)
    """
    
    
    nomfich = LogicalUnitFile.filename_from_unit(UNITE)
    
    
    pathlib.Path(os.path.dirname(nomfich)).mkdir(parents=True,  exist_ok=True)
    
    
    tab_out = [["nume_ordre","inst","lieu","entite",
              "sigma_weibull","proba_weibull","sigma_weibull**m"
              ]]
    
    MED_FILENAME_IN = os.path.join(os.getcwd(), "tmp.med")
    DEFI_FICHIER(UNITE=46,
                 ACTION="ASSOCIER",
                 TYPE="BINARY",
                 FICHIER=MED_FILENAME_IN)

    IMPR_RESU(UNITE=46,
              RESU=_F(RESULTAT=RESULTAT))

    DEFI_FICHIER(ACTION="LIBERER", UNITE=46)


    if len(GROUP_MA)>1:UTMESS('F', 'RUPTURE1_87')
    grma = GROUP_MA[0]
    data_in = get_data(RESULTAT, MED_FILENAME_IN, DEFORMATION, grma)

    if [elt[2] for elt in mc.GetAllFieldIterations(
                MED_FILENAME_IN, data_in[
                    "stress_fieldname"])] != RESULTAT.LIST_VARI_ACCES()["INST"]:
            UTMESS('F', 'RUPTURE1_80')
    for disct in mc.GetAllFieldIterations(MED_FILENAME_IN,
                                              data_in["stress_fieldname"]):

            stress = mc.ReadFieldGauss(MED_FILENAME_IN,
                                       mc.GetMeshNamesOnField(
                                           MED_FILENAME_IN,
                                           data_in["stress_fieldname"])[0], 0,
                                       data_in["stress_fieldname"],
                                       disct[0], disct[1])[data_in["grpma"]]

            sigw = stress.deepCopy()
            sigw.setName("SIGW")
            if _FB:
                data_in["MED_FILENAME_IN"] = MED_FILENAME_IN
            sigw.setArray(weibull_stress(
                grma, DEFORMATION, FILTRE_SIGM, data_in, stress, array_cmp(
                    mc.ReadFieldGauss(
                        MED_FILENAME_IN,
                        mc.GetMeshNamesOnField(
                            MED_FILENAME_IN,
                            data_in["stress_fieldname"])[0], 0,
                        data_in["plasti_fieldname"], disct[0],
                        disct[1])[data_in["grpma"]],
                    ("EPSPEQ", "INDIPLAS")), disct))
            sigw.checkConsistencyLight()

            if data_in["first"]:
                mc.WriteField(nomfich, sigw, True)
                data_in["first"] = False
            else:
                mc.WriteFieldUsingAlreadyWrittenMesh(nomfich, sigw)
            
            _1 = sigweibull_proba(grma, COEF_MULT, data_in, sigw)
            tab_out.append((
                disct[1], disct[2], grma, "GROUP_MA",
                _1[0], _1[1], _1[0]**data_in["wb_kwd"][grma]["M"]
            ))
            del _1

    #$
    class csv:
        def __init__(self, t_out):
            self.tabout = [[_, None] for _ in t_out[0]]
            for _ in range(len(self.tabout)):
                self.tabout[_][1] =  [i_[_] for i_ in t_out[1:]]
        def getname(self, col):
            return self.tabout[col][0]
        def getvals(self, col):
            return self.tabout[col][1]
    csvO = csv(tab_out)
    
    tabout = CREA_TABLE(LISTE=(
        _F(PARA=csvO.getname(0), LISTE_I=csvO.getvals(0)),
        _F(PARA=csvO.getname(1), LISTE_R=csvO.getvals(1)),
        _F(PARA=csvO.getname(2), LISTE_K=csvO.getvals(2)),
        _F(PARA=csvO.getname(3), LISTE_K=csvO.getvals(3)),
        _F(PARA=csvO.getname(4), LISTE_R=csvO.getvals(4)),
        _F(PARA=csvO.getname(5), LISTE_R=csvO.getvals(5)),
        _F(PARA=csvO.getname(6), LISTE_R=csvO.getvals(6)),
    ))

    return tabout

def weibull_stress(GROUP_MA, DEFORMATION, FILTRE_SIGM, data_in, stress, plas, disct):
    """
    Inputs :
    - GROUP_MA: mesh cells group on which Beremin post-treatment is carried out
                Only 1 mesh cells group authorized
    - DEFORMATION:  "PETIT", "PETIT_REAC", "GDEF_LOG"
    - FILTRE_SIGM: "SIGM_ELGA" or "SIGM_ELMOY".
    - data_in: dictionary of useful data
    - stress: stress field at instant disct[2]
    - plas: cumulated plastic strain field and plasticity indicator
    - disct: triplet (iteration, sub-iteration, physical instant)
    Input and output:
    - data_in: dictionary of useful data

    """
    sig = array_cmp(stress, data_in["stress_cmpnames"])
    if (data_in["dim_geom"] == 2) & (DEFORMATION in ("PETIT", "PETIT_REAC")):
        zero = mc.DataArrayDouble(len(sig))
        zero.fillWithValue(0)
        sig.meldWith(zero)
        sig.meldWith(zero)

    if data_in["first"] & (FILTRE_SIGM == "SIGM_ELMOY"):
        for idcell in range(data_in["nbofcells"]):
            data_in["l_idcell"].append(
                stress.computeTupleIdsToSelectFromCellIds(idcell))
    # sigma_u: cleavage stress. 2 cases depending on whether it depends
    #          on the temperature or not
    #   - real
    #   - nom du fichier med où est rangée la contrainte de clivage évaluée
    #     sur le maillage mécanique, en fonction de la température
    #     Le fichier MED doit contenir un champ SIEF_ELGA avec la contrainte
    #     de clivage sur le groupe "group_ma", rangée dans la composante SIX#X
    #     La discrétisation temporelle (ainsi que le maillage) de ce
    #     fichier et de MED_FILENAME_IN doivent donc être identiques

    # The following keyword is necessary only if sigma_u is a function of
    # temperature (hence a MED file here)
    # - "sigm_cnv": conventional cleavage stress
    if not _FB:
        if "sigm_cnv" in data_in and isinstance(data_in["wb_kwd"][GROUP_MA][
                "SIGM_REFE"],
                                                str):
            sig = sigma_u_t(GROUP_MA, data_in["grpma"], sig, data_in, disct[2])
    else:
        if "sigm_cnv" in data_in and "temperature" in data_in \
           and "sigma_u" in data_in:
            fb_sigma_u_t(GROUP_MA, stress, data_in, disct[2])


    if FILTRE_SIGM == "SIGM_ELGA":
        smax = maxsig1_elga(sig)
    elif FILTRE_SIGM == "SIGM_ELMOY":
        smax = maxsig1_elmoy(sig, data_in)

    # Points to be taken into account regarding plasticity state
    # and evolution
    zone = mc.DataArrayDouble(len(sig))
    zone.fillWithValue(1)
    # kappa_c: threshold under which we consider the equivalent cumulated
    #          plastic strain as 0
    zone[plas[:, 0].findIdsLowerThan(data_in["wb_kwd"][GROUP_MA][
        "SEUIL_EPSP_CUMU"])] = 0
    zone[plas[:, 1].findIdsLowerThan(0.9)] = 0
    seff = smax*zone

    if data_in["wb_stress"] is not None:
        data_in["wb_stress"].meldWith(seff)
        data_in["wb_stress"] = data_in["wb_stress"].maxPerTuple()
    else:
        data_in["wb_stress"] = seff

    data_in["wb_stress"].setInfoOnComponent(0, "W")
    return data_in["wb_stress"]

def maxsig1_elga(sig):
    """
    Calcul du maximum de la contrainte principale majeure pour tous les points
    de Gauss (option = "sigm_elga")

    input :
    - sig : tableau des contraintes

    output :
    - smax : maximum des valeurs propres pour tous les points de Gauss
    """
    sig_np = sig.toNumPyArray()
    smax = mc.DataArrayDouble(len(sig))
    
    for numline in range(len(sig)):
        sig_ptga = sig_np[numline, :]
        
        d = mc.DataArrayDouble(1, 6)
        for _1 in range(6):
            d[0,_1] = sig_ptga[_1]
        smax[numline] = d.eigenValues().maxPerTuple()[0]
    
    return smax

def maxsig1_elmoy(sig, data_in):
    """
    Calcul du maximum de la contrainte principale majeure pour tous les points
    de Gauss (option = "sigm_elmoy")

    input :
    - sig : tableau des contraintes
    _data_in : dictionnaire de données

    output :
    - smax : maximum des valeurs propres pour tous les points de Gauss
    """
    sig_np = sig.toNumPyArray()
    smax = mc.DataArrayDouble(len(sig))
    meansig = mc.DataArrayDouble(data_in["nbofcells"],
                                 sig.getNumberOfComponents())
    s1max_np = []
    for idcell in range(data_in["nbofcells"]):
        for idcmp in range(sig.getNumberOfComponents()):
            meancell = 0.
            for idptga in data_in["l_idcell"][idcell]:
                meancell += sig_np[idptga[0]][idcmp]
            meansig[idcell, idcmp] = meancell/len(
                data_in["l_idcell"][idcell])

        d = mc.DataArrayDouble(1, 6)##
        for _1 in range(6):
            d[0,_1] = meansig[idcell, _1]
        s1max_np.append(d.eigenValues().maxPerTuple()[0])
    
    for idcell in range(data_in["nbofcells"]):
        for idptga in data_in["l_idcell"][idcell]:
            smax[idptga] = s1max_np[idcell]

    return smax

def sigweibull_proba(GROUP_MA, COEF_MULT, data_in, sigw):
    """
    Calcul de la contrainte de rupture et de la probabilité de rupture
    Inputs :
    - COEF_MULT
    - data_in
    - sigw
    Outputs :
    - sig_weibull
    - proba
    """
    density = sigw.deepCopy()
    p_aux = density.getArray()**data_in["wb_kwd"][GROUP_MA]["M"]/data_in[
        "wb_kwd"][GROUP_MA]["VOLU_REFE"]

    if data_in["axis"] == "OUI":
        rayon_axis = sigw.getLocalizationOfDiscr().getValuesAsTuple()
        for ptga in range(p_aux.getNbOfElems()):
            p_aux[ptga] = p_aux[ptga]*rayon_axis[ptga][0]
    density.setArray(p_aux)
    integ = density.integral(0, True)

    if "sigm_cnv" in data_in:
        proba = 1 - math.exp(-integ*COEF_MULT/data_in["sigm_cnv"]**data_in[
            "wb_kwd"][GROUP_MA]["M"])
    else:
        proba = 1 - math.exp(-integ*COEF_MULT/data_in[
            "wb_kwd"][GROUP_MA]["SIGM_REFE"]**data_in["wb_kwd"][GROUP_MA]["M"])
    return (integ*COEF_MULT)**(1/data_in["wb_kwd"][GROUP_MA]["M"]), proba

def sigma_u_t(GROUP_MA, grpma, sig, data_in, _tfm):
    """
    Calcul de sig si sigma_u est dépendant de la température
    Inputs :
    - grpma :
    - sig :
    - data_in :
    - _tfm :
    Output :
    - sig :
    """
    compteur = 0
    for field in mc.MEDFileFields(data_in["wb_kwd"][GROUP_MA]["SIGM_REFE"],
                                  False).getFieldsNames():
        if field[8:] == "SIEF_ELGA":
            stress_fieldname = field
            compteur += 1
    if compteur != 1:
        raise PostBereminError("SIEF_ELGA non trouvé 1 fois exactement")

    if _tfm in [disct[2] for disct in mc.GetAllFieldIterations(
            data_in["wb_kwd"][GROUP_MA]["SIGM_REFE"], stress_fieldname)]:
        for (_dtt, _itt, _tft) in mc.GetAllFieldIterations(data_in["wb_kwd"][
                GROUP_MA]["SIGM_REFE"],
                                                           stress_fieldname):
            if _tft == _tfm:
                sigma_u = array_cmp(mc.ReadFieldGauss(
                    data_in["wb_kwd"][GROUP_MA]["SIGM_REFE"],
                    mc.GetMeshNamesOnField(
                        data_in["wb_kwd"][GROUP_MA]["SIGM_REFE"],
                        stress_fieldname)[0], 0,
                    stress_fieldname, _dtt, _itt)[grpma], ("SIXX",))

                for ptga in range(sig.getNumberOfTuples()):
                    if sigma_u[ptga] != 0.:
                        sig[ptga] = sig[ptga]*data_in["sigm_cnv"]/sigma_u[ptga]
                    else:
                        raise PostBereminError(
                            "sigma_u nulle au point de Gauss {}".format(ptga))

    else:
        raise PostBereminError("Discrétisations temporelles de sigma_u et "
                               "du résultat de contraintes différentes")

    return sig

def fb_sigma_u_t(GROUP_MA, stress, data_in, _tfm):
    """
    Cleavage stress depending on temperature
    """
    #inefficace pour l'instant
    if "sigm_cnv" in data_in and not isinstance(data_in["wb_kwd"][GROUP_MA][
            "SIGM_REFE"], float):
        for (_dtt, _itt, _tft) in mc.GetAllFieldIterations(
                data_in["temperature"], mc.MEDFileFields(
                    data_in["temperature"], False).getFieldsNames()[0]):
            if _tft == _tfm:
                temp = mc.ReadFieldNode(
                    data_in["temperature"],
                    mc.MEDFileMeshes(data_in["temperature"]).getMeshesNames()[
                        0], 0, mc.MEDFileFields(data_in["temperature"],
                                                False).getFieldsNames()[0],
                    _dtt, _itt)

                mesh = temp.getMesh()
                _mesht3 = mesh.deepCopy()
                _mesht3.convertQuadraticCellsToLinear()
                if _mesht3.getMeshDimension() == 2:
                    _mesht3.simplexize(0)
                elif _mesht3.getMeshDimension() == 3:
                    _mesht3.simplexize(mc.PLANAR_FACE_5)
                    #_mesht3.tetrahedrize(mc.PLANAR_FACE_5)
                #_mesht3.convertLinearCellsToQuadratic()
                _mesht3.writeVTK("/home/F82953/debug.vtu")
                _mesht3.zipCoords()
                (_ok, _subset) = mesh.getCoords().areIncludedInMe(
                    _mesht3.getCoords(), 1.e-8)
                assert _ok

                _tempt3 = temp.New(mc.ON_NODES)
                _tempt3.setMesh(_mesht3)
                _tempt3.setArray(temp.getArray()[_subset])
                temp = _tempt3
                temp.setName("temperature")
                temp.writeVTK("/home/F82953/debug_temp.vtu")
                #1. Projeter le champ temp sur le maillage mécanique. Ici c'est
                #identique mais dans le cas général, les maillages thermique et
                #mécanique ne sont pas identiques. On appelle temp_mec le champ
                #obtenu
                #2. Interpoler temp_mec aux point de Gauss. Appelons le champ
                #obtenu temp_mec_elga
                #3. Appliquer la formule sigref (concept aster)
                #   applyFuncOnThis("sqrt(x)")
                # Valeur de la temperature aux points de Gauss
                xpg = stress.getLocalizationOfDiscr()
                tpg = temp.getValueOnMulti(xpg)
                print(tpg)
                #data_in["sigme_refe"](tpg)
                #for ptga in range(sig.getNbOfElems()):
                #    sig[ptga] = sig[ptga]*data_in["sigm_cnv"]/1.
                    #sig[ptga] = sig[ptga]*data_in["sigm_cnv"]/sigm_refe[
                    #ptga]
#        return sig

#class PostBereminError(Exception):    
#    """Base class for post_beremin errors."""

def get_data(RESULTAT, MED_FILENAME_IN, DEFORMATION, GROUP_MA):
    """
    Check that the file MED_FILENAME_IN contains necessary data
    for Weibull computation
    Function of deformation

    Output:
    Dictionary of data
    """

    if mc.MEDFileData(MED_FILENAME_IN).getNumberOfMeshes() != 1:
        UTMESS('F', 'RUPTURE1_81')

    if not RESULTAT.getModel().isMechanical():
        UTMESS('F', 'RUPTURE1_82')
    else:
        axis = aster.dismoi("AXIS", RESULTAT.getModel().getName(), "MODELE",
                            "F")[-1]
        dim_geom = aster.dismoi("DIM_GEOM", RESULTAT.getModel().getName(),
                                "MODELE", "F")[-2]
    wb_kwd = dict()
    _chmater = RESULTAT.getMaterialField().getVectorOfPartOfMaterialField()
    nomres = ["SIGM_REFE", "M", "VOLU_REFE", "SEUIL_EPSP_CUMU"]
    for (indi_mate, _) in enumerate(_chmater):
        (vale_wbmat, codret_wbmat) = _chmater[indi_mate].getVectorOfMaterial()[
            0].RCVALE("WEIBULL", nomres=nomres, stop=0)
        if codret_wbmat == (0, 0, 0, 0):
            meshEntity = _chmater[indi_mate].getMeshEntity()
            entityType = meshEntity.getType()
            if entityType is EntityType.GroupOfCellsType:
                for grpname in _chmater[indi_mate].getMeshEntity().getNames():
                    wb_kwd.update({
                        grpname:dict(list(zip(nomres, list(vale_wbmat))))})
            elif entityType is EntityType.AllMeshEntitiesType:
                wb_kwd.update({
                    GROUP_MA:dict(list(zip(nomres, list(vale_wbmat))))})
    compteur = {"plas":0, "stress_s":0, "stress_t":0}

    for field in mc.MEDFileFields(MED_FILENAME_IN, False).getFieldsNames():
        if field[8:] == "VARI_ELGA_NOMME":
            plasti_fieldname = field
            compteur["plas"] += 1

        if DEFORMATION in ("PETIT", "PETIT_REAC"):
            if field[8:] == "SIEF_ELGA":
                stress_fieldname = field
                if dim_geom == 3:
                    stress_cmpnames = ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ",
                                       "SIYZ")
                elif dim_geom == 2:
                    stress_cmpnames = ("SIXX", "SIYY", "SIZZ", "SIXY")
                compteur["stress_s"] += 1
        elif DEFORMATION == "GDEF_LOG":
            if field[8:] == "VARI_ELGA_NOMME":
                stress_fieldname = field
                stress_cmpnames = ("TXX", "TYY", "TZZ", "TXY", "TXZ", "TYZ")
                compteur["stress_t"] += 1
        elif DEFORMATION == "SIMO_MIEHE":
            UTMESS('F', 'RUPTURE1_83')

    if compteur["plas"] is not 1:
        UTMESS('F', 'RUPTURE1_84')
    if compteur["stress_s"] != 1 & (DEFORMATION in ("PETIT", "PETIT_REAC")):
        UTMESS('F', 'RUPTURE1_85')
    if compteur["stress_t"] != 1 & (DEFORMATION == "GDEF_LOG"):
        UTMESS('F', 'RUPTURE1_86')

    return {"axis":axis,
            "dim_geom":dim_geom,
            "first": True,
            "grpma": mc.MEDFileMesh.New(MED_FILENAME_IN).getGroupArr(
                0, GROUP_MA),
            "l_idcell":[],
            "nbofcells": mc.MEDFileMesh.New(MED_FILENAME_IN).getGroup(
                0, GROUP_MA).getNumberOfCells(),
            "plasti_fieldname": plasti_fieldname,
            "stress_fieldname": stress_fieldname,
            "stress_cmpnames": stress_cmpnames,
            "wb_stress" : None,
            "wb_kwd":wb_kwd}

def array_cmp(field, cmp_names):
    """
    Returns the array of the field reduced to a list of component Names
    cmp_names
    """
    arr = field.getArray()
    return arr.keepSelectedComponents([arr.getInfoOnComponents().index(comp)
                                       for comp in cmp_names])
