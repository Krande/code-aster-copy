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

import string
import sys
import aster
from ..Cata.Syntax import _F
from ..Commands import CALC_FERRAILLAGE, CREA_CHAMP, CREA_RESU
from ..Messages import UTMESS


def combinaison_ferraillage_ops(self, **args):
    """Command to combine results to estimate reinforcement of the structure."""
    
    resu         = args.get('RESULTAT')
    combinaison  = args.get('COMBINAISON')
    affe         = args.get('AFFE')
    codification = args.get('CODIFICATION')

    # Retriving from RESULTAT
    modele = resu.getModel()
    maillage = resu.getMesh()
    caraelem = resu.getElementaryCharacteristics()
    
    # Counting numbers of load combinations.
    nmb_cas = countCase(combinaison)
    UTMESS('I', 'COMBFERR_10', vali=nmb_cas)

    # Build instant list (lst_inst_value of float) and instant index
    # and control about existing or duplicates NUM_CAS and NUME_ORDRE.
    #
    lst_inst_index, lst_inst_value, resu_nom_cas, resu_num_ord, lst_type_combo = \
        lstInst(nmb_cas,combinaison,resu)

    UTMESS('I', 'COMBFERR_11', valk=(
        '-'.join([str(i) for i in lst_inst_value]),
        '*'.join([str(i) for i in resu_nom_cas]),
        '*'.join([str(i) for i in resu_num_ord]),
    ))

    # Controlling overwriting results
    if 'COMB_DIME_ORDRE' in resu_nom_cas:
        UTMESS('A', 'COMBFERR_9', valk='COMB_DIME_ORDRE')

    if 'COMB_DIME_ACIER' in resu_nom_cas:
        UTMESS('A', 'COMBFERR_9', valk='COMB_DIME_ACIER')

    UTMESS('I', 'COMBFERR_12', valk='\n    '.join(lst_type_combo))

    resu = algo_ferr(resu, affe, lst_inst_index, codification, lst_type_combo)
    
    # - Build result type EVOL_ELAS from MULTI_ELAS and combo type list in order
    #   to select the right verify. This because CREA_CHAMP doesn't EXTR the
    #   MAXI from multi_elas
        
    resferr = \
       evolElasFromMulti(nmb_cas,combinaison,lst_inst_value,resu)

    # Maximum reinforcement field (elementwise, component by component)
    maxiferr = CREA_CHAMP (
          RESULTAT = resferr,
          NOM_CHAM = 'FERRAILLAGE',
          TYPE_CHAM = 'ELEM_FER2_R',
          OPERATION = 'EXTR',
          TYPE_MAXI = 'MAXI',
          TYPE_RESU = 'VALE' ,
          TOUT_ORDRE = 'OUI',)

    # Instant for which maximum reinforcement field is retrieved (elementwise,
    # component by component)
    instferr = CREA_CHAMP (
          RESULTAT = resferr,
          NOM_CHAM = 'FERRAILLAGE',
          TYPE_CHAM = 'ELEM_FER2_R',
          OPERATION = 'EXTR',
          TYPE_MAXI = 'MAXI',
          TYPE_RESU = 'INST',
          TOUT_ORDRE = 'OUI',)
    
    resu = CREA_RESU(reuse = resu,
                         RESULTAT = resu,
                         OPERATION = 'AFFE',
                         TYPE_RESU = 'MULT_ELAS' ,
                         NOM_CHAM = 'FERRAILLAGE',
                         AFFE = (
                              _F(NOM_CAS = 'COMB_DIME_ACIER',
                                 CHAM_GD = maxiferr,
                                 MODELE = modele,
                                 ),),)

    resu = CREA_RESU(reuse = resu,
                         RESULTAT = resu,
                         OPERATION = 'AFFE',
                         TYPE_RESU = 'MULT_ELAS' ,
                         NOM_CHAM = 'FERRAILLAGE',
                         AFFE = (
                              _F(NOM_CAS = 'COMB_DIME_ORDRE',
                                 CHAM_GD = instferr,
                                 MODELE = modele,
                                 ),),)
                              
    nc = resu.LIST_VARI_ACCES()['NOM_CAS']
    UTMESS('I', 'COMBFERR_13', valk='\n    '.join(nc))

    return resu

def algo_ferr (resferr, affe, lst_nume_ordre, code, type_combo):
    
    #   From physical_quantities.py :
    #   FER2_R   = PhysicalQuantity(type='R',
    #   components=(
    #          'DNSXI'    : Ferraillage longitudinal inférieur suivant X (2D)
    #          'DNSXS'    : Ferraillage longitudinal supérieur suivant X (2D)
    #          'DNSYI'    : Ferraillage longitudinal inférieur suivant Y (2D)
    #          'DNSYS'    : Ferraillage longitudinal supérieur suivant Y (2D)
    #          'DNSXT'    : Ferraillage transversal suivant X            (2D)
    #          'DNSYT'    : Ferraillage transversal suivant Y            (2D) 
    #          'AYI'      : Ferraillage longitudinal inférieur suivant Y (1D)
    #          'AYS'      : Ferraillage longitudinal supérieur suivant Y (1D)
    #          'AZI'      : Ferraillage longitudinal inférieur suivant Z (1D)
    #          'AZS'      : Ferraillage longitudinal supérieur suivant Z (1D)
    #          'AST'      : Ferraillage transversal                      (1D)
    #          'ATOT'     : Ferraillage longitudinal total               (1D)
    #          'DNSVOL'   : Densité volumique de ferraillage 
    #          'CONSTRUC' : Indice de constructibilité
    #),)
            
    for idx, nume_ordre in enumerate(lst_nume_ordre):

        dic_type_comb = {}
        if type_combo [idx]  == 'ELS_CARACTERISTIQUE':
           dic_type_comb['TYPE_COMB'] = 'ELS'
        elif type_combo [idx] == 'ELS_QUASIPERMANENT':
           dic_type_comb['TYPE_COMB'] = 'ELS_QP'        
        else:
           dic_type_comb['TYPE_COMB'] = 'ELU'
           
        dic_type_comb['CODIFICATION'] = code
        lst_tmp_affe = []
                
        for i_affe in affe :
            
            dict_i_affe = i_affe.List_F()[0]
            i_affe_for_cf = dict_i_affe.copy();

            if type_combo [idx]  == 'ELU_FONDAMENTAL':
               i_affe_for_cf.update({'GAMMA_S':i_affe_for_cf['GAMMA_S_FOND']})
               i_affe_for_cf.update({'GAMMA_C':i_affe_for_cf['GAMMA_C_FOND']})
            elif type_combo [idx]  == 'ELU_ACCIDENTEL':
               i_affe_for_cf.update({'GAMMA_S':i_affe_for_cf['GAMMA_S_ACCI']})
               i_affe_for_cf.update({'GAMMA_C':i_affe_for_cf['GAMMA_C_ACCI']})

            # adjusting affe for calc_ferraillage
            i_affe_for_cf.pop('GAMMA_C_FOND')
            i_affe_for_cf.pop('GAMMA_S_FOND')
            i_affe_for_cf.pop('GAMMA_C_ACCI')
            i_affe_for_cf.pop('GAMMA_S_ACCI')

            lst_tmp_affe.append(_F(**i_affe_for_cf),)

        dic_type_comb['AFFE']=tuple(lst_tmp_affe)

        resferr = CALC_FERRAILLAGE (reuse = resferr,
                                    RESULTAT = resferr,
                                    NUME_ORDRE = nume_ordre,
                                    **dic_type_comb,)

    return resferr

# Counting numbers of load combinations.
def countCase(comb):
    nmb_cas = 0
    for idx_i_combo, i_combo in enumerate(comb) :
        # Case number from MULTI_ELAS
        lst_nomcas = i_combo.get('NOM_CAS')
        lst_numord = i_combo.get('NUME_ORDRE')
        if lst_nomcas is None:
            nmb_cas = nmb_cas + len(lst_numord)      # combinations number
        else:
            nmb_cas = nmb_cas + len(lst_nomcas)
    return nmb_cas


def lstInst(ncas, comb, resultat):

    lst_inst_value = [None] * ncas # list of float value as time for AFFE
    lst_inst_index = [None] * ncas # list of int value as index
    type_combo     = [None] * ncas # list of string as type of combo

    # Recuperer les numéros d'ordre et les noms des cas de chargement associés
    resu_nume_ordre = resultat.LIST_VARI_ACCES()['NUME_ORDRE']
    resu_nom_cas = resultat.LIST_VARI_ACCES()['NOM_CAS']

    # Elimination du vide dans [resu_nom_cas]
    for idx, val in enumerate(resu_nom_cas):
        resu_nom_cas[idx]=val.strip()

    idx_shift = 0
    for idx_i_combo, i_combo in enumerate( comb ) :
        # Case number from MULTI_ELAS
        lst_nomcas = i_combo.get('NOM_CAS')
        lst_numord = i_combo.get('NUME_ORDRE')
        if lst_nomcas is None:
            lst_combo = lst_numord
            key_name_combo = 'NUME_ORDRE'
        else:
            lst_combo = lst_nomcas
            key_name_combo = 'NOM_CAS'

        for idx_combo, val_combo in enumerate( lst_combo ):

            # type combo list couple with instant
            type_combo [idx_shift] = i_combo.get('TYPE')

            if key_name_combo == 'NUME_ORDRE':
                inst = val_combo
                if not (inst in resu_nume_ordre):
                    #le cas de charge renseigné dans COMBINAISON_FERRAILLAGE n'a pas un équivalent dans RESULTAT (issu de MACRO_ELAS_MULT)
                    UTMESS('F', 'COMBFERR_7')
                if inst in lst_inst_index:
                    #le cas de charge renseigné dans COMBINAISON_FERRAILLAGE a déjà été renseigné
                    UTMESS('F', 'COMBFERR_4')
            else:
                if not (val_combo in resu_nom_cas):
                    UTMESS('F', 'COMBFERR_6')
                # Index needs to be [inst_index + 1] because MUME_ORDRE is index + 1
                inst = resu_nom_cas.index(val_combo) + 1
                if inst in lst_inst_index:
                    UTMESS('F', 'COMBFERR_5')

            lst_inst_index [ idx_shift ] = inst
            lst_inst_value [ idx_shift ] = float(inst)
            idx_shift = idx_shift + 1

    return lst_inst_index, lst_inst_value, resu_nom_cas, resu_nume_ordre, type_combo

# Build result type EVOL_ELAS from MULTI_ELAS
def evolElasFromMulti(ncas, comb, lst_inst_value, resu):

    modele = resu.getModel()
    caraelem = resu.getElementaryCharacteristics()
    
    __EFGE         = [None] * ncas
    lst_AFFE_EFGE  = [None] * ncas

    idx_shift = 0
    for idx_i_combo, i_combo in enumerate( comb ) :

        # Case number from MULTI_ELAS
        lst_nomcas = i_combo.get('NOM_CAS')
        lst_numord = i_combo.get('NUME_ORDRE')
        if lst_nomcas is None:
            lst_combo = lst_numord    # list with combinations
            key_name_combo = 'NUME_ORDRE'
        else:
            lst_combo = lst_nomcas
            key_name_combo = 'NOM_CAS'


        for idx_combo, val_combo in enumerate( lst_combo ):

            dic_idx_combo = { key_name_combo : val_combo }

            __EFGE[ idx_shift ] = CREA_CHAMP (
            OPERATION = 'EXTR',
            RESULTAT = resu,
            TYPE_CHAM = 'ELEM_FER2_R',
            NOM_CHAM = 'FERRAILLAGE',
            **dic_idx_combo
            )
            
            lst_AFFE_EFGE [ idx_shift ] = _F (
                 INST = lst_inst_value [ idx_shift ],
              CHAM_GD = __EFGE [ idx_shift ],
               MODELE = modele,
            CARA_ELEM = caraelem,
            )

            idx_shift = idx_shift  + 1

    resferr = CREA_RESU (
        OPERATION='AFFE',
        TYPE_RESU ='EVOL_ELAS',
        NOM_CHAM ='FERRAILLAGE',
        AFFE = lst_AFFE_EFGE,
    )

    return resferr
