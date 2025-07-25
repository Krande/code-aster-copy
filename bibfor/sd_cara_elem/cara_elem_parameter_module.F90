! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

module cara_elem_parameter_module
!
!
!   Pour l'opérateur AFFE_CARA_ELEM
!
! --------------------------------------------------------------------------------------------------
! person_in_charge: jean-luc.flejou at edf.fr
!
    implicit none
!
    integer(kind=8), parameter :: ACE_NOTHING = 0
!
! Nombre total de mots clefs facteur
    integer(kind=8), parameter :: ACE_NB_MCLEF = 16
    character(len=16), parameter :: ACE_MCLEF(ACE_NB_MCLEF) = [ &
                                    'POUTRE          ', 'COQUE           ', &
                                    'DISCRET         ', 'ORIENTATION     ', &
                                    'CABLE           ', &
                                    'BARRE           ', 'MASSIF          ', &
                                    'POUTRE_FLUI     ', 'RIGI_PARASOL    ', &
                                    'GRILLE          ', 'RIGI_MISS_3D    ', &
                                    'DISCRET_2D      ', 'MEMBRANE        ', &
                                    'MASS_AJOU       ', 'MULTIFIBRE      ', &
                                    'MASS_REP        ']
!
    integer(kind=8), parameter :: ACE_POUTRE = 1
    integer(kind=8), parameter :: ACE_COQUE = 2
    integer(kind=8), parameter :: ACE_DISCRET = 3
    integer(kind=8), parameter :: ACE_ORIENTATION = 4
    integer(kind=8), parameter :: ACE_CABLE = 5
    integer(kind=8), parameter :: ACE_BARRE = 6
    integer(kind=8), parameter :: ACE_MASSIF = 7
    integer(kind=8), parameter :: ACE_POUTRE_FLUI = 8
    integer(kind=8), parameter :: ACE_RIGI_PARASOL = 9
    integer(kind=8), parameter :: ACE_GRILLE = 10
    integer(kind=8), parameter :: ACE_RIGI_MISS_3D = 11
    integer(kind=8), parameter :: ACE_DISCRET_2D = 12
    integer(kind=8), parameter :: ACE_MEMBRANE = 13
    integer(kind=8), parameter :: ACE_MASS_AJOU = 14
    integer(kind=8), parameter :: ACE_MULTIFIBRE = 15
    integer(kind=8), parameter :: ACE_MASS_REP = 16
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: ACE_GRMA_MA(4) = [ &
                                    'GROUP_MA        ', 'MAILLE          ', &
                                    'GROUP_MA_POI1   ', 'GROUP_MA_SEG2   ']
    character(len=8), parameter :: ACE_GRMA_TY(4) = [ &
                                   'GROUP_MA', 'MAILLE  ', &
                                   'GROUP_MA', 'GROUP_MA']
    integer(kind=8), parameter :: ACE_GR_MAI = 1
    integer(kind=8), parameter :: ACE_MAILLE = 2
    integer(kind=8), parameter :: ACE_GR_PO1 = 3
    integer(kind=8), parameter :: ACE_GR_SE2 = 4
!
    integer(kind=8), parameter :: ACE_NB_GRMA_MA = 4
    integer(kind=8), parameter :: MCLEF_GRP_MA(ACE_NB_GRMA_MA*ACE_NB_MCLEF) = [ &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_GR_PO1, ACE_GR_SE2, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_PO1, ACE_GR_SE2, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_GR_PO1, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_MAILLE, ACE_NOTHING, ACE_NOTHING, &
                                  ACE_GR_MAI, ACE_GR_PO1, ACE_NOTHING, ACE_NOTHING]
!
! --------------------------------------------------------------------------------------------------
! Toutes les cartes créées par AFFE_CARA_ELEM
    integer(kind=8), parameter :: ACE_NB_CARTE = 4
    integer(kind=8), parameter :: ACE_NB_CARTE_CMP = 3
    character(len=10), parameter :: ACE_CARTE(ACE_NB_CARTE*ACE_NB_CARTE_CMP) = [ &
                                    '.CARDINFO ', 'CINFDI_R  ', 'DISCRET   ', &
                                    '.CARDISCK ', 'CADISK_R  ', 'DISCRET   ', &
                                    '.CARDISCM ', 'CADISM_R  ', 'DISCRET   ', &
                                    '.CARDISCA ', 'CADISA_R  ', 'DISCRET   ']
!
    integer(kind=8), parameter :: ACE_CAR_DINFO = 1
    integer(kind=8), parameter :: ACE_CAR_DISCK = 2
    integer(kind=8), parameter :: ACE_CAR_DISCM = 3
    integer(kind=8), parameter :: ACE_CAR_DISCA = 4

! #define ACE_NB_CARTE  14
!         '.CARCABLE', 'CACABL_R ', 'CABLE     ',
!         '.CARGENBA', 'CAGNBA_R ', 'BARRE     ',
!         '.CARMASSI', 'CAMA_R   ', 'MASSIF    ',
!         '.CARCOQUE', 'CACOQU_R ', 'COQUE     ',
!         '.CARCOQUF', 'CACOQU_F ', 'COQUE    ',
!         '.CARARCPO', 'CAARPO_R ', 'POUTRE    ',
!         '.CARGENPO', 'CAGNPO_R ', 'POUTRE    ',
!         '.CARGEOPO', 'CAGEPO_R ', 'POUTRE    ',
!         '.CARPOUFL', 'CAPOUF_R ', 'POUTREFLUI',
!         '.CVENTCXF', 'VENTCX_F ', 'VENT     ',
!         '.CARORIEN', 'CAORIE_R ', 'ORIENT   ',
!
! --------------------------------------------------------------------------------------------------
!
!       ACE_NB_ELEMENT              : Nombre de type d'élements différents
!       ACE_NU_(el)                 : Numéro du type de l'élément (POUTRE DISCRET COQUE ...)
!       ACE_NM_(el)                 : Mon des types d'élément (POUTRE DISCRET COQUE ...)
!       ACE_NB_(el)                 : Nombre de support dans le type
!       ACE_EL_(el)[ACE_NB_(el)]    : Liste des supports dans le type
!
! Nombre de type d'éléments
    integer(kind=8), parameter :: ACE_NB_ELEMENT = 9
    character(len=16), parameter :: ACE_NM_ELEMENT(ACE_NB_ELEMENT) = [ &
                                    'POUTRE          ', 'DISCRET         ', &
                                    'COQUE           ', 'CABLE           ', &
                                    'BARRE           ', 'GRILLE          ', &
                                    'MEMBRANE        ', 'MASSIF          ', &
                                    'MASSIF THM      ']
!
    integer(kind=8), parameter :: ACE_NU_POUTRE = 1
    integer(kind=8), parameter :: ACE_NU_DISCRET = 2
    integer(kind=8), parameter :: ACE_NU_COQUE = 3
    integer(kind=8), parameter :: ACE_NU_CABLE = 4
    integer(kind=8), parameter :: ACE_NU_BARRE = 5
    integer(kind=8), parameter :: ACE_NU_GRILLE = 6
    integer(kind=8), parameter :: ACE_NU_MEMBRANE = 7
    integer(kind=8), parameter :: ACE_NU_MASSIF = 8
    integer(kind=8), parameter :: ACE_NU_THHMM = 9
!
    integer(kind=8), parameter :: ACE_NB_POUTRE = 13
    character(len=16), parameter :: ACE_EL_POUTRE(ACE_NB_POUTRE) = [ &
                                    'MECA_POU_D_T    ', 'MECA_POU_D_E    ', &
                                    'MECA_POU_D_T_GD ', 'MEFS_POU_D_T    ', &
                                    'MECA_POU_D_TG   ', 'MECA_POHO_HEXA8 ', &
                                    'MECA_POHO_HEXA20', 'MET3SEG3        ', &
                                    'MET6SEG3        ', 'MET3SEG4        ', &
                                    'MECA_POU_D_EM   ', 'MECA_POU_D_TGM  ', &
                                    'MECA_POU_D_SQUE ']
!
!     integer, parameter :: ACE_MECA_POU_D_T         1
!     integer, parameter :: ACE_MECA_POU_D_E         2
!     integer, parameter :: ACE_MECA_POU_D_T_GD      3
!     integer, parameter :: ACE_MEFS_POU_D_T         4
!     integer, parameter :: ACE_MECA_POU_D_TG        5
!     integer, parameter :: ACE_MECA_POHO_HEXA8      6
!     integer, parameter :: ACE_MECA_POHO_HEXA20     7
!     integer, parameter :: ACE_MET3SEG3             8
!     integer, parameter :: ACE_MET6SEG3             9
!     integer, parameter :: ACE_MET3SEG4            10
!     integer, parameter :: ACE_MECA_POU_D_EM       11
!     integer, parameter :: ACE_MECA_POU_D_TGM      12
!     integer, parameter :: ACE_MECA_POU_D_SQUE     13
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_DISCRET = 8
    character(len=16), parameter :: ACE_EL_DISCRET(ACE_NB_DISCRET) = [ &
                                    'MECA_DIS_T_N    ', 'MECA_DIS_T_L    ', &
                                    'MECA_DIS_TR_N   ', 'MECA_DIS_TR_L   ', &
                                    'MECA_2D_DIS_T_N ', 'MECA_2D_DIS_T_L ', &
                                    'MECA_2D_DIS_TR_N', 'MECA_2D_DIS_TR_L']
!
!     integer, parameter :: ACE_MECA_DIS_T_N        1
!     integer, parameter :: ACE_MECA_DIS_T_L        2
!     integer, parameter :: ACE_MECA_DIS_TR_N       3
!     integer, parameter :: ACE_MECA_DIS_TR_L       4
!     integer, parameter :: ACE_MECA_2D_DIS_T_N     5
!     integer, parameter :: ACE_MECA_2D_DIS_T_L     6
!     integer, parameter :: ACE_MECA_2D_DIS_TR_N    7
!     integer, parameter :: ACE_MECA_2D_DIS_TR_L    8
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_COQUE = 27
    character(len=16), parameter :: ACE_EL_COQUE(ACE_NB_COQUE) = [ &
                                    'THCOTR3         ', 'THCOTR6         ', &
                                    'THCOQU4         ', 'THCOQU8         ', &
                                    'THCOTR7         ', 'THCOQU9         ', &
                                    'MEDKTR3         ', 'MEDSTR3         ', &
                                    'MET3TR3         ', 'MEDKQU4         ', &
                                    'MEDSQU4         ', 'MEQ4QU4         ', &
                                    'MECXSE3         ', 'MEDKTG3         ', &
                                    'MEDKQG4         ', 'MEQ4GG4         ', &
                                    'MET3GG3         ', 'THCASE3         ', &
                                    'THCPSE3         ', 'MEC3QU9H        ', &
                                    'MEC3TR7H        ', 'MEBODKT         ', &
                                    'MEBODST         ', 'MEBOQ4G         ', &
                                    'MEBOCQ3         ', 'THCOSE3         ', &
                                    'THCOSE2         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_CABLE = 2
    character(len=16), parameter :: ACE_EL_CABLE(ACE_NB_CABLE) = [ &
                                    'MECABL2         ', 'MEPOULI         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_BARRE = 3
    character(len=16), parameter :: ACE_EL_BARRE(ACE_NB_BARRE) = [ &
                                    'MECA_BARRE      ', 'MECA_2D_BARRE   ', &
                                    'MECGSEG3        ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_GRILLE = 6
    character(len=16), parameter :: ACE_EL_GRILLE(ACE_NB_GRILLE) = [ &
                                    'MEGCQU4         ', 'MEGMTR3         ', &
                                    'MEGMQU4         ', 'MEGMTR6         ', &
                                    'MEGMQU8         ', 'MEGCTR3         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_MEMBRANE = 6
    character(len=16), parameter :: ACE_EL_MEMBRANE(ACE_NB_MEMBRANE) = [ &
                                    'MEMBTR3         ', 'MEMBTR6         ', &
                                    'MEMBQU4         ', 'MEMBQU8         ', &
                                    'MEMBTR7         ', 'MEMBQU9         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_MASSIF = 97
    character(len=16), parameter :: ACE_EL_MASSIF(ACE_NB_MASSIF) = [ &
                                    'MECA_HEXA8      ', 'MECA_PENTA6     ', &
                                    'MECA_PENTA18    ', 'MECA_TETRA4     ', &
                                    'MECA_HEXA27     ', 'MECA_HEXA20     ', &
                                    'MECA_PENTA15    ', 'MECA_TETRA10    ', &
                                    'MECA_TETRS10    ', 'MECA_PYRAM5     ', &
                                    'MECA_PYRAM13    ', 'MECA_HEXS8      ', &
                                    'MECA_HEXS20     ', 'MEAXTR3         ', &
                                    'MEAXQU4         ', 'MEAXTR6         ', &
                                    'MEAXQU8         ', 'MEAXQU9         ', &
                                    'MEDPTR3         ', 'MEDPQU4         ', &
                                    'MEDPTR6         ', 'MEDPQU8         ', &
                                    'MEDPQU9         ', 'MECPTR3         ', &
                                    'MECPQU4         ', 'MECPTR6         ', &
                                    'MECPQU8         ', 'MECPQU9         ', &
                                    'THER_HEXA8      ', 'THER_HEXA8_D    ', &
                                    'THER_PENTA6_D   ', 'THER_TETRA4_D   ', &
                                    'THER_PENTA6     ', &
                                    'THER_TETRA4     ', 'THER_PYRAM5     ', &
                                    'THER_HEXA27     ', 'THER_HEXA20     ', &
                                    'THER_PENTA15    ', 'THER_TETRA10    ', &
                                    'THER_PYRAM13    ', 'THAXTR3         ', &
                                    'THAXQU4         ', 'THAXTR6         ', &
                                    'THAXQU8         ', 'THAXQU9         ', &
                                    'THPLTR3         ', 'THPLQU4         ', &
                                    'THPLTR6         ', 'THPLQU8         ', &
                                    'THPLQU9         ', 'MET3SEG3        ', &
                                    'MET6SEG3        ', 'MET3SEG4        ', &
                                    'MECA3DH27_HHO111', 'MECA3DT15_HHO111', &
                                    'MECA3DP21_HHO111', 'MECA3DP19_HHO111', &
                                    'MECA3DH27_HHO222', 'MECA3DT15_HHO222', &
                                    'MECA3DP21_HHO222', 'MECA3DP19_HHO222', &
                                    'MECA_DPQ9_HHO111', 'MECA_DPT7_HHO111', &
                                    'MECA_DPQ9_HHO222', 'MECA_DPT7_HHO222', &
                                    'THER3DH27_HHO000', 'THER3DT15_HHO000', &
                                    'THER3DP21_HHO000', 'THER3DP19_HHO000', &
                                    'THER3DH27_HHO111', 'THER3DT15_HHO111', &
                                    'THER3DP21_HHO111', 'THER3DP19_HHO111', &
                                    'THER3DH27_HHO222', 'THER3DT15_HHO222', &
                                    'THER3DP21_HHO222', 'THER3DP19_HHO222', &
                                    'THER2DQ9_HHO000 ', 'THER2DT7_HHO000 ', &
                                    'THER2DQ9_HHO111 ', 'THER2DT7_HHO111 ', &
                                    'THER2DQ9_HHO222 ', 'THER2DT7_HHO222 ', &
                                    'THERAXQ9_HHO000 ', 'THERAXT7_HHO000 ', &
                                    'THERAXQ9_HHO111 ', 'THERAXT7_HHO111 ', &
                                    'THERAXQ9_HHO222 ', 'THERAXT7_HHO222 ', &
                                    'MEEI_HEXA20     ', 'MEEI_HEXS20     ', &
                                    'MEEI_PENTA15    ', 'MEEI_PENTS15    ', &
                                    'EIAXQU8         ', 'EIAXQS8         ', &
                                    'EIPLQU8         ', 'EIPLQS8         ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_THHMM = 75
    character(len=16), parameter :: ACE_EL_THHMM(ACE_NB_THHMM) = [ &
                                    'HM_DPQ8S        ', 'HM_AXIS_QU8S    ', &
                                    'HM_DPTR6S       ', 'HM_AXIS_TR6S    ', &
                                    'HM_HEXA20S      ', 'HM_PENTA15S     ', &
                                    'HM_TETRA10S     ', 'THM_DPQ8S       ', &
                                    'THM_AXIS_QU8S   ', 'THM_DPTR6S      ', &
                                    'THM_AXIS_TR6S   ', 'THM_HEXA20S     ', &
                                    'THM_PENTA15S    ', 'THM_TETRA10S    ', &
                                    'H_DPQ8S         ', 'H_DPTR6S        ', &
                                    'H_HEXA20S       ', 'H_PENTA15S      ', &
                                    'H_TETRA10S      ', 'THHM_DPQ8S      ', &
                                    'THHM_AXIS_QU8S  ', 'THHM_DPTR6S     ', &
                                    'THHM_AXIS_TR6S  ', 'THHM_HEXA20S    ', &
                                    'THHM_PENTA15S   ', 'THHM_TETRA10S   ', &
                                    'HHM_DPQ8S       ', 'HHM_AXIS_QU8S   ', &
                                    'HHM_DPTR6S      ', 'HHM_AXIS_TR6S   ', &
                                    'HHM_HEXA20S     ', 'HHM_PENTA15S    ', &
                                    'HHM_TETRA10S    ', 'THH_DPQ8S       ', &
                                    'THH_AXIS_QU8S   ', 'THH_DPTR6S      ', &
                                    'THH_AXIS_TR6S   ', 'THH_HEXA20S     ', &
                                    'THH_PENTA15S    ', 'THH_TETRA10S    ', &
                                    'HH_DPQ8S        ', 'HH_AXIS_QU8S    ', &
                                    'HH_DPTR6S       ', 'HH_AXIS_TR6S    ', &
                                    'HH_HEXA20S      ', 'HH_PENTA15S     ', &
                                    'HH_TETRA10S     ', 'THH2M_DPQ8S     ', &
                                    'THH2M_AXIS_QU8S ', 'THH2M_DPTR6S    ', &
                                    'THH2M_AXIS_TR6S ', 'THH2M_HEXA20S   ', &
                                    'THH2M_PENTA15S  ', 'THH2M_TETRA10S  ', &
                                    'HH2M_DPQ8S      ', 'HH2M_AXIS_QU8S  ', &
                                    'HH2M_DPTR6S     ', 'HH2M_AXIS_TR6S  ', &
                                    'HH2M_HEXA20S    ', 'HH2M_PENTA15S   ', &
                                    'HH2M_TETRA10S   ', 'THH2_DPQ8S      ', &
                                    'THH2_AXIS_QU8S  ', 'THH2_DPTR6S     ', &
                                    'THH2_AXIS_TR6S  ', 'THH2_HEXA20S    ', &
                                    'THH2_PENTA15S   ', 'THH2_TETRA10S   ', &
                                    'HH2_DPQ8S       ', 'HH2_AXIS_QU8S   ', &
                                    'HH2_DPTR6S      ', 'HH2_AXIS_TR6S   ', &
                                    'HH2_HEXA20S     ', 'HH2_PENTA15S    ', 'HH2_TETRA10S    ']
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ACE_NB_TYPE_ELEM = ACE_NB_POUTRE+ACE_NB_DISCRET+ACE_NB_COQUE+ &
                                  ACE_NB_CABLE+ACE_NB_BARRE+ACE_NB_MASSIF+ &
                                  ACE_NB_GRILLE+ACE_NB_MEMBRANE+ACE_NB_THHMM

end module
