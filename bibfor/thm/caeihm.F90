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
! aslint: disable=W1504
!
subroutine caeihm(ds_thm, nomte, l_axi, mecani, press1, &
                  press2, tempe, dimdef, dimcon, ndim, &
                  nno1, nno2, npi, npg, dimuel, &
                  iw, ivf1, idf1, ivf2, idf2, &
                  jgano1, iu, ip, ipf, iq, &
                  inte_type)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/greihm.h"
#include "asterfort/lteatt.h"
#include "asterfort/thmGetElemIntegration.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    character(len=3), intent(out) :: inte_type
!
! --------------------------------------------------------------------------------------------------
!
! PREPARATION DU CALCUL SUR UN ELEMENT DE JOINT HM
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! IN NOMTE   : NOM DU TYPE D'ELEMENT
! IN AXI     : AXI ?
! OUT PERMAN : MODELISATION HM PERMAMENTE ?
! OUT MECANI : TABLEAU INFO SUR MECANIQUE
! OUT PRESS1 : TABLEAU INFO SUR HYDRAULIQUE CONSTITUANT 1
! OUT PRESS2 : TABLEAU INFO SUR HYDRAULIQUE CONSTITUANT 2
! OUT TEMPE  : TABLEAU INFO SUR THERMIQUE
! OUT DIMDEF : DIMENSION DES DEFORMATIONS GENERALISEES
! OUT DIMCON : DIMENSION DES CONTRAINTES GENERALISEES
! OUT NDIM   : DIMENSION DU PROBLEME (2 OU 3)
! OUT NNO1   : NOMBRE DE NOEUDS DES BORDS INF ET DUP DE L'ELEMENT
! OUT NNO2   : NOMBRE DE NOEUDS DU SEGMENT CENTRAL
! OUT NPI    : NOMBRE DE POINTS D'INTEGRATION DE L'ELEMENT
! OUT NPG    : NOMBRE DE POINTS DE GAUSS
! OUT DIMUEL : NOMBRE DE DDL TOTAL DE L'ELEMENT
! OUT IW     : ADRESSE DU TABLEAU POIDS POUR FONCTION DE FORME P2
! OUT IVF1   : ADRESSE DU TABLEAU DES FONCTIONS DE FORME P2
! OUT IDF1   : ADRESSE DU TABLEAU DES DERIVESS DES FONCTIONS DE FORME P2
! OUT IVF2   : ADRESSE DU TABLEAU DES FONCTIONS DE FORME P1
! OUT IDF2 : ADRESSE DU TABLEAU DES DERIVESS DES FONCTIONS DE FORME P1
! OUT JGANO1  : ADRESSE DANS ZR DE LA MATRICE DE PASSAGE
!              GAUSS -> NOEUDS
! OUT IU     : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! OUT IP     : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE PRESSION MILIEU
! OUT IPF    : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE PRESSION FACES
! OUT IQ     : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE HYDRO
! OUT MODINT : MODE D'INTEGRATION
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_axi, l_vf
    integer(kind=8) :: mecani(8), press1(9), press2(9), tempe(5), dimuel
    integer(kind=8) :: ndim, nnos, nno1, nno2, ntrou
    integer(kind=8) :: dimdef, dimcon
    integer(kind=8) :: npg, npi, n, i
    integer(kind=8) :: ivf1, idf1, ivf2, idf2, jgano1, jgano2, iw
    integer(kind=8) :: iu(3, 18), ip(2, 9), ipf(2, 2, 9), iq(2, 2, 9)
    integer(kind=8), parameter :: f1q8(6) = (/1, 2, 5, 4, 3, 7/)
    integer(kind=8), parameter :: f2q8(2) = (/8, 6/)
    integer(kind=8), parameter :: f3q8(2) = (/1, 2/)
    integer(kind=8), parameter :: f4q8(2) = (/4, 3/)
    character(len=8) :: lielrf(10)
    character(len=16) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
    ndim = 2
    l_axi = ASTER_FALSE
    l_vf = ASTER_FALSE
! ======================================================================
! --- INITIALISATION DES GRANDEURS GENERALISEES SELON MODELISATION -----
! ======================================================================
    call greihm(ndim, mecani, press1, press2, &
                tempe, dimdef, dimcon)

!
! - Flag for *JHMS elements
!
    ds_thm%ds_elem%l_jhms = ASTER_TRUE
!
! - Get type of integration
!
    call thmGetElemIntegration(l_vf, inte_type=inte_type)
!
    l_axi = lteatt('AXIS', 'OUI')
! ======================================================================
! --- ADAPTATION AU MODE D'INTEGRATION ---------------------------------
! --- DEFINITION DE L'ELEMENT (NOEUDS, SOMMETS, POINTS DE GAUSS) -------
! ======================================================================
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, nnos=nnos, &
                     npg=npi, jpoids=iw, jvf=ivf1, jdfde=idf1, jgano=jgano1)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nno2, nnos=nnos, &
                     npg=npi, jpoids=iw, jvf=ivf2, jdfde=idf2, jgano=jgano2)
!
    if (inte_type .eq. 'RED') then
        npg = npi-nnos
    end if
    if (inte_type .eq. 'CLA') then
        npg = npi
    end if
!
    ndim = ndim+1
!
! ======================================================================
! --- DETERMINATION DES DECALAGES D'INDICE POUR ACCEDER AUX DDL --------
! ======================================================================
!
    dimuel = 2*nno1*ndim+nno2*3*(press1(1)+press2(1))+2
    do n = 1, 5
        do i = 1, 2
            iu(i, n) = i+(f1q8(n)-1)*3
        end do
    end do
    do i = 1, 2
        iu(i, 6) = iu(i, 3)+4
    end do
    do n = 1, 2
        ip(1, n) = 16+(f2q8(n)-6)*2
    end do
    do n = 1, 2
        ipf(1, 1, n) = 3+(f4q8(n)-1)*3
    end do
    do n = 1, 2
        ipf(1, 2, n) = 3+(f3q8(n)-1)*3
    end do
    iq(1, 1, 1) = iu(2, 6)+1
    iq(1, 2, 1) = iu(2, 3)+1
!
end subroutine
