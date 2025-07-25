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

subroutine eclaty(nomte, elrefa, fapg, npg, npoini, &
                  nterm1, nsomm1, csomm1, tyma, nbno2, &
                  connx, mxnbn2, mxnbpi, mxnbte, mxnbse, &
                  nbsel, corsel, iret)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/ecla2d.h"
#include "asterfort/ecla3d.h"
#include "asterfort/elraca.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: mxnbn2, mxnbpi, mxnbte, mxnbse
    integer(kind=8) :: ndim, npg, connx(mxnbn2, mxnbse), nsomm1(mxnbpi, mxnbte)
    integer(kind=8) :: nterm1(mxnbpi), nbno2(mxnbse), npoini, tyma(mxnbse)
    integer(kind=8) :: nbsel, corsel(mxnbse), iret
    real(kind=8) :: csomm1(mxnbpi, mxnbte)
    character(len=8) :: elrefa, fapg
    character(len=16) :: nomte
!
! BUT : DECOMPOSER LES TYPE_ELEM EN AUTANT DE SOUS-ELEMENTS QUE
!       DE POINTS DE GAUSS.
! ---------------------------------------------------------------------
! IN : NOMTE  : NOM D'UN TYPE_ELEM
! IN : ELREFA : NOM DE L'ELREFA
! IN : FAPG   : FAMILLE DE POINTS DE GAUSS
!      MXNBN2 : MAX DU NOMBRE DE NOEUDS D'UN SOUS-ELEMENT (HEXA8)
!      MXNBPI : MAX DU NOMBRE DE POINT_I (HEXA A 27 POINTS DE GAUSS)
!      MXNBTE : MAX DU NOMBRE DE TERMES DE LA C.L. DEFINISSANT 1 POINT_I
!      MXNBSE : MAX DU NOMBRE DE SOUS-ELEMENTS
!
! OUT : NPG    : NOMBRE DE POINTS DE GAUSS (PG)
! OUT : NBSEL  : NOMBRE DE SOUS-ELEMENTS (SE)
! OUT : CORSEL : CORRESPONDANCE SE -> KPG
!       EN GENERAL : NBSE=NBPG ET CORSEL(KSE)=KPG
!       MAIS PARFOIS :
!           NBSEL > NBPG : PLUSIEURS SE POUR 1 SEUL PG
!           NBSEL < NBPG : IL Y A DES PG SANS SE (DEHORS)
! OUT : IRET = 0 -> le type_elem a été traité
!            = 1 -> non traité
!
! ---------------------------------------------------------------------
! DESCRIPTION DES POINTS INTERMEDIAIRES (POINT_I) :
! ------------------------------------------------
! UN POINT_I EST DEFINI COMME UNE COMBINAISON LINEAIRE DES NOEUDS
! DE LA MAILLE SOUS-JACENTE AU TYPE_ELEM :
! POINT_I = SOMME COEF(K)*NOEUD(K)  (1<= K <=NTERMES)
!           NTERMES <= 27 (HEXA27)
!
! OUT : NPOINI : NOMBRE DE POINT_I POUR LE TYPE_ELEM/FAPG
!                    (NPOINI <= 64 : (HEXA20 A 27 POINTS DE GAUSS))
! OUT : NTERM1 : NTERM1(IPOINI) : NOMBRE DE TERMES DE LA COMBINAISON
!                POUR LE POINT_I IPOINI
! OUT : NSOMM1 : NSOMM1(IPOINI,K) NUMERO DU NOEUD DU KEME TERME.
! OUT : CSOMM1 : CSOMM1(IPOINI,K) COEF. DU NOEUD DU KEME TERME.
!
! CONNECTIVITE DES SOUS-ELEMENTS :
!-------------------------------
! OUT : TYMA (I)  : TYMA(KSE)-> TYPE_MAILLE ASSOCIE AU SOUS-ELEMENT KSE
! OUT : NBNO2 (I) : NBNO2(KSE) -> NNO2
!                    NNO2 : NOMBRE DE NOEUDS DU SOUS-ELEMENT KSE
! OUT : CONNX (I) : CONNX(INO2,KSE) -> IPOINI
!                    IPOINI= NUMERO DU POINT_I ASSOCIE AU
!                    SOMMET INO2 DU SOUS-ELEMENT KSE
!                    (1<= INO2 <= NNO2(KSE))
!                    (1<= IPOINI <= NPOINI)
!
! ---------------------------------------------------------------------
    integer(kind=8) :: nbfpg, nbpg(MT_NBFAMX), nufpg
    character(len=8) :: famg(MT_NBFAMX)
! ---------------------------------------------------------------------
    call jemarq()
!
    npg = 0
    npoini = 0
    nbsel = 0
    iret = 1

! - Get list of integration schemes of geometric support
    call elraca(elrefa, &
                nbfpg_=nbfpg, fapg_=famg, nbpg_=nbpg, &
                ndim_=ndim)

! - Get index for integration scheme
    nufpg = indik8(famg, fapg, 1, nbfpg)
    ASSERT(nufpg .gt. 0)

    npg = nbpg(nufpg)
!
    if (ndim .eq. 2) then
        call ecla2d(nomte, elrefa, fapg, npg, npoini, &
                    nterm1, nsomm1, csomm1, tyma, nbno2, &
                    connx, mxnbn2, mxnbpi, mxnbte, mxnbse, &
                    nbsel, corsel)
        iret = 0
    else if (ndim .eq. 3) then
        call ecla3d(nomte, elrefa, fapg, npg, npoini, &
                    nterm1, nsomm1, csomm1, tyma, nbno2, &
                    connx, mxnbn2, mxnbpi, mxnbte, mxnbse, &
                    nbsel, corsel)
        iret = 0
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jedema()
!
end subroutine
