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
!
subroutine ef0042(nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/infdis.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
#include "asterfort/ut2mgl.h"
#include "asterfort/ut2pgl.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/utmess.h"
#include "asterfort/utppgl.h"
#include "asterfort/utpsgl.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "asterfort/vecmap.h"
!
    character(len=16) :: nomte
!     CALCUL DE EFGE_ELNO
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbterm, nno, nc, neq, irep, i, ndim, ibid
    integer(kind=8) :: ldis, lorien, jeffo, jdepl, infodi, itype
    real(kind=8) :: ulr(12), flr(12)
    real(kind=8) :: pgl(3, 3), klc(144), mat(144), r8bid
    character(len=8) :: k8bid
!     ------------------------------------------------------------------
!
!     ON VERIFIE QUE LES CARACTERISTIQUES ONT ETE AFFECTEES
!     LE CODE DU DISCRET
    call infdis('CODE', ibid, r8bid, nomte)
!     LE CODE STOKE DANS LA CARTE
    call infdis('TYDI', infodi, r8bid, k8bid)
    if (infodi .ne. ibid) then
        call utmess('F+', 'DISCRETS_25', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'F+')
    end if
!     DISCRET DE TYPE RAIDEUR
    call infdis('DISK', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_27', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    end if
!
!     MATRICE DE RAIDEUR SYMETRIQUE OU PAS
    infodi = 1
    call infdis('SYMK', infodi, r8bid, k8bid)
! --- INFORMATIONS SUR LES DISCRETS :
!        NBTERM  = NOMBRE DE COEFFICIENTS DANS K
!        NNO     = NOMBRE DE NOEUDS
!        NC      = NOMBRE DE COMPOSANTE PAR NOEUD
!        NDIM    = DIMENSION DE L'ELEMENT
!        ITYPE = TYPE DE L'ELEMENT
    call infted(nomte, infodi, nbterm, nno, nc, &
                ndim, itype)
    neq = nno*nc
!
!     --- MATRICE DE RIGIDITE LOCALE ---
    call jevech('PCADISK', 'L', ldis)
    call jevech('PCAORIE', 'L', lorien)
    call matrot(zr(lorien), pgl)
!
!     --- ABSOLU VERS LOCAL ? ---
!     --- IREP = 1 = MATRICE EN REPERE GLOBAL ==> PASSER EN LOCAL ---
    call infdis('REPK', irep, r8bid, k8bid)
    if (irep .eq. 1) then
        if (ndim .eq. 3) then
            if (infodi .eq. 1) then
                call utpsgl(nno, nc, pgl, zr(ldis), mat)
            else if (infodi .eq. 2) then
                call utppgl(nno, nc, pgl, zr(ldis), mat)
            end if
        else if (ndim .eq. 2) then
            if (infodi .eq. 1) then
                call ut2mgl(nno, nc, pgl, zr(ldis), mat)
            else if (infodi .eq. 2) then
                call ut2pgl(nno, nc, pgl, zr(ldis), mat)
            end if
        end if
    else
        do i = 1, nbterm
            mat(i) = zr(ldis+i-1)
        end do
    end if
!
!     ---- MATRICE RIGIDITE LIGNE > MATRICE RIGIDITE CARRE
    if (infodi .eq. 1) then
        call vecma(mat, nbterm, klc, neq)
    else if (infodi .eq. 2) then
        call vecmap(mat, nbterm, klc, neq)
    end if
!
!     --- CALCUL DES VECTEURS ELEMENTAIRES ----
!
    call jevech('PEFFORR', 'E', jeffo)
    call jevech('PDEPLAR', 'L', jdepl)
!
!        --- VECTEUR DEPLACEMENT LOCAL  ULR = PGL * UG  ---
    if (ndim .eq. 3) then
        call utpvgl(nno, nc, pgl, zr(jdepl), ulr)
    else if (ndim .eq. 2) then
        call ut2vgl(nno, nc, pgl, zr(jdepl), ulr)
    end if
!
!        --- VECTEUR EFFORT      LOCAL  FLR = KLC * ULR  ---
    call pmavec('ZERO', neq, klc, ulr, flr)
!
!     ON CHANGE LE SIGNE DES EFFORTS SUR LE PREMIER NOEUD, POUR LES
!     ELEMENTS A 2 NOEUDS
    if (nno .eq. 1) then
        do i = 1, neq
            zr(jeffo+i-1) = flr(i)
        end do
    else
        do i = 1, nc
            zr(jeffo+i-1) = -flr(i)
            zr(jeffo+i+nc-1) = flr(i+nc)
        end do
    end if
end subroutine
