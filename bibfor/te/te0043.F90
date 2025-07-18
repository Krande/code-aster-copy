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
subroutine te0043(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/infdis.h"
#include "asterfort/infted.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2plg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpplg.h"
#include "asterfort/utpslg.h"
!
    character(len=*) :: option, nomte
!     CALCUL DES CHARGES DE PESANTEUR DANS LES ELEMENTS DISCRETS
!     ------------------------------------------------------------------
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'CHAR_MECA_PESA_R'  : CALCUL DES CHARGES DE PESANTEUR
!        'CHAR_MECA_ROTA_R'  : CALCUL DES CHARGES DE ROTATION
! IN  NOMTE  : K16 : NOM DU TYPE D'ELEMENT DISCRET :
!         MECA_DIS_T_N
!         MECA_DIS_T_L
!         MECA_DIS_TR_N
!         MECA_DIS_TR_L
!         MECA_2D_DIS_T_N
!         MECA_2D_DIS_T_L
!         MECA_2D_DIS_TR_N
!         MECA_2D_DIS_TR_L
!     ------------------------------------------------------------------
!
!      CHARACTER*32   JEXNUM,JEXNOM,JEXR8,JEXATR
!
    real(kind=8) :: pgl(3, 3), mat(144), r8bid, omega
    integer(kind=8) :: infodi, lrota, ldis, lorien, lpesa, lvectu, itype
    integer(kind=8) :: nbterm, nno, nc, ndim, irep, i, ibid
    character(len=8) :: k8bid
    character(len=20) :: kmess(5)
!     ------------------------------------------------------------------
!
    if ((option .ne. 'CHAR_MECA_PESA_R') .and. (option .ne. 'CHAR_MECA_ROTA_R')) then
        kmess(1) = option
        kmess(2) = nomte
        kmess(3) = 'TE0043'
        call utmess('F', 'DISCRETS_14', nk=3, valk=kmess)
    end if
!
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
!     DISCRET DE TYPE MASSE
    call infdis('DISM', infodi, r8bid, k8bid)
    if (infodi .eq. 0) then
        call utmess('A+', 'DISCRETS_26', sk=nomte)
        call infdis('DUMP', ibid, r8bid, 'A+')
    end if
!
!     [M] : SYMETRIQUE ?
    call infdis('SYMM', infodi, r8bid, k8bid)
! --- INFORMATIONS SUR LES DISCRETS :
!        NBTERM   =  NOMBRE DE COEFFICIENTS DANS K
!        NNO      =  NOMBRE DE NOEUDS
!        NC       =  NOMBRE DE COMPOSANTE PAR NOEUD
!        NDIM     =  DIMENSION DE L'ELEMENT
!        ITYPE    =  TYPE DE L'ELEMENT
    call infted(nomte, infodi, nbterm, nno, nc, &
                ndim, itype)
!
!     --- CALCUL DES VECTEURS ELEMENTAIRES ----
    if (option .eq. 'CHAR_MECA_PESA_R') then
!        --- MATRICE DE RIGIDITE LOCALE ---
        call jevech('PCADISM', 'L', ldis)
        call jevech('PCAORIE', 'L', lorien)
        call matrot(zr(lorien), pgl)
!
!        --- GLOBAL VERS LOCAL ? ---
!        --- IREP = 2 = MATRICE EN REPERE LOCAL ==> PASSER EN GLOBAL ---
        call infdis('REPM', irep, r8bid, k8bid)
        if (irep .eq. 2) then
            if (ndim .eq. 3) then
                if (infodi .eq. 1) then
                    call utpslg(nno, nc, pgl, zr(ldis), mat)
                else if (infodi .eq. 2) then
                    call utpplg(nno, nc, pgl, zr(ldis), mat)
                end if
            else if (ndim .eq. 2) then
                if (infodi .eq. 1) then
                    call ut2mlg(nno, nc, pgl, zr(ldis), mat)
                else if (infodi .eq. 2) then
                    call ut2plg(nno, nc, pgl, zr(ldis), mat)
                end if
            end if
        else
            do i = 1, nbterm
                mat(i) = zr(ldis+i-1)
            end do
        end if
!
!        --- CHAMP DE PESANTEUR ---
        call jevech('PPESANR', 'L', lpesa)
!        --- VECTEUR CHARGEMENT ---
        call jevech('PVECTUR', 'E', lvectu)
!
!        --- ON Y VA ---
        lvectu = lvectu-1
        if (nomte .eq. 'MECA_DIS_TR_L') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(06)*zr(lpesa)*zr(lpesa+3)
                zr(lvectu+7) = mat(28)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+8) = mat(36)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+9) = mat(45)*zr(lpesa)*zr(lpesa+3)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(14)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(27)*zr(lpesa)*zr(lpesa+3)
                zr(lvectu+7) = mat(79)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+8) = mat(92)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+9) = mat(105)*zr(lpesa)*zr(lpesa+3)
            end if
        else if (nomte .eq. 'MECA_DIS_TR_N') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(06)*zr(lpesa)*zr(lpesa+3)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(08)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(15)*zr(lpesa)*zr(lpesa+3)
            end if
        else if (nomte .eq. 'MECA_DIS_T_L') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(06)*zr(lpesa)*zr(lpesa+3)
                zr(lvectu+4) = mat(10)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+5) = mat(15)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+6) = mat(21)*zr(lpesa)*zr(lpesa+3)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(08)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(15)*zr(lpesa)*zr(lpesa+3)
                zr(lvectu+4) = mat(22)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+5) = mat(29)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+6) = mat(36)*zr(lpesa)*zr(lpesa+3)
            end if
        else if (nomte .eq. 'MECA_DIS_T_N') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(06)*zr(lpesa)*zr(lpesa+3)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(05)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(09)*zr(lpesa)*zr(lpesa+3)
            end if
        else if (nomte .eq. 'MECA_2D_DIS_TR_L') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+4) = mat(10)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+5) = mat(15)*zr(lpesa)*zr(lpesa+2)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(08)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+4) = mat(22)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+5) = mat(29)*zr(lpesa)*zr(lpesa+2)
            end if
        else if (nomte .eq. 'MECA_2D_DIS_TR_N') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(05)*zr(lpesa)*zr(lpesa+2)
            end if
        else if (nomte .eq. 'MECA_2D_DIS_T_L') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(06)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+4) = mat(10)*zr(lpesa)*zr(lpesa+2)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(06)*zr(lpesa)*zr(lpesa+2)
                zr(lvectu+3) = mat(11)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+4) = mat(16)*zr(lpesa)*zr(lpesa+2)
            end if
        else if (nomte .eq. 'MECA_2D_DIS_T_N') then
            if (infodi .eq. 1) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(03)*zr(lpesa)*zr(lpesa+2)
            else if (infodi .eq. 2) then
                zr(lvectu+1) = mat(01)*zr(lpesa)*zr(lpesa+1)
                zr(lvectu+2) = mat(04)*zr(lpesa)*zr(lpesa+2)
            end if
        else
            kmess(1) = option
            kmess(2) = nomte
            kmess(3) = 'TE0043'
            call utmess('F', 'DISCRETS_15', nk=3, valk=kmess)
        end if
    end if
!
!
    if (option .eq. 'CHAR_MECA_ROTA_R') then
!        DISCRET DE TYPE MASSE ?
        call infdis('DISM', infodi, r8bid, k8bid)
!        SI C'EST UN DISCRET DE MASSE :
!           L'OPTION N'EST PAS DEVELOPPEE SAUF SI OMEGA=0
!        DANS TOUS LES CAS ON RENVOI UN EFFORT NUL
        if (infodi .eq. 1) then
            call jevech('PROTATR', 'L', lrota)
            omega = zr(lrota)
            if (abs(omega) .ge. r8prem()) then
                kmess(1) = nomte
                kmess(2) = option
                call utmess('F', 'CALCUL_37', nk=2, valk=kmess)
            end if
        end if
        call jevech('PVECTUR', 'E', lvectu)
        do i = 1, nno*nc
            zr(lvectu-1+i) = 0.0d0
        end do
    end if
!
end subroutine
