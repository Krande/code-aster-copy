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
subroutine te0335(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fgequi.h"
#include "asterfort/tecach.h"
    character(len=16) :: nomte, option
! person_in_charge: josselin.delmas at edf.fr
!
!     BUT:
!       CALCULER LES GRANDEURS EQUIVALENTES SUIVANTES
!       . CONTRAINTES EQUIVALENTES               (= 17 VALEURS)
!          . VON MISES                               (= 1 VALEUR)
!          . TRESCA                                  (= 1 VALEUR)
!          . CONTRAINTES PRINCIPALES                 (= 3 VALEURS)
!          . VON-MISES * SIGNE (PRESSION)            (= 1 VALEUR)
!          . DIRECTIONS DES CONTRAINTES PRINCIPALES  (= 3*3 VALEURS)
!          . TRACE                                   (= 1 VALEUR)
!          . TAUX DE TRIAXIALITE                     (= 1 VALEUR)
!
!       . DEFORMATIONS EQUIVALENTES              (= 14 VALEURS)
!          . SECOND INVARIANT                        (= 1 VALEUR)
!          . DEFORMATIONS PRINCIPALES                (= 3 VALEURS)
!          . 2EME INV. * SIGNE (1ER.INV.)            (= 1 VALEUR)
!          . DIRECTIONS DES DEFORMATIONS PRINCIPALES (= 3*3 VALEURS)
!
!       AUX POINTS DE GAUSS ET AUX NOEUDS :
!       A PARTIR DE SIGM_ELGA ET SIGM_ELNO POUR LES CONTRAINTES
!       A PARTIR DE EPSI_ELGA ET EPSI_ELNO POUR LES DEFORMATIONS
!       A PARTIR DE EPME_ELGA ET EPME_ELNO POUR LES DEF. HORS THERMIQUE
!
!       OPTION : 'SIEQ_ELGA'
!                'SIEQ_ELNO'
!                'EPEQ_ELGA'
!                'EPEQ_ELNO'
!                'EPMQ_ELGA'
!                'EPMQ_ELNO'
!                'EPGQ_ELGA'
!                'EPGQ_ELNO'
!
! ----------------------------------------------------------------------
!
!
!
    integer(kind=8) :: neeqmx, nceqmx
    parameter(neeqmx=14, nceqmx=17)
!
    integer(kind=8) :: ndim, ndim1, nno, nnos, npg, ipoids, ivf, idfde, jgano
    integer(kind=8) :: idefo, icont, iequi
    integer(kind=8) :: iret, itabin(7), itabou(7), nbcmp, ncmpeq, nbsp
    integer(kind=8) :: idec, ideceq, ipg, ino, isp
!
!
! ----------------------------------------------------------------------
!
!
    if ((nomte .eq. 'MEC3QU9H') .or. (nomte .eq. 'MEC3TR7H')) then
        call elrefe_info(fami='MASS', ndim=ndim1, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    else
        call elrefe_info(fami='RIGI', ndim=ndim1, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    end if
!
    if ((option .eq. 'EPEQ_ELGA') .or. (option .eq. 'EPEQ_ELNO') .or. (option .eq. 'EPMQ_ELGA') &
        .or. (option .eq. 'EPMQ_ELNO') .or. (option .eq. 'EPGQ_ELGA') .or. &
        (option .eq. 'EPGQ_ELNO')) then
!
        call tecach('OOO', 'PDEFORR', 'L', iret, nval=7, &
                    itab=itabin)
        idefo = itabin(1)
        call tecach('OOO', 'PDEFOEQ', 'E', iret, nval=7, &
                    itab=itabou)
        ASSERT(itabou(2)/itabou(3) .eq. neeqmx)
!
    elseif ((option .eq. 'SIEQ_ELGA') .or. (option .eq. 'SIEQ_ELNO')) &
        then
!
        call tecach('OOO', 'PCONTRR', 'L', iret, nval=7, &
                    itab=itabin)
        icont = itabin(1)
        call tecach('OOO', 'PCONTEQ', 'E', iret, nval=7, &
                    itab=itabou)
        ASSERT(itabou(2)/itabou(3) .eq. nceqmx)
!
    else
        ASSERT(.false.)
    end if
!
    iequi = itabou(1)
!
    nbsp = itabou(7)
    ASSERT(nbsp .ge. 1)
    ASSERT(nbsp .eq. itabin(7))
!
    nbcmp = itabin(2)/itabin(3)
    ASSERT((nbcmp .eq. 1) .or. (nbcmp .eq. 4) .or. (nbcmp .eq. 6))
!
    ncmpeq = itabou(2)/itabou(3)
    ASSERT((ncmpeq .eq. neeqmx) .or. (ncmpeq .eq. nceqmx))
!
    ASSERT(itabin(6) .le. 1)
    ASSERT(itabou(6) .le. 1)
!
    if (nbcmp .eq. 6) then
        ndim = 3
    else if (nbcmp .eq. 4) then
        ndim = 2
    else if (nbcmp .eq. 1) then
        ndim = 1
    end if
!
! ----------------------------------------------------------------
! --- DEFORMATIONS ET CONTRAINTES EQUIVALENTES AUX POINTS DE GAUSS
! ----------------------------------------------------------------
!
    if (option(6:9) .eq. 'ELGA') then
!
! ------ DEFORMATIONS :
! -------------------
        if ((option .eq. 'EPEQ_ELGA') .or. (option .eq. 'EPMQ_ELGA') .or. &
            (option .eq. 'EPGQ_ELGA')) then
            do ipg = 1, npg
                do isp = 1, nbsp
                    idec = idefo+(ipg-1)*nbcmp*nbsp+(isp-1)*nbcmp
                    ideceq = iequi+(ipg-1)*ncmpeq*nbsp+(isp-1)*ncmpeq
                    call fgequi(zr(idec), 'EPSI_DIR', ndim, zr(ideceq))
                end do
            end do
!
! ----- CONTRAINTES :
! -----------------
        else if (option .eq. 'SIEQ_ELGA') then
            do ipg = 1, npg
                do isp = 1, nbsp
                    idec = icont+(ipg-1)*nbcmp*nbsp+(isp-1)*nbcmp
                    ideceq = iequi+(ipg-1)*ncmpeq*nbsp+(isp-1)*ncmpeq
                    call fgequi(zr(idec), 'SIGM_DIR', ndim, zr(ideceq))
                end do
            end do
        end if
!
! -------------------------------------------------------
! --- DEFORMATIONS ET CONTRAINTES EQUIVALENTES AUX NOEUDS
! -------------------------------------------------------
!
    else if (option(6:9) .eq. 'ELNO') then
!
! ------ DEFORMATIONS :
! -------------------
        if ((option .eq. 'EPEQ_ELNO') .or. (option .eq. 'EPMQ_ELNO') .or. &
            (option .eq. 'EPGQ_ELNO')) then
            do ino = 1, nno
                do isp = 1, nbsp
                    idec = idefo+(ino-1)*nbcmp*nbsp+(isp-1)*nbcmp
                    ideceq = iequi+(ino-1)*ncmpeq*nbsp+(isp-1)*ncmpeq
                    call fgequi(zr(idec), 'EPSI_DIR', ndim, zr(ideceq))
                end do
            end do
!
! ----- CONTRAINTES :
! -----------------
        else if (option .eq. 'SIEQ_ELNO') then
            do ino = 1, nno
                do isp = 1, nbsp
                    idec = icont+(ino-1)*nbcmp*nbsp+(isp-1)*nbcmp
                    ideceq = iequi+(ino-1)*ncmpeq*nbsp+(isp-1)*ncmpeq
                    call fgequi(zr(idec), 'SIGM_DIR', ndim, zr(ideceq))
                end do
            end do
        end if
!
    end if
!
!
end subroutine
