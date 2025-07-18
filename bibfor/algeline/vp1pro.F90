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
subroutine vp1pro(optiom, lraide, lmasse, ldynam, neq, &
                  nfreq, nfreqb, tolv, nitv, iexcl, &
                  fcorig, vec, resufi, resufr, resufk, &
                  nbrssa, nbpari, nbparr, nbpark, typres, &
                  optiof, solveu)
!     CALCUL DES VECTEURS ET VALEURS PROPRES PAR LA METHODE D'ITERATION
!     INVERSE.
!     ------------------------------------------------------------------
! IN  VALP : R8 : TABLEAU DES VALEURS PROPRES INITIALES
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
! ----------------------------------------------------------------------
!
! aslint: disable=W1504
    implicit none
!
! PARAMETRES D'APPEL
#include "jeveux.h"
#include "asterfort/freqom.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/omega2.h"
#include "asterfort/rectfr.h"
#include "asterfort/utmess.h"
#include "asterfort/vp1ite.h"
#include "asterfort/vpstur.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbpari, nbparr, nbpark, nitv, iexcl(*), nbrssa, lraide, lmasse
    integer(kind=8) :: ldynam, neq, nfreq, nfreqb, resufi(nfreqb, nbpari)
    real(kind=8) :: tolv, vec(neq, *), resufr(nfreqb, nbparr), fcorig
    character(len=*) :: optiom, resufk(nfreqb, nbpark)
    character(len=16) :: typres, optiof
    character(len=19) :: solveu
!
!
! VARIABLES LOCALES
    integer(kind=8) :: idet, place, irperm, lmx, lx0, iquoti, iprec, ifreq, ier, imode
    integer(kind=8) :: jfreq, iter, nbessa, i, iperm, j, kl, naux
    real(kind=8) :: det, err, eps, valeur, rperm
    character(len=24) :: kperm, cmulti, cvect0
!     ------------------------------------------------------------------
!
    call jemarq()
    cmulti = '&&VP1PRO.MODES.MULTIPLES'
    cvect0 = '&&VP1PRO.VECTEUR.TAMPON '
!
!     --- SAUVEGARDE DES TERMES MX POUR ORTHOGONALISER MODES MULTIPLES -
    call wkvect(cmulti, 'V V R', neq*nfreq, lmx)
    call wkvect(cvect0, 'V V R', neq, lx0)
!
    iquoti = 0
    if (optiom(1:8) .eq. 'RAYLEIGH') iquoti = 1
!     INITIALISATION A UNE VALEUR NON ATTEINTE
    iprec = -(nfreqb+1)
!
    do ifreq = 1, nfreq
        valeur = resufr(ifreq, 2)
20      continue
! --- POUR OPTIMISER ON NE CALCULE PAS LE DET
        call vpstur(lraide, valeur, lmasse, ldynam, det, &
                    idet, place, ier, solveu, .false._1, &
                    .true._1)
        if (ier .gt. 1) then
            valeur = 1.01d0*valeur
            goto 20
        end if
        if (resufi(ifreq, 1) .eq. 0) then
            iprec = -(nfreqb+1)
            imode = 1
            jfreq = ifreq
        else if (resufi(ifreq, 1) .eq. iprec) then
            imode = imode+1
        else
            iprec = resufi(ifreq, 1)
            imode = 1
            jfreq = ifreq
        end if
        call vp1ite(lmasse, lraide, ldynam, vec(1, jfreq), imode, &
                    valeur, neq, nitv, tolv, iter, &
                    zr(lx0), zr(lmx), err, iexcl, place, &
                    iquoti, solveu)
        if (resufi(ifreq, 1) .eq. 0) then
            resufi(ifreq, 1) = place
        else if (imode .gt. 1) then
            place = resufi(ifreq, 1)-imode+1
            resufi(ifreq, 1) = place
        end if
        resufr(ifreq, 2) = valeur
        resufi(ifreq, 4) = iter
        resufr(ifreq, 15) = err
    end do
!
! RECALCUL DU NUME_MODE POUR CHAQUE FREQUENCE
!
    if ((typres .eq. 'DYNAMIQUE') .and. (optiof .ne. 'PROCHE')) then
!
        ifreq = 1
30      continue
        valeur = resufr(ifreq, 2)
        if (abs(valeur) .le. omega2(fcorig)) then
            ifreq = ifreq+1
            if (ifreq .gt. nfreq) then
                call utmess('A', 'ALGELINE3_53')
            else
                goto 30
            end if
        end if
        nbessa = 0
        if (valeur .ge. 0.d0) then
            valeur = 0.95d0*valeur
        else
            valeur = 1.05d0*valeur
        end if
40      continue
! --- POUR OPTIMISER ON NE GARDE PAS LA FACTO (SI MUMPS) ET ON NE
! --- CALCULE PAS LE DET.
        call vpstur(lraide, valeur, lmasse, ldynam, det, &
                    idet, place, ier, solveu, .false._1, &
                    .false._1)
        if (ier .ne. 0) then
            nbessa = nbessa+1
            if (nbessa .gt. nbrssa) then
                call utmess('F', 'ALGELINE3_54')
            else
                if (valeur .ge. 0.d0) then
                    valeur = 0.95d0*valeur
                else
                    valeur = 1.05d0*valeur
                end if
                goto 40
            end if
        end if
        do ifreq = 1, nfreq
            resufr(ifreq, 2) = resufr(ifreq, 2)-valeur
        end do
!
    end if
!
! TRI DES VALEURS PROPRES SUIVANT ORDRE CROISSANT
! EN DYNAMIQUE: EN VALEUR ABSOLUE
! EN FLAMBEMENT: EN VALEUR ALGEBRIQUE
    eps = 1.d-7
    do i = 1, nfreq
        iperm = i
        if (typres .eq. 'DYNAMIQUE') then
            rperm = abs(resufr(i, 2))
        else
            rperm = resufr(i, 2)
        end if
        do j = i+1, nfreq
            if (typres .eq. 'DYNAMIQUE') then
                if (abs(resufr(j, 2)) .lt. (rperm*(1.d0-eps))) then
                    iperm = j
                    rperm = abs(resufr(iperm, 2))
                end if
                if ((abs(resufr(j, 2))-rperm) .le. (eps*rperm)) then
                    if (((resufr(j, 2)*resufr(iperm, 2)) .ge. 0.d0) .and. &
                        (abs(resufr(j, 2)) .lt. rperm)) then
                        iperm = j
                        rperm = abs(resufr(iperm, 2))
                    end if
                    if (((resufr(j, 2)*resufr(iperm, 2)) .lt. 0.d0) .and. &
                        (resufr(j, 2) .lt. 0.d0)) then
                        iperm = j
                        rperm = abs(resufr(iperm, 2))
                    end if
                end if
            else
                if (resufr(j, 2) .lt. rperm) then
                    iperm = j
                    rperm = resufr(iperm, 2)
                end if
            end if
        end do
!
! PERMUTATION DES DONNEES LIEES AUX VALEURS PROPRES
        if (iperm .ne. i) then
            do kl = 1, nbparr
                rperm = resufr(iperm, kl)
                resufr(iperm, kl) = resufr(i, kl)
                resufr(i, kl) = rperm
            end do
            do kl = 1, nbpari
                irperm = resufi(iperm, kl)
                resufi(iperm, kl) = resufi(i, kl)
                resufi(i, kl) = irperm
            end do
            do kl = 1, nbpark
                kperm = resufk(iperm, kl)
                resufk(iperm, kl) = resufk(i, kl)
                resufk(i, kl) = kperm
            end do
            do j = 1, neq
                rperm = vec(j, i)
                vec(j, i) = vec(j, iperm)
                vec(j, iperm) = rperm
            end do
        end if
    end do
!
    if ((typres .eq. 'DYNAMIQUE') .and. (optiof .ne. 'PROCHE')) then
!
        naux = 0
        call rectfr(nfreq, nfreq, valeur, place, naux, &
                    resufr(1, 2), nfreqb, resufi, resufr, nfreqb)
    else
        do ifreq = 1, nfreq
            resufi(ifreq, 1) = ifreq
        end do
    end if
!
    do ifreq = 1, nfreq
        resufr(ifreq, 1) = freqom(resufr(ifreq, 2))
        resufk(ifreq, 2) = 'INVERSE_R'
    end do
!
!     --- DESTRUCTION ZONE DE TRAVAIL ---
    call jedetr(cmulti)
    call jedetr(cvect0)
!
    call jedema()
end subroutine
