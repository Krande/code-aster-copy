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
subroutine vpinte(option, nfreq, valp, det, idet, &
                  ieme, npas, tolf, nitf, lraide, &
                  lmasse, ldynam, resufi, resufr, nfreqb, &
                  solveu)
    implicit none
#include "asterfort/vpstur.h"
    character(len=16) :: option
    integer(kind=8) :: idet(*), ieme(*), npas(*)
    real(kind=8) :: tolf
    integer(kind=8) :: nfreqb
    integer(kind=8) :: lraide, lmasse, ldynam, resufi(nfreqb, *)
    real(kind=8) :: valp(*), det(*), resufr(nfreqb, *)
    character(len=19) :: solveu
!     METHODE D'INTERPOLATION DE VALEUR PROPRE COMPRISE ENTRE 2 BORNES
!     ------------------------------------------------------------------
! IN OPTION  : CH8 : OPTION D'INTERPOLATION
!          SEPARE ==> ON PREND LE MILIEU DE L'INTERVALLE
!          AJUSTE ==> ON EFFECTUE DES ITERATIONS DE TYPE SECANTE
! IN TOLF    : R8  : TOLERANCE POUR L'AJUSTEMENT DE FREQUENCE
! IN NITF    : R8  : NOMBRE MAXIMUM D'ITERATION DE LA METHODE
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
!     ------------------------------------------------------------------
!     METHODE D'INTERPOLATION DE DEGRE 1 (A DEUX POINTS)
!
!                   ZK * F(ZK-1) - ZK-1 * F(ZK)
!          ZK+1  =  ---------------------------
!                   WK-1 * F(ZK-1) - WK * F(ZK)
!     ------------------------------------------------------------------
!
    real(kind=8) :: det0, om0, om, preci, ftol
    real(kind=8) :: dd1, d1, f1, w1
    real(kind=8) :: dd2, d2, f2, w2
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: idec, idet0, ieme0, iemei, ier, ifreq, ii
    integer(kind=8) :: ip1, ip2, itrig, jdec, jfreq, nbiter, nfreq
    integer(kind=8) :: niter, nitf, ntrou
    real(kind=8) :: valeur
!-----------------------------------------------------------------------
    ifreq = 0
    det0 = 0.d0
    idet0 = 0
    if (option(1:6) .eq. 'AJUSTE') then
!
        do jfreq = 1, nfreq-1
            iemei = ieme(jfreq+1)-ieme(jfreq)
            if (iemei .eq. 0) goto 1
            if (valp(jfreq) .eq. 0.0d0) then
!              --- FREQUENCE NULLE  ---
                nbiter = 0
                valeur = valp(jfreq)
                preci = 0.0d0
            else if (abs(iemei) .gt. 1) then
!              --- FREQUENCES "MULTIPLES"  ---
                nbiter = 0
                valeur = (valp(jfreq)+valp(jfreq+1))*0.5d0
                preci = (valp(jfreq+1)-valp(jfreq))/valeur
                if (jfreq .eq. 1) then
                    do ii = 1, abs(iemei)-1
                        ifreq = ifreq+1
                        resufi(ifreq, 1) = ieme(jfreq+1)
                        resufi(ifreq, 2) = max(npas(jfreq), npas(jfreq+1))
                        resufi(ifreq, 3) = nbiter
                        resufi(ifreq, 4) = nbiter
                        resufr(ifreq, 2) = valeur
                        resufr(ifreq, 3) = 0.0d0
                        resufr(ifreq, 14) = preci
                    end do
                end if
            else
                om0 = valp(jfreq)
                d1 = det(jfreq)*10.d0**(idet(jfreq)-idet(jfreq+1))
                d2 = det(jfreq+1)
                f1 = valp(jfreq)
                f2 = valp(jfreq+1)
                ip1 = 1
                ip2 = 1
                w1 = 1.d0
                w2 = 1.d0
                ier = 0
!
!              --- CALCUL DU ZERO ---
                do niter = 1, nitf
                    itrig = 0
                    om = (w1*f2*d1-w2*f1*d2)/(w1*d1-w2*d2)
                    if (niter .le. 3) then
!                    --- CONTROLE DE NON SORTIE DES BORNES ---
                        ftol = (f2-f1)/20.d0
                        if (om .lt. f1+ftol) om = f1+ftol
                        if (om .gt. f2-ftol) om = f2-ftol
                    end if
25                  continue
                    idet0 = 0
! --- POUR OPTIMISER ON NE GARDE PAS LA FACTO (SI MUMPS)
                    call vpstur(lraide, om, lmasse, ldynam, det0, &
                                idet0, ieme0, ier, solveu, .true._1, &
                                .false._1)
                    preci = abs(om0-om)/om0
                    if (preci .le. tolf .or. ier .ne. 0) then
                        nbiter = niter
                        goto 60
                    end if
                    om0 = om
                    if ((det0*d1) .ge. 0.d0) then
                        dd1 = det0*10.d0**(idet0-idet(jfreq+1))
                        if (abs(dd1) .gt. abs(d1) .and. itrig .eq. 0) then
                            om = (f1+f2)*0.5d0
                            itrig = 1
                            goto 25
                        end if
                        d1 = dd1
                        f1 = om
                        ip1 = ip1+1
                        ip2 = 1
                    else
                        dd2 = det0*10.d0**(idet0-idet(jfreq+1))
                        if (abs(dd2) .gt. abs(d2) .and. itrig .eq. 0) then
                            om = (f1+f2)*0.5d0
                            itrig = 1
                            goto 25
                        end if
                        d2 = dd2
                        f2 = om
                        ip2 = ip2+1
                        ip1 = 1
                    end if
                    w1 = 2**(((ip1-1)*(ip1-2))/2)
                    w2 = 2**(((ip2-1)*(ip2-2))/2)
                end do
                nbiter = -nitf
!
!              --- SORTIE DE LA BOUCLE SUR LES FREQUENCES ---
60              continue
                if (ier .ne. 0) then
                    valeur = 1.01d0*om
                else
                    valeur = (w1*f2*d1-w2*f1*d2)/(w1*d1-w2*d2)
                end if
            end if
!
            ifreq = ifreq+1
            resufi(ifreq, 1) = ieme(jfreq+1)
            resufi(ifreq, 2) = max(npas(jfreq), npas(jfreq+1))
            resufi(ifreq, 3) = nbiter
            resufi(ifreq, 4) = nbiter
            resufr(ifreq, 2) = valeur
            resufr(ifreq, 3) = 0.0d0
            resufr(ifreq, 14) = preci
1           continue
        end do
!
    else
!
!        --- ON PREND LE MILIEU DE L'INTERVALLE ---
        nbiter = 0
        do jfreq = 1, nfreq-1
            if (abs(ieme(jfreq+1)-ieme(jfreq)) .gt. 0) then
                ifreq = ifreq+1
                valeur = (valp(jfreq)+valp(jfreq+1))*0.5d0
                preci = (valp(jfreq+1)-valp(jfreq))/valeur
                resufi(ifreq, 1) = ieme(jfreq+1)
                resufi(ifreq, 2) = max(npas(jfreq), npas(jfreq+1))
                resufi(ifreq, 3) = nbiter
                resufi(ifreq, 4) = nbiter
                resufr(ifreq, 2) = valeur
                resufr(ifreq, 3) = 0.0d0
                resufr(ifreq, 14) = preci
            end if
        end do
!
    end if
!
!     --- MISE EN EXTENSION DES VALEURS PROPRES "MULTIPLES" ---
    if (abs(resufi(ifreq, 1)-resufi(1, 1)) .gt. ifreq-1) then
        nfreq = 0
        do jfreq = 1, ifreq-1
            nfreq = nfreq+1
            ntrou = abs(resufi(nfreq+1, 1)-resufi(nfreq, 1))-1
            if (ntrou .gt. 0) then
                do idec = nfreq+ifreq-jfreq, nfreq+1, -1
                    do jdec = 1, 4
                        resufi(idec+ntrou, jdec) = resufi(idec, jdec)
                    end do
                    do jdec = 1, 14
                        resufr(idec+ntrou, jdec) = resufr(idec, jdec)
                    end do
                end do
                do idec = 2, ntrou
                    do jdec = 1, 4
                        resufi(nfreq+idec, jdec) = resufi(nfreq+1, jdec)
                    end do
                    do jdec = 1, 14
                        resufr(nfreq+idec, jdec) = resufr(nfreq+1, jdec)
                    end do
                end do
                nfreq = nfreq+ntrou
            end if
        end do
        nfreq = nfreq+1
    else
!        --- PAS DE VALEURS PROPRES "MULTIPLES" ---
        nfreq = ifreq
    end if
!
!
end subroutine
