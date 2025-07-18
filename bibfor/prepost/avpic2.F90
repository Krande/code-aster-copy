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
subroutine avpic2(method, nbvec, nbordr, jrtrv, jitrv, &
                  npoin, jvalpo, jvalor, npic, jpic, &
                  jordpi)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbvec, nbordr, npoin(nbvec), jvalor
    integer(kind=8) :: npic(nbvec), jordpi
    integer(kind=8) :: jitrv
!
    character(len=8) :: method
    integer(kind=8) :: jrtrv, jvalpo, jpic
! ----------------------------------------------------------------------
! BUT: EXTRACTION DES PICS POUR RAINFLOW <=> REARANGEMENT DES PICS,
!      PIC LE PLUS GRAND AU DEBUT ET A LA FIN.
! ----------------------------------------------------------------------
! ARGUMENTS:
! METHOD    IN   K  : METHODE D'EXTRACTION DES PICS, PAR EXEMPLE :
!                     RAINFLOW.
! NBVEC     IN   I  : NOMBRE DE VECTEURS NORMAUX.
! NBORDR    IN   I  : NOMBRE DE NUMERO D'ORDRE.
! RTRV      IN   R  : VECTEUR DE TRAVAIL REEL (POUR LES POINTS)
! ITRV      IN   I  : VECTEUR DE TRAVAIL ENTIER (POUR LES NUME_ORDRE)
! NPOIN     IN   I  : NOMBRE DE PICS DETECTES POUR TOUS LES VECTEURS
!                     NORMAUX.
! VALPOI    IN   R  : VALEUR DES PICS DETECTES POUR TOUS LES VECTEURS
!                     NORMAUX.
! VALORD    IN   I  : NUMEROS D'ORDRE ASSOCIES AUX PICS DETECTES POUR
!                     TOUS LES VECTEURS NORMAUX.
! NPIC      OUT  I  : NOMBRE DE PICS DETECTES POUR TOUS LES VECTEURS
!                     NORMAUX APRES REARANGEMENT DES PICS.
! PIC       OUT  R  : VALEUR DES PICS DETECTES POUR TOUS LES VECTEURS
!                     NORMAUX APRES REARANGEMENT DES PICS.
! ORDPIC    OUT  I  : NUMEROS D'ORDRE ASSOCIES AUX PICS DETECTES POUR
!                     TOUS LES VECTEURS NORMAUX APRES REARANGEMENT
!                     DES PICS.
!
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: ivect, adrs, i, nmax, ntrv, ointer
    real(kind=8) :: epsilo, pmax, pinter, dp1, dp2
    character(len=8) :: k8b
!-----------------------------------------------------------------------
!234567                                                              012
!-----------------------------------------------------------------------
    epsilo = 1.0d-7
!-----------------------------------------------------------------------
!
    call jemarq()
!
    if (method .ne. 'RAINFLOW') then
        k8b = method(1:8)
        call utmess('F', 'PREPOST_4', sk=k8b)
    end if
!
    do ivect = 1, nbvec
!
! LE TEST SI (NPOIN(IVECT) .EQ. 0) EST EQUIVALENT
! AU TEST SI (IFLAG(IVECT) .EQ. 3).
        if (npoin(ivect) .eq. 0) then
            goto 10
        end if
!
! ----- RECHERCHE DU POINT LE PLUS GRAND EN VALEUR ABSOLUE -----
!
        ASSERT(nbordr .ge. npoin(ivect))
        adrs = (ivect-1)*(nbordr+2)
!
        pmax = abs(zr(jvalpo+(ivect-1)*nbordr+1))
        nmax = 1
        do i = 2, npoin(ivect)
            if (abs(zr(jvalpo+(ivect-1)*nbordr+i)) .gt. pmax*(1.0d0+epsilo)) then
!
                pmax = abs(zr(jvalpo+(ivect-1)*nbordr+i))
                nmax = i
            end if
        end do
        pmax = zr(jvalpo+(ivect-1)*nbordr+nmax)
        ntrv = npoin(ivect)
!
! ----- REARANGEMENT AVEC POINT LE PLUS GRAND AU DEBUT
!       ET A LA FIN                                    -----
!
        do i = nmax, npoin(ivect)
            zr(jrtrv+i-nmax+1) = zr(jvalpo+(ivect-1)*nbordr+i)
            zi(jitrv+i-nmax+1) = zi(jvalor+(ivect-1)*nbordr+i)
        end do
        do i = 1, nmax-1
            zr(jrtrv+npoin(ivect)+i-nmax+1) = zr(jvalpo+(ivect-1)*nbordr+i)
            zi(jitrv+npoin(ivect)+i-nmax+1) = zi(jvalor+(ivect-1)*nbordr+i)
        end do
!
! ----- EXTRACTION DES PICS SUR LE VECTEUR REARANGE -----
!
! 1. LE PREMIER POINT EST UN PIC
!
        npic(ivect) = 1
        zr(jpic+adrs+1) = zr(jrtrv+1)
        pinter = zr(jrtrv+2)
        zi(jordpi+adrs+1) = zi(jitrv+1)
        ointer = zi(jitrv+2)
!
! 2. RECHERCHE DE TOUS LES PICS
!
        do i = 3, ntrv
            dp1 = pinter-zr(jpic+adrs+npic(ivect))
            dp2 = zr(jrtrv+i)-pinter
!
! 2.1 ON CONSERVE LE POINT INTERMEDIAIRE COMME UN PIC
!
            if (dp1*dp2 .lt. -epsilo) then
                npic(ivect) = npic(ivect)+1
                zr(jpic+adrs+npic(ivect)) = pinter
                zi(jordpi+adrs+npic(ivect)) = ointer
            end if
!
! 2.2 LE DERNIER POINT DEVIENT POINT INTERMEDIAIRE
!
            pinter = zr(jrtrv+i)
            ointer = zi(jitrv+i)
        end do
!
! 3. LE DERNIER POINT EST UN PIC
!
        npic(ivect) = npic(ivect)+1
        zr(jpic+adrs+npic(ivect)) = zr(jrtrv+ntrv)
        zi(jordpi+adrs+npic(ivect)) = zi(jitrv+ntrv)
!
10      continue
    end do
!
    call jedema()
!
end subroutine
