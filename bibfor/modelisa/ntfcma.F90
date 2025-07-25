! --------------------------------------------------------------------
! Copyright (C) 2019 Christophe Durand - www.code-aster.org
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

subroutine ntfcma(compo, jmat, aniso, ifon)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: imate, jmat, ifon(6)
    character(len=*) :: compo
    aster_logical :: aniso
! ----------------------------------------------------------------------
!     OBTENTION DES ADRESSES DES FONCTIONS BETA ET LAMBDA DANS LE
!     MATERIAU CODE IMATE
! IN  COMPO  : NOM DU COMPORTEMENT CHERCHE
! IN  IMATE  : ADRESSE DU MATERIAU CODE
! IN  ANISO  : LOGICAL ANISOTROPE OU ISOTROPE
! OUT IFON   : ADRESSE RELATIVE DES PARAMETRES BETA ET LAMBDA
!      IFON(1) : ADRESSE RELATIVE DU PARAMETRE BETA OU -1 SI BETA ABSENT
!      IFON(2) : ADRESSE RELATIVE DU PARAMETRE LAMBDA   (ISOTROPIE)
!      IFON(3) : ADRESSE DU NOM DE LA FONCTION AFFINITE SI THER_HYDR
!      IFON(4) : ADRESSE RELATIVE DU PARAMETRE LAMBDA_L (ORTHOTROPIE)
!      IFON(5) : ADRESSE RELATIVE DU PARAMETRE LAMBDA_T (ORTHOTROPIE)
!      IFON(6) : ADRESSE RELATIVE DU PARAMETRE LAMBDA_N (ORTHOTROPIE)
!
    integer(kind=8) :: ipi, k, nbmat
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: idf, lfct, lmat
    character(len=16) :: valk(2), compo2
!-----------------------------------------------------------------------
    parameter(lmat=9, lfct=10)
! DEB ------------------------------------------------------------------
!
!
    nbmat = zi(jmat)
    ASSERT(nbmat .eq. 1)
    imate = jmat+zi(jmat+nbmat+1)
!
    ASSERT(compo(1:7) .eq. 'THER_NL' .or. compo(1:9) .eq. 'THER_HYDR' .or. compo .eq. ' ')
!
    if (compo .eq. ' ') then
        do k = 1, zi(imate+1)
            if ('THER_NL' .eq. zk32(zi(imate)+k-1) (1:7)) then
                ipi = zi(imate+2+k-1)
                compo2 = 'THER_NL'
                goto 11
            end if
        end do
        do k = 1, zi(imate+1)
            if ('THER_HYDR' .eq. zk32(zi(imate)+k-1) (1:9)) then
                ipi = zi(imate+2+k-1)
                compo2 = 'THER_HYDR'
                goto 11
            end if
        end do
    else
        do k = 1, zi(imate+1)
            if (compo(1:7) .eq. zk32(zi(imate)+k-1) (1:7)) then
                ipi = zi(imate+2+k-1)
                compo2 = compo
                goto 11
            end if
        end do
    end if
    do k = 1, zi(imate+1)
        if ('THER_' .eq. zk32(zi(imate)+k-1) (1:5)) then
            valk(1) = zk32(zi(imate)+k-1) (1:16)
            valk(2) = compo
            if (compo .eq. ' ') then
                call utmess('F', 'ELEMENTS2_65', sk=valk(1))
            else
                call utmess('F', 'ELEMENTS2_64', nk=2, valk=valk)
            end if
        end if
    end do
    if (compo .eq. ' ') then
        call utmess('F', 'ELEMENTS2_66')
    else
        call utmess('F', 'ELEMENTS2_63', sk=compo)
    end if
11  continue
!!!
! RECUPERATION DE L ADRESSE DE BETA DANS IFON(1)
!
    idf = zi(ipi)+zi(ipi+1)
    do k = 1, zi(ipi+2)
        if ('BETA    ' .eq. zk16(zi(ipi+3)+idf+k-1)) then
            ifon(1) = ipi+lmat-1+lfct*(k-1)
            goto 25
        end if
    end do
    call utmess('F', 'MODELISA5_44')
25  continue
!!!
! RECUPERATION DE L ADRESSE DE LA CONDUCTIVITE LAMBDA
! CAS ISOTROPE   : LAMBDA DANS IFON(2)
! CAS ORTHOTROPE : LAMBDA_L DANS IFON(4)
! CAS ORTHOTROPE : LAMBDA_T DANS IFON(5)
! CAS ORTHOTROPE : LAMBDA_N DANS IFON(6)
!
    if (.not. aniso) then
        do k = 1, zi(ipi+2)
            if ('LAMBDA  ' .eq. zk16(zi(ipi+3)+idf+k-1)) then
                ifon(2) = ipi+lmat-1+lfct*(k-1)
                goto 35
            end if
        end do
        call utmess('F', 'MODELISA5_45')
35      continue
    else
        do k = 1, zi(ipi+2)
            if ('LAMBDA_L' .eq. zk16(zi(ipi+3)+idf+k-1)) then
                ifon(4) = ipi+lmat-1+lfct*(k-1)
                goto 45
            end if
        end do
        call utmess('F', 'MODELISA5_46', sk='LAMBDA_L')
45      continue
        do k = 1, zi(ipi+2)
            if ('LAMBDA_T' .eq. zk16(zi(ipi+3)+idf+k-1)) then
                ifon(5) = ipi+lmat-1+lfct*(k-1)
                goto 55
            end if
        end do
        call utmess('F', 'MODELISA5_46', sk='LAMBDA_T')
55      continue
        do k = 1, zi(ipi+2)
            if ('LAMBDA_N' .eq. zk16(zi(ipi+3)+idf+k-1)) then
                ifon(6) = ipi+lmat-1+lfct*(k-1)
                goto 65
            end if
        end do
        call utmess('F', 'MODELISA5_46', sk='LAMBDA_N')
65      continue
    end if
!!!
! TRAITEMENT DE L HYDRATATION
! RECUPERATION DE L ADRESSE DE AFFINITE DANS IFON(3)
!
    if (compo2(1:9) .eq. 'THER_HYDR') then
        do k = 1, zi(ipi+2)
            if ('AFFINITE  ' .eq. zk16(zi(ipi+3)+idf+k-1)) then
                ifon(3) = ipi+lmat-1+lfct*(k-1)
                goto 75
            end if
        end do
        call utmess('F', 'MODELISA5_47')
75      continue
    end if
!
! FIN ------------------------------------------------------------------
end subroutine
