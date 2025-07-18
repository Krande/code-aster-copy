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
subroutine rsindi(tysca, iaobj, paobj, jordr, ival, &
                  rval, kval, cval, epsi, crit, &
                  nbordr, nbtrou, nutrou, ndim)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
!
    integer(kind=8) :: nbordr, nbtrou, nutrou(*), ndim, ival, paobj
    real(kind=8) :: rval, epsi
    character(len=4) :: tysca
    character(len=*) :: kval, crit
    complex(kind=8) :: cval
!
!      TROUVER DANS LA TABLE ZX(IAOBJ-1+I),I=1,NBORDR LE SCALAIRE
!      IVAL,RVAL,CVAL...
!               AVEC LA PRECISION  / RELATIVE EPSI
!                                  / ABSOLUE  EPSI
!      LE TEST FAIT EST : ABS(V-VR).LE.ABS(EPSI*VR) EN RELATIF
!                    OU   ABS(V-VR).LE.ABS(EPSI)    EN ABSOLU
! ----------------------------------------------------------------------
! IN  : TYSCA  : R8 OU C16 OU I8 OU K8,  K16, K24, K32, K80
! IN  : IAOBJ  : ADRESSE DE LA TABLE DANS ZI, ZR OU ZC
! IN  : PAOBJ  : "PAS" dU parcours de l'objet :
!                 ZX(IAOBJ),ZX(IAOBJ+1*PAOBJ),,ZX(IAOBJ+2*PAOBJ), ...
! IN  : JORDR  : ADRESSE DU .ORDR DU RESULTAT
! IN  : IVAL   : VALEUR CHECHEE SI ENTIERE.
! IN  : RVAL   : VALEUR CHECHEE SI REELLE.
! IN  : KVAL   : VALEUR CHECHEE SI CARACTERE.
! IN  : CVAL   : VALEUR CHECHEE SI COMPLEXE.
! IN  : CRIT   : CRITERE : 'RELATIF' OU 'ABSOLU'
! IN  : EPSI   : PRECISION VOULUE.
! IN  : NBORDR : DIMENSION MAXI DE LA TABLE DE RECHERCHE.
! IN  : NDIM   : DIMENSION MAXI DE LA LISTE A REMPLIR.
! OUT : NBTROU : NOMBRE DE VALEURS CONVENABLES.
! OUT : NUTROU : LISTE DES INDICES DES VALEURS CONVENABLES.
!                   SI NBTROU EST > NDIM , ON SIGNALE L'ERREUR EN
!                   RENDANT NBTROU = - NBTROU
! ----------------------------------------------------------------------
    character(len=8) :: crit2
    aster_logical :: depass, trouve
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iaobj, jordr
!-----------------------------------------------------------------------
    crit2 = crit
    nbtrou = 0
    depass = .false.
    if (tysca(1:1) .eq. 'R') then
        do i = 1, nbordr
            if (crit2(1:4) .eq. 'RELA') then
                if (abs(zr(iaobj+(i-1)*paobj)-rval) .le. abs(epsi*rval)) then
                    trouve = .true.
                else
                    trouve = .false.
                end if
            else if (crit2(1:4) .eq. 'ABSO') then
                if (abs(zr(iaobj+(i-1)*paobj)-rval) .le. abs(epsi)) then
                    trouve = .true.
                else
                    trouve = .false.
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
            if (trouve) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca(1:1) .eq. 'I') then
        do i = 1, nbordr
            if (zi(iaobj+(i-1)*paobj) .eq. ival) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca .eq. 'K8  ') then
        do i = 1, nbordr
            if (zk8(iaobj+(i-1)*paobj) .eq. kval) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca .eq. 'K16 ') then
        do i = 1, nbordr
            if (zk16(iaobj+(i-1)*paobj) .eq. kval) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca .eq. 'K24 ') then
        do i = 1, nbordr
            if (zk24(iaobj+(i-1)*paobj) .eq. kval) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca .eq. 'K32 ') then
        do i = 1, nbordr
            if (zk32(iaobj+(i-1)*paobj) .eq. kval) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca .eq. 'K80 ') then
        do i = 1, nbordr
            if (zk80(iaobj+(i-1)*paobj) .eq. kval) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else if (tysca(1:1) .eq. 'C') then
        do i = 1, nbordr
            if (crit2(1:4) .eq. 'RELA') then
                if (abs(zc(iaobj+(i-1)*paobj)-cval) .le. abs(epsi*cval)) then
                    trouve = .true.
                else
                    trouve = .false.
                end if
            else if (crit2(1:4) .eq. 'ABSO') then
                if (abs(zc(iaobj+(i-1)*paobj)-cval) .le. abs(epsi)) then
                    trouve = .true.
                else
                    trouve = .false.
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
            if (trouve) then
                nbtrou = nbtrou+1
                if (nbtrou .le. ndim) then
                    nutrou(nbtrou) = zi(jordr+i-1)
                else
                    depass = .true.
                end if
            end if
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
!
    if (depass) nbtrou = -nbtrou
!
end subroutine
