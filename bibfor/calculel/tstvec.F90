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
subroutine tstvec(perm, iad, nlong, type, sommi, &
                  sommr, nbign)
! aslint: disable=C1002,W0405
    implicit none
!
! BUT : RECUPERER 2 NOMBRES RESUMANT UN VECTEUR JEVEUX
!
! IN: PERM  K3 : /OUI/NON
!           NON : ON FAIT LA SOMME BETE DES ELEMENTS DU VECTEUR
!                 => UNE PERMUTATION DU VECTEUR NE SE VOIT PAS !
!           OUI : ON FAIT UNE "SOMME" QUI DONNE UN RESULTAT
!                 DEPENDANT UN PEU DE L'ORDRE DES ELEMENTS DU VECTEUR
! IN: IAD   I  : ADRESSE DU VECTEUR
! IN: NLONG  I  : LONGUEUR DU VECTEUR
! IN: TYPE  K3 : TYPE DES ELEMENTS DU VECTEUR :
!                   I/L/R/C/K8/K16/K24/K32/K80
!
! OUT: SOMMI   I      : SOMME(V(I)) QUELQUE SOIT LE TYPE DE V
! OUT: SOMMR   R      : SOMME(V(I)) SI V EST DE TYPE "R/C"
! OUT: NBIGN   I      : NOMBRE DE VALEURS IGNOREES DANS SOMMR :
!                       (UNDEF OU TRES GRAND )
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/tstk2i.h"
!
    character(len=3) :: type
    real(kind=8) :: sommr
    integer(kind=8) :: sommi, nbign
    character(len=24) :: k24
    aster_logical :: l
    character(len=*) :: perm
    character(len=8) :: k8
    real(kind=8) :: x
    integer(kind=8) :: iad, ix, c1, nlong, ico, k
    integer(kind=8) :: sommi2, i8
!
    equivalence(x, ix)
    equivalence(l, ix)
    equivalence(k8, ix)
!
!
    if (perm .eq. 'NON') then
        c1 = 0
    else
        c1 = 1
    end if
!
!
!     -- CALCUL DE SOMMR :
!     --------------------
    sommr = 0.d0
    ico = 0
    nbign = 0
    if (type .eq. 'R') then
        do k = 1, nlong
            x = zr(iad-1+k)
            if (.not. isnan(x)) then
                if (abs(x) .lt. 1.d300) then
                    ico = ico+1
                    sommr = sommr+(c1*mod(k, 3)+1)*x
                end if
            end if
        end do
        nbign = nlong-ico
    end if
    if (type .eq. 'C') then
        do k = 1, nlong
            x = dble(zc(iad-1+k))
            if (.not. isnan(x)) then
                if (abs(x) .lt. 1.d300) then
                    ico = ico+1
                    sommr = sommr+(c1*mod(k, 3)+1)*x
                end if
            end if
            x = dimag(zc(iad-1+k))
            if (.not. isnan(x)) then
                if (abs(x) .lt. 1.d300) then
                    ico = ico+1
                    sommr = sommr+(c1*mod(k, 3)+1)*x
                end if
            end if
        end do
        nbign = 2*nlong-ico
    end if
!
!
!     -- CALCUL DE SOMMI :
!     --------------------
    sommi2 = 0
    if (type .eq. 'I') then
        do k = 1, nlong
            i8 = zi(iad-1+k)
            sommi2 = sommi2+(c1*mod(k, 3)+1)*i8
        end do
    else if (type .eq. 'S') then
        do k = 1, nlong
            i8 = zi4(iad-1+k)
            sommi2 = sommi2+(c1*mod(k, 3)+1)*i8
        end do
    else if (type .eq. 'L') then
        do k = 1, nlong
            l = zl(iad-1+k)
            if (l) sommi2 = sommi2+(c1*mod(k, 3)+1)
        end do
    else if (type .eq. 'R') then
        sommi2 = 0
    else if (type .eq. 'C') then
        sommi2 = 0
    else if (type .eq. 'K8') then
        do k = 1, nlong
            sommi2 = sommi2+(c1*mod(k, 3)+1)*tstk2i(8, zk8(iad-1+k))
        end do
    else if (type .eq. 'K16') then
        do k = 1, nlong
            sommi2 = sommi2+(c1*mod(k, 3)+1)*tstk2i(16, zk16(iad-1+k))
        end do
    else if (type .eq. 'K24') then
        do k = 1, nlong
            sommi2 = sommi2+(c1*mod(k, 3)+1)*tstk2i(24, zk24(iad-1+k))
        end do
    else if (type .eq. 'K32') then
        do k = 1, nlong
            sommi2 = sommi2+(c1*mod(k, 3)+1)*tstk2i(32, zk32(iad-1+k))
        end do
    else if (type .eq. 'K80') then
        do k = 1, nlong
            sommi2 = sommi2+(c1*mod(k, 3)+1)*tstk2i(80, zk80(iad-1+k))
        end do
    end if
!
!     -- ON TRONQUE SOMMI2 (9 DERNIERS CHIFFRES) POUR AVOIR
!        LE MEME RESULTAT SUR LES PLATEFORMES I4 ET I8 :
    write (k24, '(I24)') sommi2
    read (k24(16:24), '(I9)') sommi
!
!
end subroutine
