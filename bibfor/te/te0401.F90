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
subroutine te0401(optioz, nomtz)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/bsthco.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matpgl.h"
#include "asterfort/r8inir.h"
#include "asterfort/tranlg.h"
#include "asterfort/utvtsv.h"
#include "asterfort/vdxrig.h"
    character(len=*) :: optioz, nomtz
    character(len=16) :: option, nomte
!     ----------------------------------------------------------------
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE : COQUE_3D
!     ----------------------------------------------------------------
!
    integer(kind=8) :: nb1, nb2, nddlet
    integer(kind=8) :: lzr
    integer(kind=8) :: jgeom, jener
    integer(kind=8) :: i, j, kompt
    integer(kind=8) :: iu, imatuu
    real(kind=8) :: matloc(51, 51), plg(9, 3, 3)
    real(kind=8) :: vrs(1326)
    real(kind=8) :: bsigth(51), enerth
    aster_logical :: indith
! DEB
!
    option = optioz
    nomte = nomtz
!
    enerth = 0.0d0
!
    call jevech('PGEOMER', 'L', jgeom)
!
    if (option .eq. 'RIGI_MECA') call jevech('PMATUUR', 'E', imatuu)
!
    if (option .eq. 'RIGI_MECA' .or. option .eq. 'EPOT_ELEM') then
!
        call vdxrig(nomte, zr(jgeom), matloc, nb1, 0, &
                    0)
!
!     CONSTRUCTION DE LA MATRICE DE PASSAGE REPERE GLOBAL REPERE LOCAL
!
        call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
        nb2 = nb1+1
        call matpgl(nb2, zr(lzr), plg)
!
        call r8inir(1326, 0.d0, vrs, 1)
!
        nddlet = 6*nb1+3
!
        call tranlg(nb1, 51, nddlet, plg, matloc, &
                    vrs)
!
    else
        ASSERT(.false.)
    end if
!
!
    if (option .eq. 'RIGI_MECA') then
!
!--------- STOCKAGE
!
!--------- COMPTEUR DE POSITION
!
        kompt = 0
!
        do j = 1, 6*nb1+3
            do i = 1, j
                kompt = kompt+1
                zr(imatuu-1+kompt) = vrs(kompt)
            end do
        end do
!
    end if
!
!---- ENERGIES DE DEFORMATION ELASTIQUE
!
    if (option .eq. 'EPOT_ELEM') then
!
!------- LECTURE DE L'ADRESSE
!
        call jevech('PENERDR', 'E', jener)
!
!------- ADRESSE DES DEPLACEMENTS
!
        call jevech('PDEPLAR', 'L', iu)
!
!
!------- ENERGIE DE DEFORMATION TOTALE
!
        call utvtsv('ZERO', 6*nb1+3, vrs, zr(iu), zr(jener))
!
        zr(jener) = 0.5d0*zr(jener)
!
        call bsthco(nomte, bsigth, indith)
!
        if (indith) then
            do i = 1, 6*nb1+3
                enerth = enerth+bsigth(i)*zr(iu+i-1)
            end do
            zr(jener) = zr(jener)-enerth
        end if
!
        if (abs(zr(jener)) .gt. 1.d-6) then
!
!--------- ENERGIE DE DEFORMATION DE MEMBRANE
!
            call vdxrig(nomte, zr(jgeom), matloc, nb1, 1, &
                        0)
!
            call r8inir(1326, 0.d0, vrs, 1)
!
            call tranlg(nb1, 51, nddlet, plg, matloc, &
                        vrs)
!
            call utvtsv('ZERO', 6*nb1+3, vrs, zr(iu), zr(jener+1))
!
            zr(jener+1) = 0.5d0*zr(jener+1)
!
!
!--------- ENERGIE DE DEFORMATION DE FLEXION
!
            call vdxrig(nomte, zr(jgeom), matloc, nb1, 0, &
                        1)
!
            call r8inir(1326, 0.d0, vrs, 1)
!
            call tranlg(nb1, 51, nddlet, plg, matloc, &
                        vrs)
!
            call utvtsv('ZERO', 6*nb1+3, vrs, zr(iu), zr(jener+2))
!
            zr(jener+2) = 0.5d0*zr(jener+2)
!
!--------- VALEURS RELATIVES
!
            zr(jener+1) = zr(jener+1)/zr(jener)
            zr(jener+2) = zr(jener+2)/zr(jener)
!
        else
!
            call r8inir(2, 0.d0, zr(jener+1), 1)
!
        end if
!
!
    end if
!
!
end subroutine
