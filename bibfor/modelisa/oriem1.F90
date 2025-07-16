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
subroutine oriem1(ma, kdim, numa2d, numa3d)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
#include "blas/ddot.h"
    character(len=8), intent(in) :: ma
    character(len=2), intent(in) :: kdim
    integer(kind=8), intent(in) :: numa2d
    integer(kind=8), intent(inout) :: numa3d
! But :  regarder si la maille numa3d est bien du cote "-" de la normale
!        de la maille de peau numa2d.
!        Si ce n'est pas le cas, met a zero numa3d
!
! Pour determiner la disposition relative des mailles :
!        On NE se sert PAS de l'orientation des faces des mailles
!        volumiques (qui pourraient etre "retournees").
!        On se sert de la position geometrique des noeuds.
!        Si les mailles volumiques sont degenerees (sans volume), l'algorithme
!        echoue et l'on met numa3d a zero.
! ======================================================================
! in   : kdim   : '3D' / '2D'
! in   : ma     : nom du maillage
! in   : numa2d : numero d'une maille de peau
! inout: numa3d : numero d'une maille "volumique"
! ======================================================================
!
    integer(kind=8) :: ino, n1, n2, n3, ic, nbno1, nbno3, indi
    real(kind=8) :: nor1(3), n1n2(3), n1n3(3), ps1
    integer(kind=8), pointer :: lino1(:) => null()
    integer(kind=8), pointer :: lino3(:) => null()
    real(kind=8), pointer :: coor(:) => null()
    character(len=24) :: valk(2)
    blas_int :: b_incx, b_incy, b_n
!
!
! ========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
    call jeveuo(ma//'.COORDO    .VALE', 'L', vr=coor)
!
    call jeveuo(jexnum(ma//'.CONNEX', numa2d), 'L', vi=lino3)
    call jelira(jexnum(ma//'.CONNEX', numa2d), 'LONMAX', nbno3)
!
    call jeveuo(jexnum(ma//'.CONNEX', numa3d), 'L', vi=lino1)
    call jelira(jexnum(ma//'.CONNEX', numa3d), 'LONMAX', nbno1)
!
!
!   -- 1. Calcul de la normale de la maille de peau: nor1
!   ------------------------------------------------------
    n1 = lino3(1)
    n2 = lino3(2)
!
!   -- cas 3D :
    if (kdim .eq. '3D') then
        n3 = lino3(3)
        do ic = 1, 3
            n1n2(ic) = coor(3*(n2-1)+ic)-coor(3*(n1-1)+ic)
            n1n3(ic) = coor(3*(n3-1)+ic)-coor(3*(n1-1)+ic)
        end do
        call provec(n1n2, n1n3, nor1)
!
!   -- cas 2D :
    else
        do ic = 1, 3
            n1n2(ic) = coor(3*(n2-1)+ic)-coor(3*(n1-1)+ic)
        end do
        ASSERT(n1n2(3) .eq. 0.d0)
!
!       -- on tourne n1n2 de +90 degres :
        nor1(1) = -n1n2(2)
        nor1(2) = n1n2(1)
        nor1(3) = 0.d0
    end if
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    ASSERT(ddot(b_n, nor1, b_incx, nor1, b_incy) .gt. 0)
!
!
!   -- 2. position de la maille numa3d par rapport a la maille de peau  => ps1
!   --------------------------------------------------------------------------
    do ino = 1, nbno1
!        -- on cherche un noeud de lino1 (n2) qui ne soit pas un noeud
!           de la peau :
        indi = indiis(lino3, lino1(ino), 1, nbno3)
        if (indi .eq. 0) then
            n2 = lino1(ino)
            n1 = lino3(1)
            do ic = 1, 3
                n1n2(ic) = coor(3*(n2-1)+ic)-coor(3*(n1-1)+ic)
            end do
!           -- ps1 > 0 <=> la normale de la peau est orientee comme la
!                         la normale exterieure de la maille 1
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            ps1 = ddot(b_n, n1n2, b_incx, nor1, b_incy)
            goto 40
!
        end if
    end do
    ASSERT(.false.)
40  continue
!
!
!   -- si numa3d est degeneree ou du cote "+", on la supprime :
    if (ps1 .ge. 0.d0) then
        valk(1) = int_to_char8(numa3d)
        valk(2) = int_to_char8(numa2d)
        call utmess('A', 'CALCULEL3_47', nk=2, valk=valk)
        numa3d = 0
    end if
!
    call jedema()
end subroutine
