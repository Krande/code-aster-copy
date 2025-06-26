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
subroutine oriem0(kdim, type, coor, lino1, nbno1, &
                  lino2, nbno2, lino3, nbno3, ipos, &
                  indmai)
    implicit none
! person_in_charge: jacques.pellet at edf.fr
#include "asterfort/assert.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/provec.h"
#include "blas/ddot.h"
    character(len=2), intent(in) :: kdim
    character(len=8), intent(in) :: type
    real(kind=8), intent(in) :: coor(*)
    integer(kind=8), intent(in) :: lino1(*)
    integer(kind=8), intent(in) :: nbno1
    integer(kind=8), intent(in) :: lino2(*)
    integer(kind=8), intent(in) :: nbno2
    integer(kind=8), intent(in) :: lino3(*)
    integer(kind=8), intent(in) :: nbno3
    integer(kind=8), intent(out) :: ipos
    integer(kind=8), intent(out) :: indmai
! But :  Determiner la position relative de 2 mailles "volumiques"
!        (lino1 et lino2) par rapport a une maille de "peau" (lino3)
!        On suppose que les noeuds de la maille de peau appartiennent
!        aux 2 mailles "volumiques".
!
! Pour determiner la disposition relative des mailles :
!        On NE se sert PAS de l'orientation des faces des mailles
!        volumiques (qui pourraient etre "retournees").
!        On se sert de la postion geometrique des noeuds.
!        Si les mailles volumiques sont degenerees (sans volume), l'algorithme
!        echoue et l'on retourne indmai=-1/-2
! ======================================================================
! in  : kdim   : '3D' / '2D'
! in  : type   : type de la maille de peau (tria3, quad4, seg2, ...)
! in  : lino1  : liste des noeuds de la maille 1 (3d ou 2d)
! in  : nbno1  : nb de noeuds de lino1
! in  : lino2  : liste des noeuds de la maille 2 (3d ou 2d)
! in  : nbno2  : nb de noeuds de lino2
! in  : lino3  : liste des noeuds de la maille de peau (2.5d ou 1.5d)
! in  : nbno3  : nb de noeuds de de la maille de peau
! out : ipos   : = 0 : les 2 mailles 1 et 2 sont du meme cote
!                = 1  sinon
! out : indmai :
!          si ipos = 0 :  indmai = 0
!          si ipos = 1 :
!            / indmai = /1  /2
!              indmai est le numero de la maille qui est situee cote "-" de la
!              normale de la maille de peau.
!            / indmai = /-1  /-2  /-12 : on ne peut pas repondre car la maille 1
!                    ou la maille 2, ou les deux mailles 1 et 2 sont degenerees.
! ======================================================================
!
    integer(kind=8) :: ino, n1, n2, n3, ic, indi
    real(kind=8) :: nor1(3), n1n2(3), n1n3(3), ps1, ps2
    blas_int :: b_incx, b_incy, b_n
!
! ========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
    if (kdim .eq. '3D') then
        ASSERT(type(1:4) .eq. 'TRIA' .or. type(1:4) .eq. 'QUAD')
    else if (kdim .eq. '2D') then
        ASSERT(type(1:3) .eq. 'SEG')
    else
        ASSERT(.false.)
    end if
!
!
!   --- verification de la position des mailles 1 et 2
!       par rapport a la maille de peau :
!   ---------------------------------------------------
!
!   -- 1. Calcul de la normale de la maille de peau: nor1
    n1 = lino3(1)
    n2 = lino3(2)
!
!   -- cas 3D :
    if (type(1:3) .ne. 'SEG') then
        ASSERT(kdim .eq. '3D')
        n3 = lino3(3)
        do ic = 1, 3
            n1n2(ic) = coor(3*(n2-1)+ic)-coor(3*(n1-1)+ic)
            n1n3(ic) = coor(3*(n3-1)+ic)-coor(3*(n1-1)+ic)
        end do
        call provec(n1n2, n1n3, nor1)
!
!   -- cas 2D :
    else
        ASSERT(kdim .eq. '2D')
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
!
!   -- position de la maille 1 par rapport a la maille de peau  => ps1
    do ino = 1, nbno1
!        -- on cherche un noeud de lino1 (n2) qui ne soit pas un noeud
!           de la peau :
        indi = indiis(lino3(1), lino1(ino), 1, nbno3)
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
!   -- position de la maille 2 par rapport a la maille de peau  => ps2
    do ino = 1, nbno2
        indi = indiis(lino3(1), lino2(ino), 1, nbno3)
        if (indi .eq. 0) then
            n2 = lino2(ino)
            n1 = lino3(1)
            do ic = 1, 3
                n1n2(ic) = coor(3*(n2-1)+ic)-coor(3*(n1-1)+ic)
            end do
            b_n = to_blas_int(3)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            ps2 = ddot(b_n, n1n2, b_incx, nor1, b_incy)
            goto 70
!
        end if
    end do
    ASSERT(.false.)
70  continue
!
!
!   -- les mailles 1 et 2 sont elles du meme cote par rapport
!      a la maille de peau ?
!   ---------------------------------------------------------
!
    ipos = 0
    indmai = 0
!
!   -- si l'un des mailles est degeneree, on ne sait pas repondre :
    if (ps1 .eq. 0.d0) then
        ipos = 1
        indmai = -1
        if (ps2 .eq. 0.d0) then
            indmai = -12
        end if
    else
        if (ps2 .eq. 0.d0) then
            ipos = 1
            indmai = -2
        end if
    end if
!
    if (ipos .eq. 1) then
        ASSERT(indmai .lt. 0)
        goto 999
    end if
!
    if ((ps1 .gt. 0 .and. ps2 .gt. 0) .or. (ps1 .lt. 0 .and. ps2 .lt. 0)) then
        ipos = 0
        indmai = 0
    else
        ipos = 1
        if (ps1 .lt. 0) indmai = 1
        if (ps2 .lt. 0) indmai = 2
    end if
!
999 continue
    call jedema()
end subroutine
