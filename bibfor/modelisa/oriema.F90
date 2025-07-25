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
subroutine oriema(nomail, tpmail, nbnmai, lnmail, typ3d, &
                  lnm3d, ndim, coor, reorie, norien, &
                  ifm, niv)
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
    integer(kind=8) :: nbnmai, lnmail(*), lnm3d(*), ndim, norien, ifm, niv
    real(kind=8) :: coor(*)
    aster_logical :: reorie
    character(len=8) :: nomail, tpmail, typ3d
    character(len=24) :: valk(2)
!.======================================================================
! ORIEMA  --  ORIENTATION DE LA MAILLE DE NOM NOMAIL
!             DE TELLE MANIERE A CE QUE LA NORMALE A CETTE MAILLE
!             SOIT EXTERIEURE AU VOLUME.
!             CETTE MAILLE EST UNE MAILLE DE PEAU :
!               .EN 2D C'EST UNE MAILLE DE TYPE SEG2 OU SEG3
!               .EN 3D C'EST UNE MAILLE DE TYPE TRIA3, TRIA6,
!                                               QUAD4, QUAD8 OU QUAD9
!
! IN  : REORIE : = .FALSE.  ON VERIFIE L'ORIENTATION
!                = .TRUE.   ON REORIENTE
! OUT : NORIEN : = 0  LA MAILLE EST BIEN ORIENTEE
!                = 1  LA MAILLE EST A REORIENTER
!
!.========================= DEBUT DES DECLARATIONS ====================
!
    integer(kind=8), parameter :: nbnomx = 27
    integer(kind=8) :: i, ic, n1, n2, n3, nbnsm, nbns3d, ino, nbnoe
    integer(kind=8) :: lisnoe(nbnomx)
    real(kind=8) :: coon1(3), coon2(3), coon3(3), n1n2(3), n1n3(3)
    real(kind=8) :: n(3), norme, n1g(3), xg3d(3), xgm(3), xgn, zero
    blas_int :: b_incx, b_incy, b_n
!
! ========================= DEBUT DU CODE EXECUTABLE ==================
!
    norien = 0
    zero = 0.d0
!
    do i = 1, nbnmai
        lisnoe(i) = lnmail(i)
    end do
    do i = 1, 3
        n(i) = zero
        n1n2(i) = zero
        n1n3(i) = zero
        xgm(i) = zero
        xg3d(i) = zero
    end do
!
! --- NUMERO DES 2 (3 EN 3D) PREMIERS NOEUDS DE LA MAILLE :
!     ---------------------------------------------------
    n1 = lisnoe(1)
    n2 = lisnoe(2)
    if (ndim .eq. 3) n3 = lisnoe(3)
!
    do ic = 1, 3
        coon1(ic) = coor(3*(n1-1)+ic)
        coon2(ic) = coor(3*(n2-1)+ic)
        n1n2(ic) = coon2(ic)-coon1(ic)
    end do
!
    if (ndim .eq. 2) then
        n(1) = n1n2(2)
        n(2) = -n1n2(1)
    else if (ndim .eq. 3) then
        do ic = 1, 3
            coon3(ic) = coor(3*(n3-1)+ic)
            n1n3(ic) = coon3(ic)-coon1(ic)
        end do
        call provec(n1n2, n1n3, n)
    end if
    call normev(n, norme)
!
    if (tpmail(1:4) .eq. 'QUAD') then
        nbnsm = 4
    else if (tpmail(1:4) .eq. 'TRIA') then
        nbnsm = 3
    else if (tpmail(1:3) .eq. 'SEG') then
        if (ndim .eq. 3) goto 99
        nbnsm = 2
    else
        valk(1) = nomail
        valk(2) = tpmail
        call utmess('F', 'MODELISA5_94', nk=2, valk=valk)
    end if
!
! --- CENTRE DE GRAVITE DE LA MAILLE DE PEAU
!
    do i = 1, nbnsm
        ino = lisnoe(i)
        xgm(1) = xgm(1)+coor(3*(ino-1)+1)
        xgm(2) = xgm(2)+coor(3*(ino-1)+2)
        xgm(3) = xgm(3)+coor(3*(ino-1)+3)
    end do
    xgm(1) = xgm(1)/nbnsm
    xgm(2) = xgm(2)/nbnsm
    xgm(3) = xgm(3)/nbnsm
!
    if (typ3d(1:4) .eq. 'HEXA') then
        nbns3d = 8
    else if (typ3d(1:4) .eq. 'PENT') then
        nbns3d = 6
    else if (typ3d(1:4) .eq. 'PYRA') then
        nbns3d = 5
    else if (typ3d(1:4) .eq. 'TETR') then
        nbns3d = 4
    else if (typ3d(1:4) .eq. 'QUAD') then
        nbns3d = 4
    else if (typ3d(1:4) .eq. 'TRIA') then
        nbns3d = 3
    else
        ASSERT(ASTER_FALSE)
    end if
!
! --- DETERMINATION DU CENTRE DE GRAVITE DE LA MAILLE 3D
!
    do i = 1, nbns3d
        ino = lnm3d(i)
        xg3d(1) = xg3d(1)+coor(3*(ino-1)+1)
        xg3d(2) = xg3d(2)+coor(3*(ino-1)+2)
        xg3d(3) = xg3d(3)+coor(3*(ino-1)+3)
    end do
    xg3d(1) = xg3d(1)/nbns3d
    xg3d(2) = xg3d(2)/nbns3d
    xg3d(3) = xg3d(3)/nbns3d
    n1g = xg3d-xgm
!
    call normev(n1g, norme)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    xgn = ddot(b_n, n1g, b_incx, n, b_incy)
!
! --- SI XGN > 0, LA NORMALE A LA MAILLE DE PEAU
! --- EST DIRIGEE VERS L'INTERIEUR DU VOLUME, IL FAUT
! --- LA REORIENTER :
!     -------------
    if (xgn .gt. zero) then
        norien = norien+1
        if (.not. reorie) goto 99
        if (tpmail(1:5) .eq. 'QUAD9' .or. tpmail(1:5) .eq. 'TRIA7') then
            nbnoe = nbnmai-1
        else
            nbnoe = nbnmai
        end if
        ino = 0
        do i = nbnsm, 1, -1
            ino = ino+1
            lnmail(i) = lisnoe(ino)
        end do
        if (nbnsm .ne. nbnoe) then
            ino = 0
            do i = nbnoe-1, nbnsm+1, -1
                ino = ino+1
                lnmail(i) = lisnoe(ino+nbnsm)
            end do
            lnmail(nbnoe) = lisnoe(nbnoe)
        end if
        if (niv .eq. 2) then
            write (ifm, *) ' REORIENTATION MAILLE ', nomail, ' NOEUDS ', ( &
                lnmail(i), i=1, nbnmai), ' PRODUIT SCALAIRE ', xgn
        end if
    end if
!
99  continue
!
end subroutine
