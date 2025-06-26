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

subroutine te0142(option, nomte)
!
    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rccoma.h"

!
    character(len=16) :: option, nomte
!
    integer(kind=8), parameter :: sdim = 3
    integer(kind=8), parameter :: ncmp_el = 4
    integer(kind=8), parameter :: ncmp_th = 2
    integer(kind=8), parameter :: ncmptot = sdim+ncmp_el+ncmp_th

    integer(kind=8) :: igeom, i, j, nbpar, ipar
    integer(kind=8) :: ndim, npg1, kpg, spt, iret
    integer(kind=8) :: imate, ival, ivf, idecpg, idecno
    integer(kind=8) :: mater, nnos, nno, indir(sdim)
    real(kind=8) :: valres_el(ncmp_el), valres_th(ncmp_th)
    integer(kind=8) :: icodre_el(ncmp_el), icodre_th(ncmp_th), ndim2
    character(len=8) :: fami, poum, novrc
    character(len=16) :: nomres_el(ncmp_el), nomres_th(ncmp_th)
    character(len=8) :: nompar(sdim), nompar0(sdim)
    real(kind=8) :: xyzgau(sdim), xyzgau0(sdim)
    character(len=32) :: phenom_el, phenom_th, valk(2)

    aster_logical::lfound
!     ------------------------------------------------------------------
!
    if (option .eq. 'MATE_ELGA') then
        fami = 'MTGA'
    elseif (option .eq. 'MATE_ELEM') then
        fami = 'FPG1'
    else
        ASSERT(.false.)
    end if
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, jvf=ivf)
!
    if (lteatt('ABSO', 'OUI')) then
        ndim2 = ndim+1
    else
        ndim2 = ndim
    end if

    call jevech('PMATERC', 'L', imate)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERR', 'E', ival)
!
    mater = zi(imate)
    nomres_el(1) = 'E'
    nomres_el(2) = 'NU'
    nomres_el(3) = 'RHO'
    nomres_el(4) = 'ALPHA'
    valres_el(:) = 0.d0

    nomres_th(1) = 'LAMBDA'
    nomres_th(2) = 'RHO_CP'
    valres_th(:) = 0.d0

    spt = 1
    poum = '+'
!
    nompar0(1) = 'X'
    nompar0(2) = 'Y'
    nompar0(3) = 'Z'
!

    call rccoma(mater, 'ELAS', 0, phenom_el, iret)
    ASSERT((iret .eq. 0) .or. (iret .eq. 1))
    valk(1) = 'ELAS'
    valk(2) = ' '
    if (ALL(valk .ne. phenom_el)) then
        phenom_th = 'ELAS'
    end if

    call rccoma(mater, 'THER', 0, phenom_th, iret)
    ASSERT((iret .eq. 0) .or. (iret .eq. 1))
    valk(1) = 'THER'
    valk(2) = 'THER_NL'
    if (ALL(valk .ne. phenom_th)) then
        phenom_th = 'THER'
    end if

!   check if X, Y or Z are present in the command variables and store indirection in indir
    nbpar = 0
    do i = 1, ndim2
        lfound = .false.
        do ipar = 1, ca_nbcvrc_
            novrc = zk8(ca_jvcnom_-1+ipar)
            if (novrc .eq. nompar0(i)) then
                lfound = .true.
                cycle
            end if
        end do
        if (.not. lfound) then
            nbpar = nbpar+1
            indir(nbpar) = i
        end if
    end do
!   only use parameters in indir
    do i = 1, nbpar
        nompar(i) = nompar0(indir(i))
    end do

    do kpg = 1, npg1
        idecpg = nno*(kpg-1)-1
        ! ----- Coordinates for current Gauss point
        xyzgau0(:) = 0.d0
        do i = 1, nno
            idecno = ndim2*(i-1)-1
            do j = 1, ndim2
                xyzgau0(j) = xyzgau0(j)+zr(ivf+i+idecpg)*zr(igeom+j+idecno)
            end do
        end do
!           only use parameters in indir
        do i = 1, nbpar
            xyzgau(i) = xyzgau0(indir(i))
        end do

        do i = 1, sdim
            zr(ival-1+(kpg-1)*ncmptot+i) = xyzgau0(i)
        end do

        call rcvalb(fami, kpg, spt, poum, mater, &
                    ' ', phenom_el, nbpar, nompar, xyzgau, &
                    ncmp_el, nomres_el, valres_el, icodre_el, 0, 'NON')

        do i = 1, ncmp_el
            zr(ival-1+(kpg-1)*ncmptot+sdim+i) = valres_el(i)
        end do

        call rcvalb(fami, kpg, spt, poum, mater, &
                    ' ', phenom_th, nbpar, nompar, xyzgau, &
                    ncmp_th, nomres_th, valres_th, icodre_th, 0, 'NON')

        do i = 1, ncmp_th
            zr(ival-1+(kpg-1)*ncmptot+sdim+ncmp_el+i) = valres_th(i)
        end do

    end do
!
end subroutine
