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

subroutine te0066(option, nomte)
!
    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/matrot.h"
#include "asterfort/rcangm.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DE L'ENERGIE THERMIQUE A L'EQUILIOBRE
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'ETHE_ELEM'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!
    character(len=8) :: nompar(ca_nbcvrc_+1), fami, poum, novrc
    character(len=16) :: nomres(3)
    character(len=32) :: phenom
    integer(kind=8) :: icodre(3)
    real(kind=8) :: valpar(ca_nbcvrc_+1), lambda(1), poids, epot, valres(3), lambor(3)
    real(kind=8) :: dfdx(27), dfdy(27), dfdz(27), flux, fluy, fluz
    real(kind=8) :: angmas(3), point(3), fluglo(3), fluloc(3), p(3, 3)
    integer(kind=8) :: i, ipoids, ivf, idfde, igeom, imate, kpg, spt, ino
    integer(kind=8) :: ndim, jgano, nno, kp, npg1, iener, itemp, itempe, l, ipar
    aster_logical :: aniso
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpar, nnos
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PTEMPER', 'L', itempe)
    call jevech('PENERDR', 'E', iener)
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call tecach('ONO', 'PINSTR', 'L', iret, iad=itemp)
    if (itemp .eq. 0) then
        nbpar = 0
        nompar(1) = ' '
        valpar(1) = 0.d0
    else
        nbpar = 1
        nompar(1) = 'INST'
        valpar(1) = zr(itemp)
    end if
!
    do ipar = 1, ca_nbcvrc_
        novrc = zk8(ca_jvcnom_-1+ipar)
        nbpar = nbpar+1
        nompar(nbpar) = novrc
        call rcvarc(' ', nompar(nbpar), poum, fami, kpg, spt, valpar(nbpar), iret)
        ASSERT(iret .eq. 0)
    end do
!
    call rccoma(zi(imate), 'THER', 1, phenom, iret)
    if (phenom .eq. 'THER') then
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THER', nbpar, nompar, [valpar], &
                    1, 'LAMBDA', lambda, icodre, 1)
        aniso = .false.
    else if (phenom .eq. 'THER_ORTH') then
        nomres(1) = 'LAMBDA_L'
        nomres(2) = 'LAMBDA_T'
        nomres(3) = 'LAMBDA_N'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THER_ORTH', nbpar, nompar, [valpar], &
                    3, nomres, valres, icodre, 1)
        lambor(1) = valres(1)
        lambor(2) = valres(2)
        lambor(3) = valres(3)
        aniso = .true.
    else
        call utmess('F', 'ELEMENTS2_68')
    end if
!
!====
! PREALABLES LIES A L'ANISOTROPIE
!====
!
    epot = 0.d0
    do kp = 1, npg1
        l = (kp-1)*nno
        call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids, dfdx, dfdy, dfdz)
        flux = 0.d0
        fluy = 0.d0
        fluz = 0.d0
        do i = 1, nno
            flux = flux+zr(itempe-1+i)*dfdx(i)
            fluy = fluy+zr(itempe-1+i)*dfdy(i)
            fluz = fluz+zr(itempe-1+i)*dfdz(i)
        end do
!
        if (.not. aniso) then
            fluglo(1) = lambda(1)*flux
            fluglo(2) = lambda(1)*fluy
            fluglo(3) = lambda(1)*fluz
        else
            point(1) = 0.d0
            point(2) = 0.d0
            point(3) = 0.d0
            do ino = 1, nno
                point(1) = point(1)+zr(ivf+l+ino-1)*zr(igeom+3*ino-3)
                point(2) = point(2)+zr(ivf+l+ino-1)*zr(igeom+3*ino-2)
                point(3) = point(3)+zr(ivf+l+ino-1)*zr(igeom+3*ino-1)
            end do
            call rcangm(ndim, point, angmas)
            call matrot(angmas, p)
            fluglo(1) = flux
            fluglo(2) = fluy
            fluglo(3) = fluz
            call utpvgl(1, 3, p, fluglo, fluloc)
            fluloc(1) = lambor(1)*fluloc(1)
            fluloc(2) = lambor(2)*fluloc(2)
            fluloc(3) = lambor(3)*fluloc(3)
            call utpvlg(1, 3, p, fluloc, fluglo)
        end if
!
        epot = epot-(flux*fluglo(1)+fluy*fluglo(2)+fluz*fluglo(3))*poids
!
    end do
    zr(iener) = epot/2.d0
end subroutine
