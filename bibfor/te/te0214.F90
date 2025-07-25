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
subroutine te0214(nomopt, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecach.h"
#include "blas/ddot.h"
!
    character(len=16) :: nomopt, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'AMOR_MECA'
!                                OU 'RIGI_MECA_HYST'
!        POUR TOUS LES TYPES D'ELEMENTS JOINTS 2D
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer(kind=8) :: nbres, nbpar
    parameter(nbres=4)
    parameter(nbpar=3)
!
    integer(kind=8) :: iret, nbval, nbddl
    integer(kind=8) :: i, j, irigi, jma, i2, i3
    integer(kind=8) :: iresu, ins, irns, ivari
    integer(kind=8) :: idresu(5), idrigi(2)
    integer(kind=8) :: igeom
!
    real(kind=8) :: valres(nbres), valpar(nbpar)
!
    real(kind=8) :: x(4), y(4), z(4), c1(3), c2(3)
    real(kind=8) :: surf
    integer(kind=8) :: icodre(nbres)
    character(len=16) :: nomres(nbres)
    aster_logical :: ljfr
    blas_int :: b_incx, b_incy, b_n
!
!     -- RECUPERATION DES CHAMPS PARAMETRES ET DE LEURS LONGUEURS:
!     ------------------------------------------------------------
    ins = 0
    irns = 0
    if (nomopt .eq. 'AMOR_MECA') then
        call tecach('NNO', 'PRIGIEL', 'L', ins, iad=idrigi(1))
        if (ins .eq. 0) then
            call tecach('ONO', 'PMATUUR', 'E', iret, nval=5, &
                        itab=idresu)
        else
            call tecach('NNO', 'PMATUNS', 'E', irns, nval=5, &
                        itab=idresu)
            if (irns .ne. 0) call tecach('ONO', 'PMATUUR', 'E', iret, 5, &
                                         itab=idresu)
        end if
    else if (nomopt .eq. 'RIGI_MECA_HYST') then
        call tecach('ONO', 'PMATUUC', 'E', iret, nval=5, &
                    itab=idresu)
    else
        ASSERT(.false.)
    end if
    nbval = idresu(2)
!
    ljfr = .false.
    call tecach('NNO', 'PMATERC', 'L', iret, iad=jma)
    if ((jma .eq. 0) .or. (iret .ne. 0)) goto 1
    nomres(1) = 'K_N'
    nomres(2) = 'AMOR_NOR'
    nomres(3) = 'AMOR_TAN'
    nomres(4) = 'COEF_AMOR'
    valres(1) = 0.d0
    valres(2) = 0.d0
    valres(3) = 0.d0
    valres(4) = 0.d0
    call rcvala(zi(jma), ' ', 'JOINT_MECA_FROT', 0, ' ', &
                [valpar], 4, nomres, valres, icodre, &
                0)
    if (icodre(1) .eq. 0) then
        ljfr = .true.
    end if
    do i = 2, 4
        if (icodre(i) .ne. 0) then
            valres(i) = 0.d0
        end if
    end do
!
    call jevech('PGEOMER', 'L', igeom)
!
    do i = 1, 4
        x(i) = zr(igeom+2*(i-1)+1-1)
        y(i) = zr(igeom+2*(i-1)+2-1)
        z(i) = 0.d0
    end do
    c2(1) = x(2)-x(1)
    c2(2) = y(2)-y(1)
    b_n = to_blas_int(2)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    surf = ddot(b_n, c2, b_incx, c2, b_incy)
    c2(1) = c2(1)/sqrt(surf)
    c2(2) = c2(2)/sqrt(surf)
    surf = sqrt(surf)
    c1(1) = -c2(2)
    c1(2) = c2(1)
!
!
!     -- CALCUL PROPREMENT DIT :
!     --------------------------
    iresu = idresu(1)
    irigi = idrigi(1)
    if (nomopt .eq. 'AMOR_MECA') then
        if (ljfr) then
            if (ins .eq. 0) then
                call tecach('ONO', 'PRIGINS', 'L', irns, iad=idrigi(1))
                call tecach('ONO', 'PRIGIEL', 'L', iret, nval=2, &
                            itab=idrigi)
                nbddl = int(-1.0d0+sqrt(1.0d0+8.d0*dble(idrigi(2))))/2
                nbval = idrigi(2)
                call tecach('ONO', 'PMATUUR', 'E', iret, nval=2, &
                            itab=idresu)
            else
                call tecach('ONO', 'PRIGINS', 'L', iret, nval=2, &
                            itab=idrigi)
                call tecach('NNO', 'PMATUNS', 'E', irns, nval=5, &
                            itab=idresu)
                if (irns .ne. 0) then
                    call tecach('ONO', 'PMATUUR', 'E', iret, 5, &
                                itab=idresu)
                    nbddl = int(-1.0d0+sqrt(1.0d0+8.d0*dble(idresu(2))))/2
                else
                    nbddl = int(sqrt(dble(idresu(2))))
                end if
                nbval = idresu(2)
            end if
            irigi = idrigi(1)
            iresu = idresu(1)
            call jevech('PVARIPG', 'L', ivari)
            if (irigi .ne. 0) then
                if (ins .ne. 0 .and. irns .ne. 0) then
                    do i3 = 1, 2
                        do j = 1, 2
                            i = 2*(i3-1)+j
                            i2 = i+10-4*i3
                            zr(iresu-1+i*(i+1)/2) = valres(2)*c1(j)**2+valres(3)*c2(j)**2
                            zr(iresu-1+i*(i+1)/2) = zr(iresu-1+i*(i+1)/2)*surf/2.0d0
                            if (zr(ivari-1+7) .ge. 0.d0) then
                                zr(iresu-1+i*(i+1)/2) = zr(iresu-1+i*(i+1)/2)*valres(4)
                            end if
                            zr(iresu-1+i2*(i2+1)/2) = zr(iresu-1+i*(i+1)/2)
                            zr(iresu-1+i+i2*(i2-1)/2) = -zr(iresu-1+i*(i+1)/2)
                        end do
                    end do
                else
                    do i = 1, nbval
                        zr(iresu-1+i) = zr(irigi-1+i)*valres(2)/valres(1)
                    end do
                end if
            end if
            goto 1
        end if
    end if
1   continue
end subroutine
