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
subroutine te0165(option, nomte)
!
    use Behaviour_module, only: behaviourOption
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/fpouli.h"
#include "asterfort/jevech.h"
#include "asterfort/kpouli.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "blas/ddot.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: CABLE_POULIE
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: icodre(2)
    real(kind=8) :: valres(2)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: nomres(2)
    real(kind=8) :: aire, w(9), nx, l1(3), l2(3), l10(3), l20(3)
    real(kind=8) :: e
    real(kind=8) :: norml1, norml2, norl10, norl20, l0, allong
    real(kind=8) :: preten, epsthe
    integer(kind=8) :: imatuu, ivectu, icontp
    integer(kind=8) :: lsect, igeom, imate, idepla, ideplp
    integer(kind=8) :: i, icoret, kc
    character(len=16) :: defo_comp, rela_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: codret
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    icontp = 1
    imatuu = 1
    ivectu = 1
    codret = 0
!
! - Get input fields
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCACABL', 'L', lsect)
    call jevech('PDEPLMR', 'L', idepla)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCOMPOR', 'L', vk16=compor)
!
! - Properties of behaviour
!
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    if (rela_comp(1:4) .ne. 'ELAS') then
        call utmess('F', 'CALCULEL4_92', sk=rela_comp)
    end if
    if (defo_comp .ne. 'GROT_GDEP') then
        call utmess('F', 'CALCULEL4_93', sk=defo_comp)
    end if
!
! - Get material properties
!
    nomres(1) = 'E'
    call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                ' ', 'ELAS', 0, '  ', [0.d0], &
                1, nomres, valres, icodre, 1)
    e = valres(1)
!
! - Get section properties
!
    aire = zr(lsect)
    preten = zr(lsect+1)
!
! - Thermal dilation
!
    call verift('RIGI', 1, 1, '+', zi(imate), &
                epsth_=epsthe)
!
! - Select objects to construct from option name
!
    call behaviourOption(option, compor, lMatr, lVect, lVari, &
                         lSigm, codret)
!
! - Get output fields
!
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
!
! - Update displacements
!
    do i = 1, 9
        w(i) = zr(idepla-1+i)+zr(ideplp-1+i)
    end do
!
    do kc = 1, 3
        l1(kc) = w(kc)+zr(igeom-1+kc)-w(6+kc)-zr(igeom+5+kc)
        l10(kc) = zr(igeom-1+kc)-zr(igeom+5+kc)
    end do
    do kc = 1, 3
        l2(kc) = w(3+kc)+zr(igeom+2+kc)-w(6+kc)-zr(igeom+5+kc)
        l20(kc) = zr(igeom+2+kc)-zr(igeom+5+kc)
    end do
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norml1 = ddot(b_n, l1, b_incx, l1, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norml2 = ddot(b_n, l2, b_incx, l2, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norl10 = ddot(b_n, l10, b_incx, l10, b_incy)
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    norl20 = ddot(b_n, l20, b_incx, l20, b_incy)
    norml1 = sqrt(norml1)
    norml2 = sqrt(norml2)
    norl10 = sqrt(norl10)
    norl20 = sqrt(norl20)
    l0 = norl10+norl20
    allong = (norml1+norml2-l0)/l0
    nx = e*aire*allong
!
    if (abs(nx) .le. 1.d-6) then
        nx = preten
    else
        nx = nx-e*aire*epsthe
    end if
!
    if (lMatr) then
        call kpouli(e, aire, nx, l0, l1, &
                    l2, norml1, norml2, zr(imatuu))
    end if
    if (lVect) then
        call fpouli(nx, l1, l2, norml1, norml2, &
                    zr(ivectu))
    end if
    if (lSigm) then
        zr(icontp) = nx
    end if
!
! - Save return code
!
    if (lSigm) then
        call jevech('PCODRET', 'E', icoret)
        zi(icoret) = codret
    end if
!
end subroutine
