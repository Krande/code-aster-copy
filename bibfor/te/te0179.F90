! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine te0179(option, nomte)
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/getFluidPara.h"
#include "asterfort/vff2dn.h"
!
character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: ACOU / PLAN (boundary)
!
! Options: CHAR_ACOU_VFAC
!
! --------------------------------------------------------------------------------------------------
!
    integer :: jv_geom, jv_mate, jv_speed, jv_vect
    real(kind=8) :: nx, ny
    real(kind=8) :: rho, poids
    complex(kind=8) :: vnor
    integer :: ipoids, ivf, idfde
    integer :: nno, npg, ndim, ndof
    integer :: ldec
    integer :: i, ipg
    aster_logical :: l_axis
    real(kind=8) :: r
    integer :: j_mater
!
! --------------------------------------------------------------------------------------------------
!

! - Input fields
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
    call jevech('PVITEFC', 'L', jv_speed)

! - Get element parameters
    l_axis = (lteatt('AXIS','OUI'))
    call elrefe_info(fami='RIGI',&
                     nno=nno, npg=npg, ndim=ndim,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ndof = nno

! - Get material properties
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho)

! - Output field
    call jevech('PVECTTC', 'E', jv_vect)
    do i = 1, ndof
        zc(jv_vect+i-1) = (0.d0, 0.d0)
    end do

! - Loop on Gauss points
    do ipg = 1, npg
        ldec = (ipg-1)*nno

! ----- Compute normal
        nx = 0.d0
        ny = 0.d0
        call vff2dn(ndim, nno, ipg, ipoids, idfde,&
                    zr(jv_geom), nx, ny, poids)
        if (l_axis) then
            r = 0.d0
            do i = 1, nno
                r = r + zr(jv_geom+2*(i-1))*zr(ivf+ldec+i-1)
            end do
            poids = poids*r
        endif

! ----- Get value of normal speed
        vnor = zc(jv_speed+ipg-1)

! ----- Compute vector
        do i = 1, nno
            zc(jv_vect+i-1) = zc(jv_vect+i-1) +&
                              poids *&
                              zr(ivf+ldec+i-1) * vnor * rho
        end do

    end do
!
end subroutine
