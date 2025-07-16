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
subroutine te0381(option, nomte)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/vff2dn.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
#include "asterfort/getFluidPara.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: AXIS_FLUI_STRU, 2D_FLUI_STRU
!
! Option: AMOR_MECA
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: mmat(9, 9)
    real(kind=8) :: nx, ny, norm(2)
    real(kind=8) :: poids, rho
    integer(kind=8) :: jv_geom, jv_mate, jv_matr
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: nno, npg, ndim
    integer(kind=8) :: ij
    integer(kind=8) :: ino1, ino2, ipg, ind1, ind2, idim
    integer(kind=8) :: ldec
    integer(kind=8) :: j_mater, iret
    character(len=16) :: FEForm
    aster_logical :: l_axis
    real(kind=8) :: r
!
! --------------------------------------------------------------------------------------------------
!
    mmat = 0.d0
!
! - Input fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mate)
!
! - Get element parameters
!
    call teattr('S', 'FORMULATION', FEForm, iret)
    l_axis = (lteatt('AXIS', 'OUI'))
    call elrefe_info(fami='RIGI', &
                     nno=nno, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. 3)
!
! - Get material properties for fluid
!
    j_mater = zi(jv_mate)
    call getFluidPara(j_mater, rho_=rho)
!
! - Loop on Gauss points
!
    do ipg = 1, npg
        ldec = (ipg-1)*nno
        call vff2dn(ndim, nno, ipg, ipoids, idfde, &
                    zr(jv_geom), nx, ny, poids)
        norm(1) = nx
        norm(2) = ny
        if (l_axis) then
            r = 0.d0
            do ino1 = 1, nno
                r = r+zr(jv_geom+2*(ino1-1))*zr(ivf+ldec+ino1-1)
            end do
            poids = poids*r
        end if
        if (FEForm .eq. 'U_PSI') then
            do ino1 = 1, nno
                do ino2 = 1, nno
                    do idim = 1, 2
                        ind1 = 3*(ino1-1)+idim
                        ind2 = 3*(ino2-1)+3
                        mmat(ind2, ind1) = mmat(ind2, ind1)+ &
                                           poids*norm(idim)*rho* &
                                           zr(ivf+ldec+ino1-1)*zr(ivf+ldec+ino2-1)
                        mmat(ind1, ind2) = mmat(ind1, ind2)+ &
                                           poids*norm(idim)*rho* &
                                           zr(ivf+ldec+ino1-1)*zr(ivf+ldec+ino2-1)
                    end do
                end do
            end do
        else
            call utmess('F', 'FLUID1_2', sk=FEForm)
        end if
    end do
!
! - Output field
!
    if (FEForm .eq. 'U_PSI') then
        call jevech('PMATUUR', 'E', jv_matr)
        do ino2 = 1, 3*nno
            do ino1 = 1, ino2
                ij = (ino2-1)*ino2/2+ino1
                zr(jv_matr+ij-1) = mmat(ino1, ino2)
            end do
        end do
    else
        call utmess('F', 'FLUID1_2', sk=FEForm)
    end if
!
end subroutine
