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
subroutine thmGetParaElas(j_mater, kpi, temp, ndim, ds_thm)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/get_elasth_para.h"
#include "asterfort/utmess.h"
#include "asterfort/THM_type.h"
!
    integer(kind=8), intent(in) :: j_mater
    integer(kind=8), intent(in) :: kpi
    real(kind=8), intent(in) :: temp
    integer(kind=8), intent(in) :: ndim
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get elastic parameters
!
! --------------------------------------------------------------------------------------------------
!
! In  j_mater          : coded material address
! In  kpi              : current Gauss point
! In  temp             : current temperature
! In  ndim             : dimension of element (2 ou 3)
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: g, alpha(2)
    character(len=8) :: fami
!
! --------------------------------------------------------------------------------------------------
!
    fami = 'RIGI'
!
! - Get type of elasticity
!
    call get_elas_id(j_mater, ds_thm%ds_material%elas%id, ds_thm%ds_material%elas%keyword)
!
! - Read parameters
!
    call get_elas_para(fami, j_mater, '+', kpi, 1, &
                       ds_thm%ds_material%elas%id, ds_thm%ds_material%elas%keyword, &
                       temp=temp, &
                       e_=ds_thm%ds_material%elas%e, &
                       nu_=ds_thm%ds_material%elas%nu, &
                       e1_=ds_thm%ds_material%elas%e_l, &
                       e2_=ds_thm%ds_material%elas%e_t, &
                       e3_=ds_thm%ds_material%elas%e_n, &
                       nu12_=ds_thm%ds_material%elas%nu_lt, &
                       nu13_=ds_thm%ds_material%elas%nu_ln, &
                       nu23_=ds_thm%ds_material%elas%nu_tn, &
                       g1_=ds_thm%ds_material%elas%g_lt, &
                       g2_=ds_thm%ds_material%elas%g_ln, &
                       g3_=ds_thm%ds_material%elas%g_tn, &
                       g_=g)
    if (ds_thm%ds_material%elas%id .eq. 3) then
        ds_thm%ds_material%elas%g_ln = g
        ds_thm%ds_material%elas%g = g
    end if
!
! - Read parameters (dilatation)
!
    if (ds_thm%ds_elem%l_dof_ther) then
        call get_elasth_para(fami, j_mater, '+', kpi, 1, &
                             ds_thm%ds_material%elas%id, ds_thm%ds_material%elas%keyword, &
                             temp_vale_=temp, &
                             alpha=alpha, &
                             alpha_l=ds_thm%ds_material%ther%alpha_l, &
                             alpha_t=ds_thm%ds_material%ther%alpha_t, &
                             alpha_n=ds_thm%ds_material%ther%alpha_n)
        ds_thm%ds_material%ther%alpha = alpha(1)
    else
        ds_thm%ds_material%ther%alpha = 0.d0
    end if
!
! - Some checks: compatibility of elasticity with diffusion
!
    if (ds_thm%ds_material%biot%type .eq. BIOT_TYPE_ISOT) then
        if (ds_thm%ds_material%elas%id .ne. 1) then
            call utmess('F', 'THM1_2', sk=ds_thm%ds_material%elas%keyword)
        end if
    elseif (ds_thm%ds_material%biot%type .eq. BIOT_TYPE_ISTR) then
        if (ds_thm%ds_material%elas%id .ne. 3) then
            call utmess('F', 'THM1_2', sk=ds_thm%ds_material%elas%keyword)
        end if
    elseif (ds_thm%ds_material%biot%type .eq. BIOT_TYPE_ORTH) then
        if (ds_thm%ds_material%elas%id .ne. 2) then
            call utmess('F', 'THM1_2', sk=ds_thm%ds_material%elas%keyword)
        end if
    else
        ASSERT(.false.)
    end if
!
! - Some checks: anisotropy
!
    if (ds_thm%ds_material%elas%id .eq. 3) then
        if (ndim .ne. 3) then
            call utmess('F', 'THM1_4')
        end if
    end if
    if (ds_thm%ds_material%elas%id .eq. 2) then
        if (ndim .ne. 2) then
            call utmess('F', 'THM1_3')
        end if
    end if
!
end subroutine
