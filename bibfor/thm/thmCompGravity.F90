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
subroutine thmCompGravity(ds_thm)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/thmGetElemDime.h"
#include "asterfort/thmGetGeneDime.h"
#include "asterfort/thmGetElemInfo.h"
#include "asterfort/thmGetElemRefe.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/thmGetGene.h"
#include "asterfort/thmGetElemIntegration.h"
!
    type(THM_DS), intent(inout) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! Gravity (CHAR_MECA_PESA)
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: phenom
    real(kind=8) :: rho(1), coef, poids, rx
    integer(kind=8) :: icodre(1)
    integer(kind=8) :: jv_geom, jv_mater, jv_pesa, jv_vect
    character(len=8) :: elrefe, elref2
    integer(kind=8) :: nno, nnos, nnom
    integer(kind=8) :: npi, npi2, npg
    integer(kind=8) :: dimdep, dimdef, dimcon, dimuel
    integer(kind=8) :: nddls, nddlm
    integer(kind=8) :: nddl_meca, nddl_p1, nddl_p2, nddl_2nd
    integer(kind=8) :: jv_poids, jv_poids2
    integer(kind=8) :: jv_func, jv_func2, jv_dfunc, jv_dfunc2, jv_gano
    integer(kind=8) :: kpg, l, i, j, k, ii
    aster_logical :: l_vf, l_axi
    character(len=3) :: inte_type
    integer(kind=8) :: ndim
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
!
! --------------------------------------------------------------------------------------------------
!
!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm, l_axi, l_vf, ndim)
!
! - Cannot compute for finite volume
!
    ASSERT(.not. l_vf)
!
! - Get type of integration
!
    call thmGetElemIntegration(l_vf, inte_type)
!
! - Get generalized coordinates
!
    call thmGetGene(ds_thm, l_vf, ndim, &
                    mecani, press1, press2, tempe, second)
!
! - Get input/output fields
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PMATERC', 'L', jv_mater)
    call jevech('PPESANR', 'L', jv_pesa)
    call jevech('PVECTUR', 'E', jv_vect)
!
! - Get volumic mass
!
    call rccoma(zi(jv_mater), 'THM_DIFFU', 1, phenom, icodre(1))
    call rcvalb('FPG1', 1, 1, '+', zi(jv_mater), &
                ' ', phenom, 0, ' ', [0.d0], &
                1, 'RHO', rho, icodre, 1)
!
! - Get reference elements
!
    call thmGetElemRefe(l_vf, elrefe, elref2)
!
! - Get informations about element
!
    call thmGetElemInfo(l_vf, elrefe, elref2, &
                        nno, nnos, nnom, &
                        jv_gano, jv_poids, jv_poids2, &
                        jv_func, jv_func2, jv_dfunc, jv_dfunc2, &
                        inte_type, npi, npi2, npg)
    ASSERT(npi .le. 27)
    ASSERT(nno .le. 20)
!
! - Get dimensions of generalized vectors
!
    call thmGetGeneDime(ndim, &
                        mecani, press1, press2, tempe, second, &
                        dimdep, dimdef, dimcon)
!
! - Get dimensions about element
!
    call thmGetElemDime(ndim, nnos, nnom, &
                        mecani, press1, press2, tempe, second, &
                        nddls, nddlm, &
                        nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                        dimdep, dimdef, dimcon, dimuel)
!
! - Initializations
!
    do i = 1, dimuel
        zr(jv_vect+i-1) = 0.0d0
    end do
!
! - Compute
!
    nnom = nno-nnos
    if (ndim .eq. 3) then
        do kpg = 1, npg
            l = (kpg-1)*nno
            call dfdm3d(nno, kpg, jv_poids, jv_dfunc, zr(jv_geom), poids)
            coef = rho(1)*poids*zr(jv_pesa)
            do i = 1, nnos
                ii = nddls*(i-1)
                do j = 1, 3
                    zr(jv_vect+ii+j-1) = zr(jv_vect+ii+j-1)+ &
                                         coef*zr(jv_func+l+i-1)*zr(jv_pesa+j)
                end do
            end do
            do i = 1, nnom
                ii = nnos*nddls+nddlm*(i-1)
                do j = 1, 3
                    zr(jv_vect+ii+j-1) = zr(jv_vect+ii+j-1)+ &
                                         coef*zr(jv_func+l+i+nnos-1)*zr(jv_pesa+j)
                end do
            end do
        end do
    elseif (ndim .eq. 2) then
        do kpg = 1, npg
            k = (kpg-1)*nno
            call dfdm2d(nno, kpg, jv_poids, jv_dfunc, zr(jv_geom), poids)
            poids = poids*rho(1)*zr(jv_pesa)
            if (l_axi) then
                rx = 0.d0
                do i = 1, nno
                    rx = rx+zr(jv_geom+2*i-2)*zr(jv_func+k+i-1)
                end do
                poids = poids*rx
                do i = 1, nnos
                    zr(jv_vect+nddls*(i-1)-1+2) = zr(jv_vect+nddls*(i-1)-1+2)+ &
                                                  poids*zr(jv_pesa+2)*zr(jv_func+k+i-1)
                end do
                do i = 1, nnom
                    zr(jv_vect+nddls*nnos+nddlm*(i-1)-1+2) = zr(jv_vect+nddls*nnos+nddlm*(i-1)+1)+ &
                                                             poids* &
                                                             zr(jv_pesa+2)*zr(jv_func+k+i+nnos-1)
                end do
            else
                do i = 1, nnos
                    zr(jv_vect+nddls*(i-1)-1+1) = zr(jv_vect+nddls*(i-1)-1+1)+ &
                                                  poids*zr(jv_pesa+1)*zr(jv_func+k+i-1)
                    zr(jv_vect+nddls*(i-1)-1+2) = zr(jv_vect+nddls*(i-1)-1+2)+ &
                                                  poids*zr(jv_pesa+2)*zr(jv_func+k+i-1)
                end do
                do i = 1, nnom
                    zr(jv_vect+nddls*nnos+nddlm*(i-1)) = zr(jv_vect+nddls*nnos+nddlm*(i-1))+ &
                                                         poids* &
                                                         zr(jv_pesa+1)*zr(jv_func+k+i+nnos-1)
                    zr(jv_vect+nddls*nnos+nddlm*(i-1)+1) = zr(jv_vect+nddls*nnos+nddlm*(i-1)+1)+ &
                                                           poids* &
                                                           zr(jv_pesa+2)*zr(jv_func+k+i+nnos-1)
                end do
            end if
        end do
    else
        ASSERT(.false.)
    end if
!
end subroutine
