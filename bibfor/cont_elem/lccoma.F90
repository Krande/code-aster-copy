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
subroutine lccoma(elem_dime, nb_node_mast, nb_node_slav, nb_lagr, &
                  l_norm_smooth, &
                  indi_lagc, poidpg, shape_mast_func, &
                  jaco_upda, dist_vect, &
                  mmat)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/jevech.h"
!
    integer(kind=8), intent(in) :: elem_dime
    integer(kind=8), intent(in) :: nb_node_mast
    integer(kind=8), intent(in) :: nb_node_slav
    integer(kind=8), intent(in) :: nb_lagr
    aster_logical, intent(in) :: l_norm_smooth
    integer(kind=8), intent(in) :: indi_lagc(10)
    real(kind=8), intent(in) :: poidpg
    real(kind=8), intent(in) :: shape_mast_func(9)
    real(kind=8), intent(in) :: jaco_upda, dist_vect(3)
    real(kind=8), intent(inout) :: mmat(55, 55)
!
! --------------------------------------------------------------------------------------------------
!
! Contact (LAC) - Elementary computations
!
! Compute contact matrix (master side)
!
! --------------------------------------------------------------------------------------------------
!
! In  elem_dime        : dimension of elements
! In  nb_node_mast     : number of nodes of for master side from contact element
! In  nb_node_slav     : number of nodes of for slave side from contact element
! In  nb_lagr          : total number of Lagrangian dof on contact element
! In  l_norm_smooth    : indicator for normals smoothing
! In  indi_lagc        : PREVIOUS node where Lagrangian dof is present (1) or not (0)
! In  poidspg          : weight at integration point
! In  shape_mast_func  : shape functions at integration point
! In  jaco_upda        : updated jacobian at integration point
! In  dist_vect        : distance vector between slave and master
! IO  mmat             : matrix
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_node_lagc, i_node_mast, i_dime, jj, indlgc, shift
    real(kind=8) :: r_nb_lagr
    integer(kind=8) :: jv_norm
!
! --------------------------------------------------------------------------------------------------
!
    shift = 0
    jj = 0
    r_nb_lagr = real(nb_lagr, kind=8)
!
    if (l_norm_smooth) then
        call jevech('PSNO', 'L', jv_norm)
        do i_node_lagc = 1, nb_node_slav
            shift = shift+indi_lagc(i_node_lagc)
            if (indi_lagc(i_node_lagc+1) .eq. 1) then
                indlgc = (i_node_lagc-1)*elem_dime+shift+elem_dime+1
                do i_node_mast = 1, nb_node_mast
                    do i_dime = 1, elem_dime
                        jj = (i_node_mast-1)*elem_dime+nb_node_slav*elem_dime+nb_lagr+i_dime
                        mmat(jj, indlgc) = mmat(jj, indlgc)+ &
                                           (zr(jv_norm+nb_node_slav*elem_dime+ &
                                               (i_node_mast-1)*elem_dime+i_dime-1)* &
                                            jaco_upda*poidpg*shape_mast_func(i_node_mast)) &
                                           /(r_nb_lagr)
                        mmat(indlgc, jj) = mmat(indlgc, jj)+ &
                                           (zr(jv_norm+nb_node_slav*elem_dime+ &
                                               (i_node_mast-1)*elem_dime+i_dime-1)* &
                                            jaco_upda*poidpg*shape_mast_func(i_node_mast)) &
                                           /(r_nb_lagr)
                    end do
                end do
            end if
        end do
    else
        do i_node_lagc = 1, nb_node_slav
            shift = shift+indi_lagc(i_node_lagc)
            if (indi_lagc(i_node_lagc+1) .eq. 1) then
                indlgc = (i_node_lagc-1)*elem_dime+shift+elem_dime+1
                do i_node_mast = 1, nb_node_mast
                    do i_dime = 1, elem_dime
                        jj = (i_node_mast-1)*elem_dime+nb_node_slav*elem_dime+nb_lagr+i_dime
                        mmat(jj, indlgc) = mmat(jj, indlgc)+ &
                                           (-dist_vect(i_dime)* &
                                            jaco_upda*poidpg* &
                                            shape_mast_func(i_node_mast))/(r_nb_lagr)
                        mmat(indlgc, jj) = mmat(indlgc, jj)+ &
                                           (-dist_vect(i_dime)* &
                                            jaco_upda*poidpg* &
                                            shape_mast_func(i_node_mast))/(r_nb_lagr)
                    end do
                end do
            end if
        end do
    end if
!
end subroutine
