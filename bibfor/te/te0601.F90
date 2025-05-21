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
subroutine te0601(option, nomte)
!
    use raco3d_utils
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/writeMatrix.h"
#include "asterfort/rco3d_elem.h"
#include "asterfort/rco3d_infos.h"
#include "asterfort/rco3d_calcmat.h"
#include "jeveux.h"


!
    character(len=16) :: nomte, option
!
!
! --------------------------------------------------------------------------------------------------
!
! SHELL-3D link 
!
! Link elementary matrix
! --------------------------------------------------------------------------------------------------
!
    integer :: jv_geom, jv_cacoqu, nddl
    integer :: nnco, nn3d, dim, index
    integer :: i_node_co, i_node_3d, i
    integer :: nb_gauss, ncols, nrows
    real(kind=8) :: jac_det(NB_GAUSS_MAX), gauss_weight(NB_GAUSS_MAX) 
    real(kind=8) :: gauss_coor(2, NB_GAUSS_MAX)
    real(kind=8) :: ff_co(3, NB_GAUSS_MAX), epai, ff_3d(8,NB_GAUSS_MAX)
    real(kind=8) :: t(3,NB_GAUSS_MAX), n(3,NB_GAUSS_MAX), s(3)
    character(len=8):: typmaco, typma3d
    real(kind=8), allocatable :: mat(:, :)
    aster_logical :: skip(NB_GAUSS_MAX)


    ! RECUPERER GEOMETRIE
    call jevech('PGEOMER', 'L', jv_geom)

    ! RECUPERER EPAISSEUR
    call jevech('PCACOQU', 'L', jv_cacoqu)
    epai = zr(jv_cacoqu-1+1)


    ! RECUPERER LES INFOS SUR L'ELEMENT
    call rco3d_elem(nomte, dim, nddl, typmaco, nnco, typma3d, nn3d)
    nrows = 6 * nnco
    ncols = 6 * nnco + 3 * nn3d 
    !
    
    ! RECUPERER LES PTS GAUSS 
    call rco3d_infos(typmaco, typma3d, epai, jv_geom, nb_gauss, gauss_coor, &
            gauss_weight, jac_det, ff_co, ff_3d, s, t, n, skip)

    ! allocation and calculation of the matrix 
    allocate(mat(nrows,ncols))
    ! initialization
    mat = 0.0d0
    !
    call  rco3d_calcmat(nb_gauss, gauss_weight, gauss_coor, jac_det, &
                    ff_co, ff_3d, s, t, n, epai, & 
                        nnco, nn3d, skip, mat )
    

    ! copy the elementary matrix
    call writeMatrix('PMATUNS', nrows, ncols, ASTER_FALSE, mat)

    deallocate(mat)


end subroutine