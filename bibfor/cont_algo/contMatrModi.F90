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
subroutine contMatrModi(modelZ, ds_contact, matrAsse)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/cfdisl.h"
#include "asterfort/dismoi.h"
#include "asterfort/echmat.h"
#include "asterfort/exisd.h"
#include "asterfort/isParallelMatrix.h"
!
    character(len=*), intent(in) :: modelZ
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=19), intent(in) :: matrAsse
!
! --------------------------------------------------------------------------------------------------
!
! Modification of contact matrix at prediction
!
! --------------------------------------------------------------------------------------------------
!
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iexi
    character(len=8) :: partit
    aster_logical :: l_contact_adapt, ldist, l_parallel_matrix
    real(kind=8) :: minmat, maxmat, exponent_val
!
! --------------------------------------------------------------------------------------------------
!
    minmat = 0.d0
    maxmat = 0.d0
    exponent_val = 0.d0
!   -- Avant la factorisation et pour le cas ou il y a du contact continu avec adaptation de
!      coefficient
!   -- On cherche le coefficient optimal pour eviter une possible singularite de matrice
!   -- La valeur est estimee une seule fois a la premiere prediction du premier pas de
!      temps pour l'etape de calcul
!   -- Cette valeur estimee est passee directement a mmchml_c sans passer par mmalgo car
!   -- a la premiere iteration on ne passe pas par mmalgo
    l_contact_adapt = cfdisl(ds_contact%sdcont_defi, 'EXIS_ADAP')
    if ((nint(ds_contact%update_init_coefficient) .eq. 0) .and. l_contact_adapt) then
        l_parallel_matrix = isParallelMatrix(matrAsse)
        call dismoi('PARTITION', modelZ, 'MODELE', repk=partit)
        call exisd('PARTITION', partit, iexi)
        ldist = iexi .ne. 0
        call echmat(matrAsse, ldist, l_parallel_matrix, minmat, maxmat)
        ds_contact%max_coefficient = maxmat
        if (abs(log(minmat)) .ge. r8prem()) then
            if (abs(log(maxmat))/abs(log(minmat)) .lt. 4.0d0) then
!                     Le rapport d'arete max/min est
                !  un bon compromis pour initialiser le coefficient
                ds_contact%estimated_coefficient = &
                    ((1.D3*ds_contact%arete_max)/(1.D-2*ds_contact%arete_min))
                ds_contact%update_init_coefficient = 1.0d0
            else
                exponent_val = min(abs(log(minmat)), abs(log(maxmat)))/10.d0
                ds_contact%estimated_coefficient = 10.d0**(exponent_val)
                ds_contact%update_init_coefficient = 1.0d0
            end if
        else
            ds_contact%estimated_coefficient = 1.d16*ds_contact%arete_min
            ds_contact%update_init_coefficient = 1.0d0
        end if
    end if
!
end subroutine
