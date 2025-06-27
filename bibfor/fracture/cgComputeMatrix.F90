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
! person_in_charge: tanguy.mathieu at edf.fr
!
subroutine cgComputeMatrix(cgField, cgTheta, cgStat)
!
    use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/gmatc3.h"
#include "asterfort/gmatr1.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    type(CalcG_field), intent(in)    :: cgField
    type(CalcG_theta), intent(inout) :: cgTheta
    type(CalcG_stat), intent(inout)  :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Compute A Matrix from equation A*G(s)=g(theta) in 2D and 3D
!
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)          ::  i, j
    character(len=24) :: chabsfon
    real(kind=8) :: start, finish
    real(kind=8), pointer :: matr(:) => null()
!
! --------------------------------------------------------------------------------------------------
    call cpu_time(start)
!
    call jemarq()

    cgTheta%matrix = '&&OP0060.MATRIX'

    if (cgField%ndim == 2) then
!       EN 2D, L'EQUATION EST SCALAIRE. A = 1
        call wkvect(cgTheta%matrix, 'V V R8', 1, vr=matr)
        matr(1) = 1.d0
    else if (cgField%ndim == 3) then
!
        if (cgTheta%discretization == 'LINEAIRE') then
!       Calcul de A dans le cas LINERAIRE
            call cgTheta%getAbsfonName(chabsfon)
            call gmatc3(cgTheta%nnof, cgTheta%milieu, cgTheta%l_closed, &
                        chabsfon, cgTheta%matrix)
        elseif (cgTheta%discretization == 'LEGENDRE') then
!
!       Dans le cas LEGENDRE, A = Identité. On laisse en commentaires pour l'instant l'appel
!       à la fonction historique qui calculait cette matrice au cas où ce ne serait pas l'idendité
!       dans un cas particulier non encore identifié.
!
!             call gmatr1(cgTheta%nb_fondNoeud,cgTheta%degree,cgTheta%absfond, &
!                         lonfis,cgTheta%matrix)
!
            call wkvect(cgTheta%matrix, 'V V R8', (cgTheta%degree+1)*(cgTheta%degree+1), vr=matr)
            do i = 1, cgTheta%degree+1
                do j = 1, cgTheta%degree+1
                    if (i == j) then
                        matr((i-1)*(cgTheta%degree+1)+j) = 1.d0
                    else
                        matr((i-1)*(cgTheta%degree+1)+j) = 0.d0
                    end if
                end do
            end do

        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jedema()
!
    call cpu_time(finish)
    cgStat%cgCmpMat = finish-start
end subroutine
