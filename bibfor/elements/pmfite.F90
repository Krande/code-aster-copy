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

subroutine pmfite(typfib, nf, ncarf, vf, ve, vs)
!
!
! --------------------------------------------------------------------------------------------------
!
!       INTEGRATIONS SUR LA SECTION (TENANT COMPTE DU MODULE DE CHAQUE FIBRE)
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!
!   IN
!       typfib  : type des fibres : 1 ou 2
!       nf      : nombre de fibres
!       ncarf   : nombre de caracteristiques sur chaque fibre
!       vf(*)   : positions des fibres
!           Types 1 et 2
!               vf(1,*) : Y fibres
!               vf(2,*) : Z fibres
!               vf(3,*) : Aire fibres
!           Types 2
!               vf(4,*) : Yp groupes de fibres
!               vf(5,*) : Zp groupes de fibres
!               vf(6,*) : num du groupe
!       ve(*)   : module des fibres
!
!   OUT
!       vs(1) : int(e.ds)
!       vs(2) : int(e.y.ds)
!       vs(3) : int(e.z.ds)
!       vs(4) : int(e.y.y.ds)
!       vs(5) : int(e.z.z.ds)
!       vs(6) : int(e.y.z.ds)
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterfort/utmess.h"
!
    integer(kind=8) :: typfib, nf, ncarf
    real(kind=8) :: vf(ncarf, nf), ve(nf), vs(6)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii
    real(kind=8) :: eds, yy, zz, aire
!
! --------------------------------------------------------------------------------------------------
!
    vs(:) = 0.0d0
!
    if (typfib .eq. 1) then
!       3 caractéristiques utiles par fibre : y z aire
        do ii = 1, nf
            yy = vf(1, ii)
            zz = vf(2, ii)
            aire = vf(3, ii)
!
            eds = aire*ve(ii)
            vs(1) = vs(1)+eds
            vs(2) = vs(2)+yy*eds
            vs(3) = vs(3)+zz*eds
            vs(4) = vs(4)+yy*yy*eds
            vs(5) = vs(5)+zz*zz*eds
            vs(6) = vs(6)+yy*zz*eds
        end do
    else if (typfib .eq. 2) then

!       6 caractéristiques utiles par fibre : y z aire yp zp numgr
        do ii = 1, nf
            yy = vf(1, ii)
            zz = vf(2, ii)
            aire = vf(3, ii)
!
            eds = aire*ve(ii)
            vs(1) = vs(1)+eds
            vs(2) = vs(2)+yy*eds
            vs(3) = vs(3)+zz*eds
            vs(4) = vs(4)+yy*yy*eds
            vs(5) = vs(5)+zz*zz*eds
            vs(6) = vs(6)+yy*zz*eds
        end do
    else
        call utmess('F', 'ELEMENTS2_40', si=typfib)
    end if
!
end subroutine
