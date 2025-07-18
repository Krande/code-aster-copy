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
subroutine modelCheck(model, lCheckJacobian, lCheckFSINorms, lCheckPlaneity)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/modexi.h"
#include "asterfort/utmess.h"
#include "asterfort/taxis.h"
#include "asterfort/modelCheckFSINormals.h"
#include "asterfort/modelCheckFluidFormulation.h"
#include "asterfort/modelCheckPlaneity.h"
!
    character(len=8), intent(in) :: model
    aster_logical, intent(in) :: lCheckJacobian, lCheckFSINorms, lCheckPlaneity
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_MODELE
!
! Check model
!
! --------------------------------------------------------------------------------------------------
!
! In  model           : name of the model
! In  lCheckJacobian  : .true. if check jacobian (element quality)
! In  lCheckFSINorms  : .true. if check normals for FSI elements
! In  lCheckPlaneity  : .true. if check planeity for plate elements
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_dim_geom, nb_dim_geom2, nb_dim_geom3
    character(len=16) :: repk
    integer(kind=8) :: i_disc_2d, i_disc_3d
    character(len=8) :: mesh
    aster_logical :: lAxis, lHHO
    integer(kind=8) :: nbCell
    character(len=19) :: modelLigrel
    integer(kind=8), pointer :: modelCells(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    modelLigrel = model//'.MODELE'

! - Get mesh support
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)

! - Parameters about model
    call dismoi('AXIS', model, 'MODELE', repk=repk)
    lAxis = repk .eq. 'OUI'
    call dismoi('DIM_GEOM', model, 'MODELE', repi=nb_dim_geom)
    call dismoi('EXI_HHO', modelLigrel, 'LIGREL', repk=repk)
    lHHO = repk .eq. 'OUI'
!
! - HHO should be alone
!
    if (lHHO) then
        call dismoi('EXI_NO_HHO', modelLigrel, 'LIGREL', repk=repk)
        if (repk .eq. 'OUI') then
            call utmess('F', 'MODELE1_10')
        end if
    end if
!
! - Check topoaster_logical dimensions
!
    if (nb_dim_geom .gt. 3) then
        nb_dim_geom2 = 0
        call utmess('A', 'MODELE1_14')
    else
        nb_dim_geom2 = 3
        nb_dim_geom3 = 3
        call dismoi('Z_CST', mesh, 'MAILLAGE', repk=repk)
        if (repk .eq. 'OUI') then
            nb_dim_geom2 = 2
            call dismoi('Z_ZERO', mesh, 'MAILLAGE', repk=repk)
            if (repk .eq. 'OUI') then
                nb_dim_geom3 = 2
            end if
        end if
        if ((nb_dim_geom .eq. 3) .and. (nb_dim_geom2 .eq. 2)) then
! --------- Correct: shells elements with Z=Constant
        else if ((nb_dim_geom .eq. 2) .and. (nb_dim_geom2 .eq. 3)) then
! --------- Warning: 2D model with 3D mesh
            call utmess('A', 'MODELE1_53')
        elseif ((nb_dim_geom .eq. 2) .and. &
                (nb_dim_geom2 .eq. 2) .and. &
                (nb_dim_geom3 .eq. 3)) then
! --------- Something strange: 2D with Z=Constant but Z<>0
            call utmess('A', 'MODELE1_58')
        end if
    end if
!
! - DISCRET elements: only 2D OR 3D
!
    call modexi(model, 'DIS_', i_disc_3d)
    call modexi(model, '2D_DIS_', i_disc_2d)
    if (nb_dim_geom2 .eq. 2 .and. i_disc_3d .eq. 1 .and. i_disc_2d .eq. 1) then
        call utmess('F', 'MODELE1_54')
    end if
!
! - Check if X>0 for axis elements
!
    if (lAxis) then
        call jeveuo(model//'.MAILLE', 'L', vi=modelCells)
        call taxis(mesh, modelCells, nbCell)
    end if
!
! - ON VERIFIE QUE LA GEOMETRIE DES MAILLES N'EST PAS TROP CHAHUTEE
!
    if (lCheckJacobian) then
        call calcul('C', 'VERI_JACOBIEN', modelLigrel, 1, mesh//'.COORDO', &
                    'PGEOMER', 1, '&&OP0018.CODRET', 'PCODRET', 'V', &
                    'OUI')
    end if

! - Check FSI norms
    if (lCheckFSINorms) then
        call modelCheckFSINormals(model)
    end if
! - Check that we have the same formulation in fluid modelisation
    call modelCheckFluidFormulation(model)

! - Check planeity of plate elements
    if (lCheckPlaneity) then
        call modelCheckPlaneity(mesh, model)
    end if
!
end subroutine
