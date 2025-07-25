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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine cfmxsd(mesh_, model_, nume_dof, sddyna, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfcrsd.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmap.h"
#include "asterfort/cfmmma.h"
#include "asterfort/cfmxme.h"
#include "asterfort/cfmxr0_lac.h"
#include "asterfort/cfmxr0.h"
#include "asterfort/infdbg.h"
#include "asterfort/lac_crsd.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: mesh_
    character(len=*), intent(in) :: model_
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: sddyna
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Prepare contact solving datastructure
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  list_func_acti   : list of active functionnalities
! In  sddyna           : name of dynamic solving datastructure
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=8) :: model, mesh
    integer(kind=8) :: nb_cont_zone
    aster_logical :: l_cont_disc, l_cont_cont, l_cont_allv, l_cont_lac
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_27')
    end if
!
! - Initializations
!
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    model = model_
    mesh = mesh_

! - Contact method
    l_cont_cont = cfdisl(ds_contact%sdcont_defi, 'FORMUL_CONTINUE')
    l_cont_disc = cfdisl(ds_contact%sdcont_defi, 'FORMUL_DISCRETE')
    l_cont_lac = cfdisl(ds_contact%sdcont_defi, 'FORMUL_LAC')
    l_cont_allv = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')

! - Create CONT_NOEU datastructure
    if (l_cont_cont .or. l_cont_disc) then
        call cfmxr0(mesh, ds_contact)
    end if

! - Create pairing datastructure
    if (l_cont_cont .or. l_cont_disc .or. l_cont_lac) then
        call cfmmap(mesh, ds_contact)
    end if
!
! - Create datastructures for solving
!
    if (.not. l_cont_allv) then

! ----- Create datastructures for DISCRETE/CONTINUE methods
        if (l_cont_cont .or. l_cont_disc) then
            call cfmmma(ds_contact)
        end if

! ----- Create datastructures for DISCRETE method
        if (l_cont_disc) then
            call cfcrsd(mesh, nume_dof, ds_contact)
        end if

! ----- Create datastructures for CONTINUE method
        if (l_cont_cont) then
            call cfmxme(nume_dof, sddyna, ds_contact)
        end if

! ----- Create datastructures for LAC method
        if (l_cont_lac) then
            call lac_crsd(nume_dof, ds_contact)
        end if

    end if

! - Create CONT_ELEM datastructure
    if (l_cont_lac) then
        call cfmxr0_lac(mesh, ds_contact)
    end if
!
end subroutine
