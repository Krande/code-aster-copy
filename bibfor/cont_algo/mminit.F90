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
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
subroutine mminit(mesh, ds_contact, sddyna, hat_valinc, ds_measure, &
                  sdnume, nume_inst, list_func_acti)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/copisd.h"
#include "asterfort/jeveuo.h"
#include "asterfort/misazl.h"
#include "asterfort/mm_cycl_init.h"
#include "asterfort/mmbouc.h"
#include "asterfort/mmopti.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/mmapin.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
    character(len=19), intent(in) :: hat_valinc(*)
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: sddyna
    character(len=19), intent(in) :: sdnume
    integer(kind=8), intent(in) :: nume_inst
    integer(kind=8), intent(in) :: list_func_acti(*)
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Initializations
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! IO  ds_contact       : datastructure for contact management
! In  hat_valinc       : hat variable for algorithm fields
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddyna           : datastructure for dynamic
! In  sdnume           : name of dof positions datastructure
! In  nume_inst        : index of current step time
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_dyna, l_cont_allv, l_step_first, l_geom_sans, l_pair
    character(len=19) :: sdcont_depgeo, sdcont_deplam, sdcont_depini
    character(len=19) :: sdcont_vitini, sdcont_accini
    character(len=19) :: disp_prev, acce_curr, vite_curr
    character(len=24) :: sdcont_tabfin
    real(kind=8), pointer :: v_sdcont_tabfin(:) => null()
    character(len=24) :: sdcont_etatct
    real(kind=8), pointer :: v_sdcont_etatct(:) => null()
    integer(kind=8) :: ztabf, zetat
    integer(kind=8) :: ipc, nb_inte_poin, ifm, niv
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'CONTACT5_13')
    end if
!
! - Initializations
!
    l_cont_allv = cfdisl(ds_contact%sdcont_defi, 'ALL_VERIF')
    l_geom_sans = cfdisl(ds_contact%sdcont_defi, 'REAC_GEOM_SANS')
    l_dyna = ndynlo(sddyna, 'DYNAMIQUE')
    ztabf = cfmmvd('ZTABF')
    zetat = cfmmvd('ZETAT')
    nb_inte_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
!
! - Using *_INIT options (like SEUIL_INIT)
!
    l_step_first = nume_inst .eq. 1
!
! - Geometric loop counter initialization
!
    call mmbouc(ds_contact, 'Geom', 'Init_Counter')
!
! - First geometric loop counter
!
    call mmbouc(ds_contact, 'Geom', 'Incr_Counter')
!
! - Get field names in hat-variables
!
    call nmchex(hat_valinc, 'VALINC', 'DEPMOI', disp_prev)
    call nmchex(hat_valinc, 'VALINC', 'VITPLU', vite_curr)
    call nmchex(hat_valinc, 'VALINC', 'ACCPLU', acce_curr)
!
! - Lagrangians initialized (LAMBDA TOTAUX)
!
    sdcont_depini = ds_contact%sdcont_solv(1:14)//'.INIT'
    call copisd('CHAMP_GD', 'V', disp_prev, sdcont_depini)
    call misazl(ds_contact, sdnume, disp_prev)
    if (l_dyna) then
        call misazl(ds_contact, sdnume, acce_curr)
        call misazl(ds_contact, sdnume, vite_curr)
    end if
!
! - Management of status for time cut
!
    sdcont_tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    sdcont_etatct = ds_contact%sdcont_solv(1:14)//'.ETATCT'
    call jeveuo(sdcont_tabfin, 'E', vr=v_sdcont_tabfin)
    call jeveuo(sdcont_etatct, 'L', vr=v_sdcont_etatct)

    do ipc = 1, nb_inte_poin
        v_sdcont_tabfin(ztabf*(ipc-1)+23) = v_sdcont_etatct(zetat*(ipc-1)+1)
        v_sdcont_tabfin(ztabf*(ipc-1)+17) = v_sdcont_etatct(zetat*(ipc-1)+2)
        v_sdcont_tabfin(ztabf*(ipc-1)+18) = v_sdcont_etatct(zetat*(ipc-1)+3)
    end do
!
! - Save speed and acceleration
!
    sdcont_vitini = ds_contact%sdcont_solv(1:14)//'.VITI'
    sdcont_accini = ds_contact%sdcont_solv(1:14)//'.ACCI'
    if (l_dyna) then
        call copisd('CHAMP_GD', 'V', vite_curr, sdcont_vitini)
        call copisd('CHAMP_GD', 'V', acce_curr, sdcont_accini)
    end if
!
! - Save displacements for geometric loop
!
    sdcont_depgeo = ds_contact%sdcont_solv(1:14)//'.DEPG'
    call copisd('CHAMP_GD', 'V', disp_prev, sdcont_depgeo)
!
! - Save displacements for friction loop
!
    sdcont_deplam = ds_contact%sdcont_solv(1:14)//'.DEPF'
    call copisd('CHAMP_GD', 'V', disp_prev, sdcont_deplam)
!
! - Update pairing ?
!
    l_pair = .not. l_geom_sans .or. (l_geom_sans .and. l_step_first)
!
! - Initial pairing
!
    if (l_pair) then
        call mmapin(mesh, ds_contact, ds_measure)
    end if
!
! - Initial options
!
    if (.not. l_cont_allv .and. l_step_first) then
        call mmopti(mesh, ds_contact, list_func_acti)
    end if
!
! - Cycling initialization
!
    call mm_cycl_init(ds_contact)
!
end subroutine
