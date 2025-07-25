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

subroutine cfprep(mesh, matr_asse, disp_iter, disp_cumu_inst, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfjein.h"
#include "asterfort/cfliin.h"
#include "asterfort/cfprch.h"
#include "asterfort/cfrsmu.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8), intent(in) :: mesh
    character(len=19), intent(in) :: matr_asse
    character(len=19), intent(in) :: disp_iter
    character(len=19), intent(in) :: disp_cumu_inst
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Preparation of contact solving
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  matr_asse        : matrix
! In  disp_iter        : displacement iteration
! In  disp_cumu_inst   : displacement increment from beginning of current time
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbliai, nt_cont_poin, model_ndim
    integer(kind=8) :: iliai
    integer(kind=8) :: lmat
    character(len=19) :: sdcont_liot
    integer(kind=8), pointer :: v_sdcont_liot(:) => null()
    character(len=19) :: sdcont_mu, sdcont_copo
    real(kind=8), pointer :: v_sdcont_mu(:) => null()
    real(kind=8), pointer :: v_sdcont_copo(:) => null()
    aster_logical :: l_pena_cont, l_frot, l_pena_frot
    aster_logical :: l_first_geom, l_pair
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ...... PREPARATION DU CALCUL'
    end if
!
! - Get parameters
!
    nbliai = cfdisd(ds_contact%sdcont_solv, 'NBLIAI')
    model_ndim = cfdisd(ds_contact%sdcont_solv, 'NDIM')
    nt_cont_poin = cfdisi(ds_contact%sdcont_defi, 'NTPC')
    l_pena_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_PENA')
    l_pena_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_PENA')
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROT_DISCRET')
!
! - Access to contact datastructures
!
    sdcont_liot = ds_contact%sdcont_solv(1:14)//'.LIOT'
    sdcont_mu = ds_contact%sdcont_solv(1:14)//'.MU'
    sdcont_copo = ds_contact%sdcont_solv(1:14)//'.COPO'
    call jeveuo(sdcont_liot, 'E', vi=v_sdcont_liot)
    call jeveuo(sdcont_mu, 'E', vr=v_sdcont_mu)
    call jeveuo(sdcont_copo, 'E', vr=v_sdcont_copo)
!
! - Get matrix access
!
    call jeveuo(matr_asse//'.&INT', 'E', lmat)
!
! - Get geometric loop state
!
    l_first_geom = ds_contact%l_first_geom
    l_pair = ds_contact%l_pair
!
! - Null pivot detection
!
    v_sdcont_liot(4*nt_cont_poin+1) = 0
    v_sdcont_liot(4*nt_cont_poin+2) = 0
    v_sdcont_liot(4*nt_cont_poin+3) = 0
    v_sdcont_liot(4*nt_cont_poin+4) = 0
!
! - Lagrange multipliers
!
    if (l_pena_cont .and. l_pena_frot .and. l_first_geom) then
        do iliai = 1, nt_cont_poin
            v_sdcont_mu(2*nt_cont_poin+iliai) = 0.d0
            v_sdcont_mu(nt_cont_poin+iliai) = 0.d0
        end do
    end if
!
    if (l_pena_cont) then
        do iliai = 1, nt_cont_poin
            v_sdcont_mu(iliai) = 0.d0
            if (l_pena_frot) then
                v_sdcont_mu(3*nt_cont_poin+iliai) = 0.d0
            end if
        end do
    end if
!
! - Restore Lagrange multiplier after pairing
!
    if (l_pair) then
        call cfrsmu(ds_contact, l_first_geom)
    end if
!
! - Prepare algorithm fields
!
    call cfprch(ds_contact, disp_iter, disp_cumu_inst)
!
! - Compute initial gaps
!
    call cfjein(mesh, ds_contact, disp_cumu_inst)
!
! - Set initial links
!
    call cfliin(mesh, ds_contact)
!
end subroutine
