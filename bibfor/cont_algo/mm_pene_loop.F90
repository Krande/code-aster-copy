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

subroutine mm_pene_loop(ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfdisl.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mminfi.h"
#include "asterfort/mminfl.h"
#include "asterfort/mminfr.h"
#include "asterfort/mminfm.h"
!
! person_in_charge: ayaovi-dzifa.kudawoo at edf.fr
!
    type(NL_DS_Contact), intent(inout) :: ds_contact
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Continue method - Management of contact loop
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  iter_newt        : index of current Newton iteration
! In  nume_inst        : index of current time step
! IO  ds_measure       : datastructure for measure and statistics management
! In  sddisc           : datastructure for time discretization
! In  disp_curr        : current displacements
! In  disp_cumu_inst   : displacement increment from beginning of current time
! IO  ds_contact       : datastructure for contact management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ztabf
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: jdecme, elem_slav_indx
    integer(kind=8) :: indi_cont_curr
    integer(kind=8) :: i_zone, i_elem_slav, i_cont_poin, i_poin_elem
    integer(kind=8) :: model_ndim, nb_cont_zone, loop_cont_vali
    integer(kind=8) ::  nb_poin_elem, nb_elem_slav
    real(kind=8) :: gap
    real(kind=8) :: time_curr
    aster_logical :: l_glis
    aster_logical :: l_veri, l_exis_glis
    aster_logical :: loop_cont_conv, l_loop_cont
    aster_logical :: l_frot_zone, l_pena_frot
    aster_logical :: l_frot, l_pena_cont
    integer(kind=8) :: type_adap, continue_calcul
    character(len=24) :: sdcont_cyccoe, sdcont_cyceta
    real(kind=8), pointer :: v_sdcont_cyccoe(:) => null()
    integer(kind=8), pointer :: v_sdcont_cyceta(:) => null()
    character(len=24) :: sdcont_tabfin, sdcont_jsupco, sdcont_apjeu
    real(kind=8), pointer :: v_sdcont_tabfin(:) => null()
    real(kind=8), pointer :: v_sdcont_jsupco(:) => null()
    real(kind=8), pointer :: v_sdcont_apjeu(:) => null()
    character(len=24) :: sdcont_pene
    real(kind=8), pointer :: v_sdcont_pene(:) => null()
    real(kind=8) ::  dist_max, vale_pene
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<CONTACT> ... ACTIVATION/DESACTIVATION'
    end if
!
! - Initializations
!
    l_glis = ASTER_FALSE
    l_veri = ASTER_FALSE
    l_exis_glis = ASTER_FALSE
    loop_cont_conv = ASTER_FALSE
    l_loop_cont = ASTER_FALSE
    l_frot_zone = ASTER_FALSE
    l_pena_frot = ASTER_FALSE
    l_frot = ASTER_FALSE
    l_pena_cont = ASTER_FALSE
    loop_cont_conv = ASTER_TRUE
    loop_cont_vali = 0
    continue_calcul = 0
    indi_cont_curr = 0
    gap = 0.0
    time_curr = 0.0
    dist_max = 1.0
    vale_pene = 1.0
!
! - Parameters
!
    l_exis_glis = cfdisl(ds_contact%sdcont_defi, 'EXIS_GLISSIERE')
    l_loop_cont = cfdisl(ds_contact%sdcont_defi, 'CONT_BOUCLE')
    type_adap = cfdisi(ds_contact%sdcont_defi, 'TYPE_ADAPT')
    model_ndim = cfdisi(ds_contact%sdcont_defi, 'NDIM')
    nb_cont_zone = cfdisi(ds_contact%sdcont_defi, 'NZOCO')
    l_frot = cfdisl(ds_contact%sdcont_defi, 'FROTTEMENT')
!
! - Acces to contact objects
!
    ztabf = cfmmvd('ZTABF')
    sdcont_tabfin = ds_contact%sdcont_solv(1:14)//'.TABFIN'
    sdcont_jsupco = ds_contact%sdcont_solv(1:14)//'.JSUPCO'
    sdcont_apjeu = ds_contact%sdcont_solv(1:14)//'.APJEU'
    call jeveuo(sdcont_tabfin, 'L', vr=v_sdcont_tabfin)
    call jeveuo(sdcont_jsupco, 'L', vr=v_sdcont_jsupco)
    call jeveuo(sdcont_apjeu, 'L', vr=v_sdcont_apjeu)
!
! - Acces to cycling objects
!
    sdcont_cyceta = ds_contact%sdcont_solv(1:14)//'.CYCETA'
    sdcont_cyccoe = ds_contact%sdcont_solv(1:14)//'.CYCCOE'
    call jeveuo(sdcont_cyceta, 'L', vi=v_sdcont_cyceta)
    call jeveuo(sdcont_cyccoe, 'L', vr=v_sdcont_cyccoe)

!
! - Access to penetration objects
!
    sdcont_pene = ds_contact%sdcont_solv(1:14)//'.PENETR'
    call jeveuo(sdcont_pene, 'E', vr=v_sdcont_pene)
!
!
! - Get current time
!
    time_curr = ds_contact%time_curr

    ds_contact%calculated_penetration = r8prem()
!
! - Loop on contact zones
!
    i_cont_poin = 1
    do i_zone = 1, nb_cont_zone
!
! ----- Parameters of zone
!
        l_glis = mminfl(ds_contact%sdcont_defi, 'GLISSIERE_ZONE', i_zone)
        l_veri = mminfl(ds_contact%sdcont_defi, 'VERIF', i_zone)
        nb_elem_slav = mminfi(ds_contact%sdcont_defi, 'NBMAE', i_zone)
        jdecme = mminfi(ds_contact%sdcont_defi, 'JDECME', i_zone)
        l_frot_zone = mminfl(ds_contact%sdcont_defi, 'FROTTEMENT_ZONE', i_zone)
        l_pena_frot = mminfl(ds_contact%sdcont_defi, 'ALGO_FROT_PENA', i_zone)
        l_pena_cont = mminfl(ds_contact%sdcont_defi, 'ALGO_CONT_PENA', i_zone)
        vale_pene = mminfr(ds_contact%sdcont_defi, 'PENE_MAXI', i_zone)

!
! ----- Loop on slave elements
!
        do i_elem_slav = 1, nb_elem_slav

!
! --------- Slave element index in contact datastructure
!
            elem_slav_indx = jdecme+i_elem_slav
!
! --------- Number of integration points on element
!
            call mminfm(elem_slav_indx, ds_contact%sdcont_defi, 'NPTM', nb_poin_elem)
!
! --------- Loop on integration points
!
            do i_poin_elem = 1, nb_poin_elem
!
! ------------- Get informations from
!
                gap = v_sdcont_apjeu(i_cont_poin)
                indi_cont_curr = nint(v_sdcont_tabfin(ztabf*(i_cont_poin-1)+23))
!
! Information for convergence : PENE_MAXI
                dist_max = vale_pene*ds_contact%arete_min
                if (l_pena_cont .and. indi_cont_curr .eq. 1) then
                    ds_contact%continue_pene = 0
                    v_sdcont_pene((i_cont_poin-1)+1) = abs(gap)
                    ! On cherche le max des penetrations
                    if (i_cont_poin .gt. 1) then
                        if (v_sdcont_pene((i_cont_poin-1)+1) .gt. &
                            v_sdcont_pene((i_cont_poin-1)+1-1)) then
                            ds_contact%calculated_penetration = &
                                v_sdcont_pene((i_cont_poin-1)+1)
                        else
                            ds_contact%calculated_penetration = &
                                v_sdcont_pene((i_cont_poin-1)+1-1)
                        end if
                    else
                        ds_contact%calculated_penetration = v_sdcont_pene((i_cont_poin-1)+1)
                    end if
                    if (ds_contact%iteration_newton .ge. ds_contact%it_adapt_maxi-6) then
                        ds_contact%continue_pene = 1.0
                    end if
!                    elseif (l_pena_cont .and. (indi_cont_curr .eq. 0) .and.&
!                            (ds_contact%calculated_penetration .le. 1.d-299)) then
!                    ! PENALISATION ADAPTATIF ACTIF mais pas encore de contact,
!                    ! la penetration n'a pas de sens
!                        ds_contact%continue_pene = 5
                end if

!
! ------------- Next contact point
!
                i_cont_poin = i_cont_poin+1
            end do
        end do
    end do
!    if (continue_calcul .eq. 10 ) ds_contact%calculated_penetration = r8prem()
!
! - Cleaning
!
    call jedema()
end subroutine
