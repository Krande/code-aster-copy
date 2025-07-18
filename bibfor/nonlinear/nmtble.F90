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
subroutine nmtble(loop_exte, model, mesh, ds_contact, &
                  list_func_acti, ds_print, &
                  sderro, ds_conv, sddisc, nume_inst, hval_incr, &
                  hval_algo, ds_algorom)
!
    use NonLin_Datastructure_type
    use Rom_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/isfonc.h"
#include "asterfort/mmbouc.h"
#include "asterfort/mm_cycl_erase.h"
#include "asterfort/nmaffi.h"
#include "asterfort/nmctcc.h"
#include "asterfort/nmctcf.h"
#include "asterfort/nmctgo.h"
#include "asterfort/nmevcv.h"
#include "asterfort/nmimci.h"
#include "asterfort/nmleeb.h"
#include "asterfort/nmcrel.h"
!
    integer(kind=8), intent(inout) :: loop_exte
    character(len=24), intent(in) :: model
    character(len=8), intent(in) :: mesh
    type(NL_DS_Contact), intent(inout) :: ds_contact
    integer(kind=8), intent(in) :: list_func_acti(*)
    type(NL_DS_Print), intent(inout) :: ds_print
    character(len=24), intent(in) :: sderro
    type(NL_DS_Conv), intent(in) :: ds_conv
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: hval_algo(*)
    type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algo
!
! External loop management - END
!
! --------------------------------------------------------------------------------------------------
!
! IO  loop_exte        : level of external loop (see nmible.F90)
!                        0 - Not use (not contact)
!                        1 - Loop for contact status
!                        2 - Loop for friction triggers
!                        3 - Loop for geometry
!                       10 - External loop for HROM
! In  model            : name of model
! In  mesh             : name of mesh
! In  ds_material      : datastructure for material parameters
! IO  ds_contact       : datastructure for contact management
! In  list_func_acti   : list of active functionnalities
! IO  ds_print         : datastructure for printing parameters
! In  sderro           : datastructure for errors during algorithm
! In  ds_conv          : datastructure for convergence management
! In  sddisc           : datastructure for time discretization
! In  nume_inst        : index of current step time
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  ds_constitutive  : datastructure for constitutive laws management
! IO  ds_algorom       : datastructure for ROM parameters
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: loop_cont_conv, loop_fric_conv, loop_geom_conv
    aster_logical :: l_loop_frot, l_loop_geom, l_loop_cont, l_cont_cont
    integer(kind=8) :: loop_geom_count, loop_fric_count, loop_cont_count, loop_cont_vali
    character(len=4) :: state_newt
    real(kind=8) :: loop_cont_vale
!
! --------------------------------------------------------------------------------------------------
!
    if ((loop_exte .ge. 1) .and. (loop_exte .le. 3)) then
!
! ----- State of Newton loop
!
        call nmleeb(sderro, 'NEWT', state_newt)
!
! ----- To evaluate contact loops: Newton has covnerged
!
        if (state_newt .ne. 'CONV') then
            goto 999
        end if
!
! ----- Contact loops
!
        l_cont_cont = isfonc(list_func_acti, 'CONT_CONTINU')
        l_loop_frot = isfonc(list_func_acti, 'BOUCLE_EXT_FROT')
        l_loop_geom = isfonc(list_func_acti, 'BOUCLE_EXT_GEOM')
        l_loop_cont = isfonc(list_func_acti, 'BOUCLE_EXT_CONT')
!
! ----- Initializations
!
        loop_cont_vali = 0
        loop_cont_conv = .false.
        loop_fric_conv = .false.
        loop_geom_conv = .false.
!
! ----- <1> - Contact loop
!
        if (loop_exte .le. 1) then
            if (l_loop_cont) then
                loop_exte = 1
                call nmctcc(mesh, model, nume_inst, &
                            sderro, sddisc, hval_incr, hval_algo, &
                            ds_contact)
                call mmbouc(ds_contact, 'Cont', 'Is_Convergence', loop_state_=loop_cont_conv)
                call mmbouc(ds_contact, 'Cont', 'Get_Vale', loop_vale_=loop_cont_vale)
                loop_cont_vali = nint(loop_cont_vale)
                if (.not. loop_cont_conv) then
                    loop_exte = 1
                    goto 500
                end if
            end if
        end if
!
! ----- <2> - Friction loop
!
        if (loop_exte .le. 2) then
            if (l_loop_frot) then
                loop_exte = 2
                call nmctcf(mesh, sderro, hval_incr, ds_print, ds_contact)
                call mmbouc(ds_contact, 'Fric', 'Is_Convergence', loop_state_=loop_fric_conv)
                if (.not. loop_fric_conv) then
                    loop_exte = 2
                    goto 500
                end if
            end if
        end if
!
! ----- <3> - Geometric loop
!
        if (loop_exte .le. 3) then
            if (l_loop_geom) then
                loop_exte = 3
                call nmctgo(mesh, sderro, hval_incr, ds_print, ds_contact)
                call mmbouc(ds_contact, 'Geom', 'Is_Convergence', loop_state_=loop_geom_conv)
                if (.not. loop_geom_conv) then
                    loop_exte = 3
                    goto 500
                end if
            end if
        end if
!
500     continue
!
! ----- Initialization of data structures for cycling detection and treatment
!
        if ((loop_cont_conv .or. loop_fric_conv .or. loop_geom_conv) .and. l_cont_cont) then
            call mm_cycl_erase(ds_contact, 0, 0)
        end if
!
! ----- Print line
!
        call nmaffi(list_func_acti, ds_conv, ds_print, sderro, sddisc, &
                    'FIXE')
!
! ----- New iteration in loops
!
        if (.not. loop_cont_conv .and. loop_exte .eq. 1) then
            call mmbouc(ds_contact, 'Cont', 'Incr_Counter')
        end if
        if (.not. loop_fric_conv .and. loop_exte .eq. 2) then
            call mmbouc(ds_contact, 'Fric', 'Incr_Counter')
        end if
        if (.not. loop_geom_conv .and. loop_exte .eq. 3) then
            call mmbouc(ds_contact, 'Geom', 'Incr_Counter')
        end if
!
! ----- Update loops index
!
        call mmbouc(ds_contact, 'Cont', 'Read_Counter', loop_cont_count)
        call mmbouc(ds_contact, 'Fric', 'Read_Counter', loop_fric_count)
        call mmbouc(ds_contact, 'Geom', 'Read_Counter', loop_geom_count)
!
! ----- Print management
!
        call nmimci(ds_print, 'CONT_NEWT', loop_cont_vali, .true._1)
        call nmimci(ds_print, 'BOUC_CONT', loop_cont_count, .true._1)
        call nmimci(ds_print, 'BOUC_FROT', loop_fric_count, .true._1)
        call nmimci(ds_print, 'BOUC_GEOM', loop_geom_count, .true._1)
    elseif (loop_exte .eq. 10) then
        if (ds_algorom%phase .eq. 'HROM') then
            call nmaffi(list_func_acti, ds_conv, ds_print, sderro, sddisc, &
                        'FIXE')
            ds_algorom%phase = 'CORR_EF'
            call nmcrel(sderro, 'DIVE_FIXG', .true._1)
        else if (ds_algorom%phase .eq. 'CORR_EF') then
            call nmaffi(list_func_acti, ds_conv, ds_print, sderro, sddisc, &
                        'FIXE')
            ds_algorom%phase = 'HROM'
            call nmcrel(sderro, 'DIVE_FIXG', .false._1)
        end if
    end if
!
! - Set loop state
!
999 continue
    call nmevcv(sderro, list_func_acti, 'FIXE')
!
end subroutine
