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
subroutine merimp(l_xfem, l_dyna, &
                  model, cara_elem, sddyna, iter_newt, &
                  ds_constitutive, ds_material, &
                  hval_incr, hval_algo, caco3d, &
                  mxchin, lpain, lchin, nbin)
!
    use NonLin_Datastructure_type
    use HHO_precalc_module, only: hhoAddInputField
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cesvar.h"
#include "asterfort/copisd.h"
#include "asterfort/exisd.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/ndynkk.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmvcex.h"
#include "asterfort/vtzero.h"
#include "asterfort/xajcin.h"
!
    aster_logical, intent(in) :: l_xfem, l_dyna
    character(len=24), intent(in) :: model, cara_elem
    character(len=19), intent(in) :: sddyna
    integer(kind=8), intent(in) :: iter_newt
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=24), intent(in) :: caco3d
    integer(kind=8), intent(in) :: mxchin
    character(len=8), intent(inout) :: lpain(mxchin)
    character(len=19), intent(inout) :: lchin(mxchin)
    integer(kind=8), intent(out) :: nbin
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Computation of rigidity matrix and internal forces - Input fields
!
! --------------------------------------------------------------------------------------------------
!
! In  l_xfem           : flag for XFEM elements
! In  l_dyna           : flag for dynamic
! In  l_hho            : flag for HHO elements
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  sddyna           : datastructure for dynamic
! In  iter_newt        : index of current Newton iteration
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_material      : datastructure for material parameters
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hhoField         : datastructure for HHO method
! In  caco3d           : name of field for COQUE_3D (field of normals)
! In  mxchin           : maximum number of input fields
! IO  lpain            : list of input parameters
! IO  lchin            : list of input fields
! Out nbin             : number of input fields
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iret
    character(len=24) :: chgeom, chcara(18), chiter
    character(len=19) :: stadyn, depent, vitent
    character(len=16) :: option
    character(len=19) :: disp_prev, sigm_prev, vari_prev, varc_prev, stru_prev
    character(len=19) :: disp_curr, sigm_curr, vari_curr, varc_curr, stru_curr
    character(len=19) :: time_prev, vrcmoi
    character(len=19) :: time_curr, vrcplu
    character(len=24) :: vrcref
    character(len=24) :: vari_iter, stru_iter
    character(len=19) :: depkm1, vitkm1, acckm1
    character(len=19) :: vite_curr, acce_curr, vite_prev, acce_prev
    character(len=19) :: romkm1, romk
    character(len=24) :: ligrmo
    character(len=19) :: disp_iter, disp_cumu_inst
!
! --------------------------------------------------------------------------------------------------
!
    ligrmo = model(1:8)//'.MODELE'
    option = 'FULL_MECA'
    chiter = '&&MERIMO.CH_ITERAT'
    vari_iter = '&&MERIMO.VARMOJ'
    stru_iter = '&&MERIMO.STRMOJ'
    nbin = 0
!
! - Get fields from hat-variables - Begin of time step
!
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', disp_prev)
    call nmchex(hval_incr, 'VALINC', 'VITMOI', vite_prev)
    call nmchex(hval_incr, 'VALINC', 'ACCMOI', acce_prev)
    call nmchex(hval_incr, 'VALINC', 'SIGMOI', sigm_prev)
    call nmchex(hval_incr, 'VALINC', 'VARMOI', vari_prev)
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varc_prev)
    call nmchex(hval_incr, 'VALINC', 'STRMOI', stru_prev)
!
! - Get fields from hat-variables - End of time step
!
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', disp_curr)
    call nmchex(hval_incr, 'VALINC', 'VITPLU', vite_curr)
    call nmchex(hval_incr, 'VALINC', 'ACCPLU', acce_curr)
    call nmchex(hval_incr, 'VALINC', 'SIGPLU', sigm_curr)
    call nmchex(hval_incr, 'VALINC', 'VARPLU', vari_curr)
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varc_curr)
    call nmchex(hval_incr, 'VALINC', 'STRPLU', stru_curr)
!
    call nmchex(hval_incr, 'VALINC', 'DEPKM1', depkm1)
    call nmchex(hval_incr, 'VALINC', 'VITKM1', vitkm1)
    call nmchex(hval_incr, 'VALINC', 'ACCKM1', acckm1)
    call nmchex(hval_incr, 'VALINC', 'ROMKM1', romkm1)
    call nmchex(hval_incr, 'VALINC', 'ROMK  ', romk)
!
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', disp_cumu_inst)
    call nmchex(hval_algo, 'SOLALG', 'DDEPLA', disp_iter)
!
! - Dynamic fields
!
    if (l_dyna) then
        call ndynkk(sddyna, 'DEPENT', depent)
        call ndynkk(sddyna, 'VITENT', vitent)
        call ndynkk(sddyna, 'STADYN', stadyn)
    end if
!
! - Get external state variables
!
    call nmvcex('TOUT', varc_prev, vrcmoi)
    call nmvcex('INST', varc_prev, time_prev)
    call nmvcex('TOUT', varc_curr, vrcplu)
    call nmvcex('INST', varc_curr, time_curr)
    call nmvcex('TOUT', ds_material%varc_refe, vrcref)
!
! - Get internal variables from previous iteration
!
    call exisd('CHAMP_GD', vari_curr(1:19), iret)
    if (iret .ne. 0) then
        call copisd('CHAMP_GD', 'V', vari_curr(1:19), vari_iter(1:19))
    else
        call copisd('CHAMP_GD', 'V', vari_prev(1:19), vari_iter(1:19))
    end if
!
! - Get structural variables from previous iteration
!
    call exisd('CHAMP_GD', stru_prev(1:19), iret)
    if (iret .ne. 0 .and. iter_newt .lt. 2) then
        call copisd('CHAMP_GD', 'V', stru_prev(1:19), stru_iter(1:19))
        call vtzero(stru_iter(1:19), 'CHAM_ELEM')
    elseif (iter_newt .ge. 2) then
        call copisd('CHAMP_GD', 'V', stru_curr(1:19), stru_iter(1:19))
    end if
!
! - Extend elementary field for internal variables
!
    call exisd('CHAM_ELEM_S', ds_constitutive%compor, iret)
    if (iret .eq. 0) then
        call cesvar(cara_elem, ds_constitutive%compor, ligrmo, ds_constitutive%compor)
    end if
    call copisd('CHAM_ELEM_S', 'V', ds_constitutive%compor, vari_curr)
    call copisd('CHAM_ELEM_S', 'V', ds_constitutive%compor, sigm_curr)
!
! - Geometry field
!
    call megeom(model, chgeom)
!
! - Elementary characteristics
!
    call mecara(cara_elem, chcara)
!
! - Field for iteration number
!
    call mecact('V', chiter, 'MODELE', ligrmo, 'NEUT_I', &
                ncmp=1, nomcmp='X1', si=iter_newt)
!
! - Input fields
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = ds_material%mateco(1:19)
    lpain(3) = 'PCONTMR'
    lchin(3) = sigm_prev(1:19)
    lpain(4) = 'PVARIMR'
    lchin(4) = vari_prev(1:19)
    lpain(5) = 'PCOMPOR'
    lchin(5) = ds_constitutive%compor(1:19)
    lpain(6) = 'PDEPLMR'
    lchin(6) = disp_prev(1:19)
    lpain(7) = 'PDEPLPR'
    lchin(7) = disp_cumu_inst(1:19)
    lpain(8) = 'PCACABL'
    lchin(8) = chcara(10) (1:19)
    lpain(9) = 'PINSTMR'
    lchin(9) = time_prev(1:19)
    lpain(10) = 'PINSTPR'
    lchin(10) = time_curr(1:19)
    lpain(11) = 'PCARCRI'
    lchin(11) = ds_constitutive%carcri(1:19)
    lpain(12) = 'PCAGNPO'
    lchin(12) = chcara(6) (1:19)
    lpain(13) = 'PCAORIE'
    lchin(13) = chcara(1) (1:19)
    lpain(14) = 'PCADISK'
    lchin(14) = chcara(2) (1:19)
    lpain(15) = 'PCACOQU'
    lchin(15) = chcara(7) (1:19)
    lpain(16) = 'PITERAT'
    lchin(16) = chiter(1:19)
    lpain(17) = 'PDDEPLA'
    lchin(17) = disp_iter(1:19)
    lpain(18) = 'PDEPKM1'
    lchin(18) = depkm1(1:19)
    lpain(19) = 'PVITKM1'
    lchin(19) = vitkm1(1:19)
    lpain(20) = 'PACCKM1'
    lchin(20) = acckm1(1:19)
    lpain(21) = 'PROMKM1'
    lchin(21) = romkm1(1:19)
    lpain(22) = 'PROMK'
    lchin(22) = romk(1:19)
    lpain(23) = 'PVARIMP'
    lchin(23) = vari_iter(1:19)
    lpain(24) = 'PCAGNBA'
    lchin(24) = chcara(11) (1:19)
    lpain(25) = 'PCAMASS'
    lchin(25) = chcara(12) (1:19)
    lpain(26) = 'PCAGEPO'
    lchin(26) = chcara(5) (1:19)
    lpain(27) = 'PVARCMR'
    lchin(27) = vrcmoi(1:19)
    lpain(28) = 'PVARCPR'
    lchin(28) = vrcplu(1:19)
    lpain(29) = 'PNBSP_I'
    lchin(29) = chcara(16) (1:19)
    lpain(30) = 'PFIBRES'
    lchin(30) = chcara(17) (1:19)
    lpain(31) = 'PCINFDI'
    lchin(31) = chcara(15) (1:19)
    lpain(32) = 'PVARCRR'
    lchin(32) = vrcref(1:19)
    lpain(33) = 'PCACO3D'
    lchin(33) = caco3d(1:19)
    lpain(34) = 'PCAARPO'
    lchin(34) = chcara(9) (1:19)
    lpain(35) = 'PSTRXMR'
    lchin(35) = stru_prev(1:19)
    lpain(36) = 'PSTRXMP'
    lchin(36) = stru_iter(1:19)
    lpain(37) = 'PMULCOM'
    lchin(37) = ds_constitutive%mult_comp(1:19)
    nbin = 37
!
! - XFEM fields
!
    if (l_xfem) then
        call xajcin(model, option, mxchin, lchin, lpain, nbin)
    end if
!
! - Dynamic
!
    if (l_dyna) then
        nbin = nbin+1
        lpain(nbin) = 'PDEPENT'
        lchin(nbin) = depent(1:19)
        nbin = nbin+1
        lpain(nbin) = 'PVITENT'
        lchin(nbin) = vitent(1:19)
        nbin = nbin+1
        lpain(nbin) = 'PSTADYN'
        lchin(nbin) = stadyn(1:19)
        nbin = nbin+1
        lpain(nbin) = 'PVITPLU'
        lchin(nbin) = vite_curr(1:19)
        nbin = nbin+1
        lpain(nbin) = 'PACCPLU'
        lchin(nbin) = acce_curr(1:19)
        nbin = nbin+1
        lpain(nbin) = 'PVITMOI'
        lchin(nbin) = vite_prev(1:19)
        nbin = nbin+1
        lpain(nbin) = 'PACCMOI'
        lchin(nbin) = acce_prev(1:19)
    end if
!
! - HHO
!
    call hhoAddInputField(model, mxchin, lchin, lpain, nbin)
!
end subroutine
