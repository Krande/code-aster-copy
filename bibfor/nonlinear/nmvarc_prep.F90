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

subroutine nmvarc_prep(type_comp, model, cara_elem, mateco, varc_refe, &
                       compor, exis_temp, mxchin, nbin, lpain, &
                       lchin, mxchout, nbout, lpaout, lchout, &
                       sigm_prev, vari_prev, varc_prev, varc_curr, nume_harm)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/detrsd.h"
#include "asterfort/exixfe.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/nmvcex.h"
#include "asterfort/xajcin.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=1), intent(in) :: type_comp
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: mateco
    character(len=24), intent(in) :: varc_refe
    character(len=24), intent(in) :: cara_elem
    character(len=24), intent(in) :: compor
    aster_logical, intent(in) :: exis_temp
    integer(kind=8), intent(in) :: mxchin
    character(len=8), intent(inout) :: lpain(mxchin)
    character(len=19), intent(inout) :: lchin(mxchin)
    integer(kind=8), intent(out) :: nbin
    integer(kind=8), intent(in) :: mxchout
    character(len=8), intent(inout) :: lpaout(mxchout)
    character(len=19), intent(inout) :: lchout(mxchout)
    integer(kind=8), intent(out) :: nbout
    character(len=19), intent(in) :: sigm_prev
    character(len=19), intent(in) :: vari_prev
    character(len=19), intent(in) :: varc_prev
    character(len=19), intent(in) :: varc_curr
    integer(kind=8), intent(in) :: nume_harm
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Command variables - Fields preparation
!
! --------------------------------------------------------------------------------------------------
!
! In  type_comp      : type of computation
!                      '-' - Previous step
!                      '+' - Current step
! In  model          : name of model
! In  mateco         : name of coded material
! In  cara_elem      : name of elementary characteristics (field)
! In  varc_refe      : name of reference command variables vector
! In  compor         : name of comportment definition (field)
! In  exis_temp      : .true. if temperature variable command exists
! In  mxchin         : maximum number of input fields
! IO  lpain          : list of input parameters
! IO  lchin          : list of input fields
! IO  nbin           : number of input fields
! In  mxchout        : maximum number of output fields
! IO  lpaout         : list of output parameters
! IO  lchout         : list of output fields
! IO  nbout          : number of output fields
! In  sigm_prev      : stress at previous step
! In  vari_prev      : internal variables at previous step
! In  varc_prev      : command variables at previous step
! In  varc_curr      : command variables at current step
! In  nume_harm      : Fourier harmonic number
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lxfem
    integer(kind=8) :: iret
    character(len=19) :: chsith
    character(len=19) :: vrcmoi, vrcplu, time_curr, time_prev
    character(len=24) :: chgeom, chcara(18), chvref, ligrmo, chharm
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    call exixfe(model, iret)
    lxfem = iret .ne. 0
    chvref = '&&NMVCPR.VREF'
    ligrmo = model(1:8)//'.MODELE'
    chsith = '&&NMVCPR.CHSITH'
!
! - Init fields
!
    call inical(mxchin, lpain, lchin, mxchout, lpaout, &
                lchout)
!
! - Get command variables
!
    call nmvcex('TOUT', varc_refe, chvref)
    call nmvcex('INST', varc_curr, time_curr)
    call nmvcex('TOUT', varc_curr, vrcplu)
    call nmvcex('TOUT', varc_prev, vrcmoi)
    call nmvcex('INST', varc_prev, time_prev)
!
! - Geometry field
!
    call megeom(model, chgeom)
!
! - Elementary characteristics
!
    call mecara(cara_elem, chcara)
!
! - Input fields
!
    lpain(1) = 'PVARCRR'
    lchin(1) = chvref(1:19)
    lpain(2) = 'PGEOMER'
    lchin(2) = chgeom(1:19)
    lpain(3) = 'PMATERC'
    lchin(3) = mateco(1:19)
    lpain(4) = 'PCACOQU'
    lchin(4) = chcara(7) (1:19)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6) (1:19)
    lpain(6) = 'PCADISM'
    lchin(6) = chcara(3) (1:19)
    lpain(7) = 'PCAORIE'
    lchin(7) = chcara(1) (1:19)
    lpain(8) = 'PCAGNBA'
    lchin(8) = chcara(11) (1:19)
    lpain(9) = 'PCAARPO'
    lchin(9) = chcara(9) (1:19)
    lpain(10) = 'PCAMASS'
    lchin(10) = chcara(12) (1:19)
    lpain(11) = 'PCAGEPO'
    lchin(11) = chcara(5) (1:19)
    lpain(12) = 'PCONTMR'
    lchin(12) = sigm_prev
    lpain(13) = 'PVARIPR'
    lchin(13) = vari_prev
    lpain(14) = 'PCOMPOR'
    lchin(14) = compor(1:19)
    lpain(15) = 'PNBSP_I'
    lchin(15) = chcara(1) (1:8)//'.CANBSP'
    lpain(16) = 'PFIBRES'
    lchin(16) = chcara(1) (1:8)//'.CAFIBR'
    lpain(17) = 'PCINFDI'
    lchin(17) = chcara(15) (1:19)
    lpain(18) = 'PCADISK'
    lchin(18) = chcara(2) (1:19)
    call meharm(model, nume_harm, chharm)
    lpain(19) = 'PHARMON'
    lchin(19) = chharm(1:19)
    nbin = 19
!
! - Computation of elementary vectors - Previous
    if (type_comp .eq. '-') then
        nbin = nbin+1
        lpain(nbin) = 'PINSTR'
        lchin(nbin) = time_prev
        nbin = nbin+1
        lpain(nbin) = 'PVARCPR'
        lchin(nbin) = vrcmoi
! - Computation of elementary vectors - Current
    else if (type_comp .eq. '+') then
        nbin = nbin+1
        lpain(nbin) = 'PINSTR'
        lchin(nbin) = time_curr
        nbin = nbin+1
        lpain(nbin) = 'PVARCPR'
        lchin(nbin) = vrcplu
        nbin = nbin+1
        lpain(nbin) = 'PVARCMR'
        lchin(nbin) = vrcmoi
    end if
!
! - XFEM input fields
!
    if (lxfem .and. exis_temp) then
        call xajcin(model, 'CHAR_MECA_TEMP_R', mxchin, lchin, lpain, nbin)
    end if
!
! - Output fields
!
    lpaout(1) = 'PVECTUR'
    nbout = 1
!
! - XFEM output fields
!
    if (lxfem .and. exis_temp) then
        call detrsd('CHAM_ELEM', chsith)
        call alchml(ligrmo, 'SIEF_ELGA', 'PCONTRR', 'V', chsith, &
                    iret, ' ')
        lpaout(2) = 'PCONTRT'
        lchout(2) = chsith
        nbout = nbout+1
    end if
!
    call jedema()
end subroutine
