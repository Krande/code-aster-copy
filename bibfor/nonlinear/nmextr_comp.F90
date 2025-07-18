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
subroutine nmextr_comp(field, field_disc, field_type, meshz, modelz, &
                       cara_elemz, ds_material, ds_constitutive, disp_curr, strx_curr, &
                       varc_curr, time, ligrelz)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/inical.h"
#include "asterfort/megeom.h"
#include "asterfort/mecara.h"
#include "asterfort/meharm.h"
#include "asterfort/mecact.h"
!
    character(len=19), intent(in) :: field
    character(len=24), intent(in) :: field_type
    character(len=4), intent(in) :: field_disc
    character(len=*), intent(in) :: modelz
    character(len=*), intent(in) :: meshz
    character(len=*), intent(in) :: cara_elemz
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=*), intent(in) :: disp_curr
    character(len=*), intent(in) :: strx_curr
    character(len=*), intent(in) :: varc_curr
    real(kind=8), intent(in) :: time
    character(len=*), optional, intent(in) :: ligrelz
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Field extraction datastructure
!
! Compute fields when not a default in nonlinear operator
!
! ONLY EPSI_ELGA !
!
! --------------------------------------------------------------------------------------------------
!
! In  field            : name of field
! In  field_disc       : localization of field (discretization: NOEU or ELGA)
! In  field_type       : type of field (name in results datastructure)
! In  model            : name of model
! In  mesh             : name of mesh
! In  cara_elem        : name of datastructure for elementary parameters (CARTE)
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  disp_curr        : current displacements
! In  varc_curr        : command variable for current time
! In  time             : current time
! In  strx_curr        : fibers information for current time
! In  ligrel           : current LIGREL (if not present: on all model)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=15)
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=24) :: chgeom, chcara(18), chharm, chtime
    integer(kind=8) :: n_harm
    character(len=19) :: ligrel
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(field_type .eq. 'EPSI_ELGA')
    ASSERT(field_disc .eq. 'ELGA')
    option = 'EPSI_ELGA'
    chtime = '&&NMEXTR_COMP.CHTIME'
    chharm = '&&NMEXTR_COMP.CHHARM'
    if (present(ligrelz)) then
        ligrel = ligrelz
    else
        ligrel = modelz(1:8)//'.MODELE'
    end if
    n_harm = 0
!
! - Time field
!
    call mecact('V', chtime, 'MAILLA', meshz, 'INST_R', &
                ncmp=1, nomcmp='INST', sr=time)
!
! - Geometry field
!
    call megeom(modelz, chgeom)
!
! - Elementary characteristics fields
!
    call mecara(cara_elemz, chcara)
!
! - Fourier field
!
    call meharm(modelz, n_harm, chharm)
!
! - Init fields
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! - Input fields
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PDEPLAR'
    lchin(2) = disp_curr(1:19)
    lpain(3) = 'PMATERC'
    lchin(3) = ds_material%mateco(1:19)
    lpain(4) = 'PINSTR'
    lchin(4) = chtime(1:19)
    lpain(5) = 'PVARCPR'
    lchin(5) = varc_curr(1:19)
    lpain(6) = 'PVARCRR'
    lchin(6) = ds_material%varc_refe(1:19)
    lpain(7) = 'PCACOQU'
    lchin(7) = chcara(7) (1:19)
    lpain(8) = 'PCOMPOR'
    lchin(8) = ds_constitutive%compor(1:19)
    lpain(9) = 'PCAGEPO'
    lchin(9) = chcara(5) (1:19)
    lpain(10) = 'PCAORIE'
    lchin(10) = chcara(1) (1:19)
    lpain(11) = 'PNBSP_I'
    lchin(11) = chcara(16) (1:19)
    lpain(12) = 'PFIBRES'
    lchin(12) = chcara(17) (1:19)
    lpain(13) = 'PCAMASS'
    lchin(13) = chcara(12) (1:19)
    lpain(14) = 'PHARMON'
    lchin(14) = chharm(1:19)
    lpain(15) = 'PSTRXMR'
    lchin(15) = strx_curr(1:19)
!
! - Output field
!
    lpaout(1) = 'PDEFOPG'
    lchout(1) = field
!
! - Computation
!
    call calcul('S', option, ligrel, nbin, lchin, &
                lpain, nbout, lchout, lpaout, 'V', &
                'OUI')
!
end subroutine
