! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
subroutine nmrefe(model, compor, materCode, caraElem, nume_dof, &
                  ds_conv, hatValinc, hatVeelem, hatVeasse)
!
    use NonLin_Datastructure_type
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assmiv.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exixfe.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/nmchex.h"
#include "asterfort/reajre.h"
#include "asterfort/vemare.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: compor
    character(len=24), intent(in) :: materCode
    character(len=24), intent(in) :: caraElem
    character(len=24), intent(in) :: nume_dof
    type(NL_DS_Conv), intent(in) :: ds_conv
    character(len=19), intent(in) :: hatValinc(*)
    character(len=19), intent(in) :: hatVeelem(*)
    character(len=19), intent(in) :: hatVeasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Computation
!
! Compute reference vector for RESI_REFE_RELA
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  compor           : name of comportment definition (field)
! In  materCode        : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  nume_dof         : name of numbering (NUME_DDL)
! In  ds_conv          : datastructure for convergence management
! In  hatValinc        : hat variable for algorithm fields
! In  hatVeelem        : hat variable for elementary vectors
! In  hatVeasse        : hat variable for vectors
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    character(len=16), parameter :: option = 'REFE_FORC_NODA'
    character(len=19) :: vectElem, vectAsse, dispPrev
    character(len=19) :: modelLigrel
    character(len=19), parameter :: resuElem = '&&NMREFE.VEREFE'
    character(len=19), parameter :: chrefe = '&&NMREFE.SIGERE'
    character(len=24) :: chgeom
    integer(kind=8) :: i_refe, nb_refe, nbFieldIn, ier
    character(len=8), pointer :: list_cmp(:) => null()
    real(kind=8), pointer :: list_vale(:) => null()
    aster_logical :: lXFEM
!
! --------------------------------------------------------------------------------------------------
!

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call exixfe(model, ier)
    lXFEM = ier .ne. 0

! - Get names of fields
    call nmchex(hatValinc, 'VALINC', 'DEPMOI', dispPrev)
    call nmchex(hatVeelem, 'VEELEM', 'CNREFE', vectElem)
    call nmchex(hatVeasse, 'VEASSE', 'CNREFE', vectAsse)

! - Get parameters from convergence datastructure
    nb_refe = ds_conv%nb_refe
    AS_ALLOCATE(vk8=list_cmp, size=nb_refe)
    AS_ALLOCATE(vr=list_vale, size=nb_refe)
    do i_refe = 1, nb_refe
        list_cmp(i_refe) = ds_conv%list_refe(i_refe)%cmp_name
        list_vale(i_refe) = ds_conv%list_refe(i_refe)%user_para
    end do

! - Create field for reference values
    call mecact('V', chrefe, 'MODELE', modelLigrel, 'PREC_R', &
                ncmp=nb_refe, lnomcmp=list_cmp, vr=list_vale)

! - Get geometry field
    call megeom(model, chgeom)

! - Preparation of VECT_ELEM
    call detrsd('VECT_ELEM', vectElem)
    call vemare('V', vectElem, model)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PREFCO'
    lchin(2) = chrefe
    lpain(3) = 'PCOMPOR'
    lchin(3) = compor(1:19)
    lpain(4) = 'PMATERC'
    lchin(4) = materCode(1:19)
    lpain(5) = 'PDEPLMR'
    lchin(5) = dispPrev
    nbFieldIn = 5

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add XFEM fields
    if (lXFEM) then
        call xajcin(model, 'REFE_FORC_NODA', nbFieldInMax, lchin, lpain, &
                    nbFieldIn)
    end if

! - Add HHO field
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Set output field
    lpaout(1) = 'PVECTUR'
    lchout(1) = resuElem

! - Computation
    call calcul('S', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

! - Copying output field
    call reajre(vectElem, lchout(1), 'V')

! - Assembly
    call assmiv('V', vectAsse, 1, vectElem, [1.d0], &
                nume_dof, 1)
!
    AS_DEALLOCATE(vk8=list_cmp)
    AS_DEALLOCATE(vr=list_vale)
!
end subroutine
