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
subroutine nonlinDSDynamicInit(hval_incr, sddyna, nlDynaDamping, ds_constitutive)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/ndynlo.h"
#include "asterfort/mginfo.h"
#include "asterfort/nmchex.h"
#include "asterfort/dismoi.h"
#include "asterfort/iden_nume.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndynkk.h"
#include "asterfort/ndynin.h"
#include "asterfort/trmult.h"
#include "asterfort/zerlag.h"
!
    character(len=19), intent(in) :: hval_incr(*)
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Dynamic management
!
! Initializations for dynamic
!
! --------------------------------------------------------------------------------------------------
!
! In  hval_incr        : hat-variable for incremental values fields
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! in  ds_constitutive  : datastructure for constitutive laws management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: nbMode, nbEqua
    integer(kind=8) :: iExci, nbExci
    integer(kind=8) :: jvMultSuppProj, jvNumeDofDEEQ
    aster_logical :: lDampModal, lMultiSupport
    character(len=19) :: pfcn1, pfcn2, dispPrev, multSuppProj
    character(len=8) :: mesh
    character(len=14) :: numeDof
    character(len=24) :: numeDofDEEQ, matrix, multSuppMode, dampMode
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE15_3')
    end if

! - Active functionnalities
    lDampModal = nlDynaDamping%lDampModal
    lMultiSupport = ndynlo(sddyna, 'MULTI_APPUI')

! - Modal damping: check numbering of equation
    if (lDampModal) then
        dampMode = nlDynaDamping%modalDamping%dampMode
        call mginfo(dampMode, numeDof, nbMode, nbEqua)
        call nmchex(hval_incr, 'VALINC', 'DEPMOI', dispPrev)
        call dismoi('NUME_EQUA', dispPrev, 'CHAM_NO', repk=pfcn1)
        call dismoi('NUME_EQUA', numeDof, 'NUME_DDL', repk=pfcn2)
        if (.not. iden_nume(pfcn1, pfcn2)) then
            call utmess('F', 'DYNAMIQUE_54')
        end if
    end if

! - Multi-support analysis
    if (lMultiSupport) then
! ----- Get parameters about modes for multi support
        call ndynkk(sddyna, 'multSuppMode', multSuppMode)
        call dismoi('REF_RIGI_PREM', multSuppMode, 'RESU_DYNA', repk=matrix)
        call dismoi('NB_EQUA', matrix, 'MATR_ASSE', repi=nbEqua)
        call dismoi('NOM_MAILLA', matrix, 'MATR_ASSE', repk=mesh)
        call dismoi('NOM_NUME_DDL', matrix, 'MATR_ASSE', repk=numeDof)
        numeDofDEEQ = numeDof//'.NUME.DEEQ'
        call jeveuo(numeDofDEEQ, 'L', jvNumeDofDEEQ)
! ----- Prepare loads for multi-support by projection
        nbExci = ndynin(sddyna, 'NBRE_EXCIT')
        call ndynkk(sddyna, 'MUAP_MAPSID', multSuppProj)
        call jeveuo(multSuppProj, 'E', jvMultSuppProj)
        do iExci = 1, nbExci
            call trmult(multSuppMode, iExci, mesh, nbEqua, jvNumeDofDEEQ, &
                        zr(jvMultSuppProj+(iExci-1)*nbEqua), numeDof)
            call zerlag(nbEqua, zi(jvNumeDofDEEQ), vectr=zr(jvMultSuppProj+(iExci-1)*nbEqua))
        end do
! ----- Check linear/DIS_* only !
        if (.not. ds_constitutive%lLinear .and. .not. ds_constitutive%lDisCtc) then
            call utmess('F', 'DYNAMIQUE1_1')
        end if
    end if
!
end subroutine
