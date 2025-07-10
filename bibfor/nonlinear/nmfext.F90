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
subroutine nmfext(eta, listFuncActi, veasse, cnfext, &
                  ds_contact_, sddyna_, nlDynaDamping_)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cfdisl.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ndasva.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmasfi.h"
#include "asterfort/nmasva.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
#include "asterfort/utmess.h"
!
    real(kind=8), intent(in) :: eta
    integer(kind=8), intent(in) :: listFuncActi(*)
    character(len=19) :: veasse(*)
    type(NL_DS_Contact), optional, intent(in) :: ds_contact_
    character(len=19), intent(in) :: cnfext
    character(len=19), optional, intent(in) :: sddyna_
    type(NLDYNA_DAMPING), optional, intent(in) :: nlDynaDamping_
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! RESULTANTE DES EFFORTS EXTERIEURS
!
! --------------------------------------------------------------------------------------------------
!
! In  listFuncActi     : list of active functionnalities
! IN  ETA    : COEFFICIENT DE PILOTAGE
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! In  ds_contact       : datastructure for contact management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! OUT CNFEXT : CHARGEMENT EXTERIEUR RESULTANT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: cnffdo, cnffpi, cnfvdo, cnvady, sddyna
    aster_logical :: lctcd, lunil, l_unil_pena
    real(kind=8) :: coeequ
    aster_logical :: ldyna, lallv, l_pilo
    integer(kind=8) :: ifdo, n
    type(NLDYNA_DAMPING) :: nlDynaDamping
    character(len=19) :: vect(20)
    real(kind=8) :: coef(20)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_44')
    end if

! - Active functionnalities
    ldyna = isfonc(listFuncActi, 'DYNAMIQUE')
    lctcd = isfonc(listFuncActi, 'CONT_DISCRET')
    lunil = isfonc(listFuncActi, 'LIAISON_UNILATER')
    lallv = isfonc(listFuncActi, 'CONT_ALL_VERIF')
    l_pilo = isfonc(listFuncActi, 'PILOTAGE')

! - Initializations
    sddyna = ' '
    if (present(sddyna_)) then
        sddyna = sddyna_
        nlDynaDamping = nlDynaDamping_
    end if
    ifdo = 0
    cnffdo = '&&CNCHAR.FFDO'
    cnffpi = '&&CNCHAR.FFPI'
    cnfvdo = '&&CNCHAR.FVDO'
    cnvady = '&&CNCHAR.FVDY'
    call vtzero(cnfext)
!
! --- FORCES DE CONTACT DISCRET
!
    if (lctcd .and. (.not. lallv)) then
        ifdo = ifdo+1
        coef(ifdo) = -1.d0
        vect(ifdo) = ds_contact_%cnctdc
    end if
!
! --- FORCES DE LIAISON_UNILATER
!
    if (lunil) then
!    On desactive pour l'instant en penalisation
        l_unil_pena = cfdisl(ds_contact_%sdcont_defi, 'UNIL_PENA')
        if (.not. l_unil_pena) then
            ifdo = ifdo+1
            coef(ifdo) = -1.d0
            vect(ifdo) = ds_contact_%cnunil
        end if
    end if
!
! - Get dead Neumann loads and multi-step dynamic schemes forces
!
    call nmasfi(listFuncActi, veasse, cnffdo, sddyna)
!
! - Get undead Neumann loads and multi-step dynamic schemes forces
!
    call nmasva(listFuncActi, veasse, cnfvdo, sddyna)
!
! - Get undead Neumann loads for dynamic
!
    if (ldyna) then
        coeequ = ndynre(sddyna, 'COEF_MPAS_EQUI_COUR')
        call ndasva(sddyna, nlDynaDamping, veasse, cnvady)
    end if
!
! --- CHARGEMENTS EXTERIEURS DONNEES
!
    ifdo = ifdo+1
    coef(ifdo) = 1.d0
    vect(ifdo) = cnffdo
    ifdo = ifdo+1
    coef(ifdo) = 1.d0
    vect(ifdo) = cnfvdo
!
! - Get dead Neumann loads (for PILOTAGE)
!
    if (l_pilo) then
        call nmchex(veasse, 'VEASSE', 'CNFEPI', cnffpi)
        ifdo = ifdo+1
        coef(ifdo) = eta
        vect(ifdo) = cnffpi
    end if
!
! --- TERMES DE RAPPEL DYNAMIQUE
!
    if (ldyna) then
        ifdo = ifdo+1
        coef(ifdo) = coeequ
        vect(ifdo) = cnvady
    end if
!
! --- VECTEUR RESULTANT
!
    if (ifdo .gt. 20) then
        ASSERT(.false.)
    end if
    do n = 1, ifdo
        call vtaxpy(coef(n), vect(n), cnfext)
    end do
!
! --- AFFICHAGE
!
    if (niv .ge. 2) then
        call nmdebg('VECT', cnfext, ifm)
    end if
!
    call jedema()
end subroutine
