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
subroutine nmchsv(fonact, veasse, sddyna, ds_system, ds_contact)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/ndynkk.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
!
    integer(kind=8) :: fonact(*)
    character(len=19) :: sddyna
    character(len=19) :: veasse(*)
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_Contact), intent(in) :: ds_contact
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (DYNAMIQUE - ALGORITHME)
!
! SAUVEGARDE DES VECTEURS SECONDS MEMBRES ET EFFORTS INTERIEURS
!
! ----------------------------------------------------------------------
!
!
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  SDDYNA : SD DYNAMIQUE
! In  ds_contact       : datastructure for contact management
! In  ds_system        : datastructure for non-linear system management
!
!
    character(len=19) :: olfedo, olfsdo, oldido, oldidi, olfint
    character(len=19) :: olondp, olcine, olviss, olhyst, olsstf
    character(len=19) :: cnfedo, cnfsdo, cndido, cndidi, cnfint
    character(len=19) :: cnondp, cncine, cnviss, cnhyst, cnsstf
    character(len=19) :: olsstr, cnsstr
    character(len=19) :: oleltc, oleltf
    aster_logical :: londe, ldidi, lviss, lsstf, l_macr
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- FONCTIONNALITES ACTIVEES
!
    londe = ndynlo(sddyna, 'ONDE_PLANE')
    lviss = ndynlo(sddyna, 'VECT_ISS')
    lsstf = isfonc(fonact, 'SOUS_STRUC')
    ldidi = isfonc(fonact, 'DIDI')
    l_macr = isfonc(fonact, 'MACR_ELEM_STAT')
!
! --- NOM DES CHAMPS PAS PRECEDENT
!
    call ndynkk(sddyna, 'OLDP_CNFEDO', olfedo)
    call ndynkk(sddyna, 'OLDP_CNFSDO', olfsdo)
    call ndynkk(sddyna, 'OLDP_CNDIDO', oldido)
    call ndynkk(sddyna, 'OLDP_CNDIDI', oldidi)
    call ndynkk(sddyna, 'OLDP_CNFINT', olfint)
    call ndynkk(sddyna, 'OLDP_CNONDP', olondp)
    call ndynkk(sddyna, 'OLDP_CNCINE', olcine)
    call ndynkk(sddyna, 'OLDP_CNVISS', olviss)
    call ndynkk(sddyna, 'OLDP_CNHYST', olhyst)
    call ndynkk(sddyna, 'OLDP_CNSSTF', olsstf)
    call ndynkk(sddyna, 'OLDP_CNSSTR', olsstr)
    call ndynkk(sddyna, 'OLDP_CNELTC', oleltc)
    call ndynkk(sddyna, 'OLDP_CNELTF', oleltf)
!
! --- NOM DES CHAMPS PAS COURANT
!
    call nmchex(veasse, 'VEASSE', 'CNFEDO', cnfedo)
    call nmchex(veasse, 'VEASSE', 'CNFSDO', cnfsdo)
    call nmchex(veasse, 'VEASSE', 'CNDIDO', cndido)
    call nmchex(veasse, 'VEASSE', 'CNDIDI', cndidi)
    cnfint = ds_system%cnfint
    call nmchex(veasse, 'VEASSE', 'CNONDP', cnondp)
    call nmchex(veasse, 'VEASSE', 'CNCINE', cncine)
    call nmchex(veasse, 'VEASSE', 'CNVISS', cnviss)
    call nmchex(veasse, 'VEASSE', 'CNHYST', cnhyst)
    call nmchex(veasse, 'VEASSE', 'CNSSTF', cnsstf)
    call nmchex(veasse, 'VEASSE', 'CNSSTR', cnsstr)
!
! --- RECOPIE DES CHAMPS
!
    call copisd('CHAMP_GD', 'V', cnfint, olfint)
    call copisd('CHAMP_GD', 'V', cnfedo, olfedo)
    call copisd('CHAMP_GD', 'V', cnfsdo, olfsdo)
    call copisd('CHAMP_GD', 'V', cndido, oldido)
    call copisd('CHAMP_GD', 'V', cncine, olcine)
    if (londe) then
        call copisd('CHAMP_GD', 'V', cnondp, olondp)
    end if
    if (ldidi) then
        call copisd('CHAMP_GD', 'V', cndidi, oldidi)
    end if
    if (lviss) then
        call copisd('CHAMP_GD', 'V', cnviss, olviss)
    end if
    call copisd('CHAMP_GD', 'V', cnhyst, olhyst)
    if (lsstf) then
        call copisd('CHAMP_GD', 'V', cnsstf, olsstf)
    end if
    if (l_macr) then
        call copisd('CHAMP_GD', 'V', cnsstr, olsstr)
    end if
    if (ds_contact%l_cneltc) then
        call copisd('CHAMP_GD', 'V', ds_contact%cneltc, oleltc)
    end if
    if (ds_contact%l_cneltf) then
        call copisd('CHAMP_GD', 'V', ds_contact%cneltf, oleltf)
    end if
!
    call jedema()
end subroutine
