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

subroutine nmevev(sddisc, nume_inst, valinc, sderro, ds_contact, &
                  loop_name)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmerge.h"
#include "asterfort/nmevel.h"
#include "asterfort/nmlecv.h"
#include "asterfort/nmltev.h"
#include "asterfort/nmcrel.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=24), intent(in) :: sderro
    character(len=19), intent(in) :: sddisc
    character(len=4), intent(in) :: loop_name
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: valinc(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! DETECTION DU PREMIER EVENEMENT DECLENCHE
!
! ----------------------------------------------------------------------
!
! NB: DES QU'UN EVENT-DRIVEN EST SATISFAIT, ON SORT
! ON NE CHERCHE PAS A VERIFIER LES AUTRES EVENEMENTS
!
! In  sddisc           : datastructure for time discretization TEMPORELLE
! In  ds_contact       : datastructure for contact management
! IN  NUMINS : NUMERO D'INSTANT
! IN  SDERRO : SD GESTION DES ERREURS
! IN  NOMBCL : NOM DE LA BOUCLE
!               'NEWT' - BOUCLE DE NEWTON
!               'FIXE' - BOUCLE DE POINT FIXE
! IN  VALINC : VARIABLE CHAPEAU INCREMENTS DES VARIABLES
!
! ----------------------------------------------------------------------
!
    aster_logical :: lsvimx, ldvres, linsta, lresmx
    aster_logical :: conver, lerror, lerrcv
!
! ----------------------------------------------------------------------
!
    call nmerge(sderro, 'SOLV_ITMX', lsvimx)
    call nmerge(sderro, 'DIVE_RESI', ldvres)
    call nmerge(sderro, 'RESI_MAXI', lresmx)
    call nmerge(sderro, 'CRIT_STAB', linsta)
!
! --- LA BOUCLE COURANTE A-T-ELLE CONVERGE ?
!
    call nmlecv(sderro, loop_name, conver)
!
! --- DECLENCHEMENT DE L'ERREUR A CONVERGENCE
!
    call nmltev(sderro, 'ERRC', loop_name, lerror)
    if (conver .and. lerror) then
        lerrcv = .true.
        call nmcrel(sderro, 'INTE_NPHY', ASTER_FALSE)
        call nmcrel(sderro, 'ERRE_NPHY', ASTER_TRUE)
    else
        lerrcv = .false.
    end if
!
! --- ERREUR IMMEDIATE AU NIVEAU DE LA BOUCLE COURANTE ?
!
    call nmltev(sderro, 'ERRI', loop_name, lerror)
!
! --- PREMIER EVENT DECLENCHE
!
    call nmevel(sddisc, nume_inst, valinc, loop_name, lsvimx, &
                ldvres, lresmx, linsta, lerrcv, lerror, &
                conver, ds_contact)
!
end subroutine
