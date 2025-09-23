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
!
subroutine nmevel(sddisc, vale, loop_name, lsvimx, &
                  ldvres, lresmx, linsta, lerrcv, lerror, &
                  conver, ds_contact_)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/eneven.h"
#include "asterfort/nmevdg.h"
#include "asterfort/nmevin.h"
#include "asterfort/utdidt.h"
#include "asterfort/getFailEvent.h"
!
    character(len=19), intent(in) :: vale(*)
    character(len=19), intent(in) :: sddisc
    character(len=4), intent(in) :: loop_name
    aster_logical, intent(in) :: lsvimx
    aster_logical, intent(in) :: ldvres
    aster_logical, intent(in) :: lresmx
    aster_logical, intent(in) :: linsta
    aster_logical, intent(in) :: lerrcv
    aster_logical, intent(in) :: lerror
    aster_logical, intent(in) :: conver
    type(NL_DS_Contact), optional, intent(in) :: ds_contact_
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
! IN  VALE   : INCREMENTS DES VARIABLES
!               OP0070: VARIABLE CHAPEAU
!               OP0033: TABLE
! IN  NOMBCL : NOM DE LA BOUCLE
!               'RESI' - RESIDUS D'EQUILIBRE
!               'NEWT' - BOUCLE DE NEWTON
!               'FIXE' - BOUCLE DE POINT FIXE
!               'INST' - BOUCLE SUR LES PAS DE TEMPS
! IN  LSVIMX : .TRUE. SI ITERATIONS MAX ATTEINT DANS SOLVEUR ITERATIF
! IN  LDVRES : .TRUE. SI DIVERGENCE DU RESIDU
! IN  LRESMX : .TRUE. SI DIVERGENCE DU RESIDU (TROP GRAND)
! IN  LINSTA : .TRUE. SI INSTABILITE DETECTEE
! IN  LERRCV : .TRUE. SI ERREUR A CONVERGENCE DECLENCHEE
! IN  LERROR : .TRUE. SI ERREUR DECLENCHEE
! IN  CONVER : .TRUE. SI BOUCLE COURANTE A CONVERGE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbFail, iFail, iFailActi
    integer(kind=8) :: eventType
!
! ----------------------------------------------------------------------
!
    iFailActi = 0
!
! --- NOMBRE D'EVENT-DRIVEN : NECHEC
!
    call utdidt('L', sddisc, 'LIST', 'NECHEC', vali_=nbFail)
!
! --- DETECTION DU _PREMIER_ EVENEMENT DECLENCHE
! --- DES QU'UN EVENT-DRIVEN EST SATISFAIT, ON SORT
! --- ON NE CHERCHE PAS A VERIFIER LES AUTRES EVENT
!
    do iFail = 1, nbFail
! ----- Get event type
        call getFailEvent(sddisc, iFail, eventType)

! ----- PAR DEFAUT: EVENEMENT NON ACTIVE
        call eneven(sddisc, iFail, ASTER_FALSE)
        if (eventType .eq. FAIL_EVT_ERROR) then
            if (lsvimx .or. lerrcv .or. lerror) then
                iFailActi = iFail
                goto 99
            end if
        else if (eventType .eq. FAIL_EVT_DIVE_RESI) then
            if (ldvres) then
                iFailActi = iFail
                if (iFailActi .ne. 0) then
                    goto 99
                end if
            end if
        else if (eventType .eq. FAIL_EVT_RESI_MAXI) then
            if (lresmx) then
                iFailActi = iFail
                if (iFailActi .ne. 0) then
                    goto 99
                end if
            end if
        else if (eventType .eq. FAIL_EVT_INCR_QUANT) then
            if (conver) then
                call nmevdg(sddisc, vale, iFail, iFailActi)
                if (iFailActi .ne. 0) then
                    goto 99
                end if
            end if
        else if (eventType .eq. FAIL_EVT_INTERPENE) then
            if (loop_name .eq. 'INST') then
                call nmevin(sddisc, ds_contact_, iFail, iFailActi)
                if (iFailActi .ne. 0) then
                    goto 99
                end if
            end if
        else if (eventType .eq. FAIL_EVT_INSTABILITY) then
            if (linsta) then
                iFailActi = iFail
            end if
            if (iFailActi .ne. 0) then
                goto 99
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end do
!
99  continue

! - DECLENCHEMENT DE L'EVENEMENT
    if (iFailActi .ne. 0) then
        call eneven(sddisc, iFailActi, ASTER_TRUE)
    end if
!
end subroutine
