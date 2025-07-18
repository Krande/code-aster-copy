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

subroutine utmess_core(typ, idmess, nk, valk, ni, &
                       vali, nr, valr, nexcep, fname)
! person_in_charge: mathieu.courtois at edf.fr
!
    use message_module, only: Message, init_message, free_message
    use superv_module, only: superv_after
    use parameters_module, only: ST_OK
    implicit none
#include "asterc/getres.h"
#include "asterc/jdcget.h"
#include "asterc/isjvup.h"
#include "asterc/uexcep.h"
#include "asterc/utprin.h"
#include "asterf_types.h"
#include "asterfort/asmpi_warn.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jefini.h"
#include "asterfort/jemarq.h"
#include "asterfort/jevema.h"
#include "asterfort/lxlgut.h"
#include "asterfort/onerrf.h"
#include "asterfort/ststat.h"
#include "asterfort/trabck.h"
    character(len=*), intent(in) :: typ
    character(len=*), intent(in) :: idmess
    integer(kind=8), intent(in) :: nk
    character(len=*), intent(in) :: valk(*)
    integer(kind=8), intent(in) :: ni
    integer(kind=8), intent(in) :: vali(*)
    integer(kind=8), intent(in) :: nr
    real(kind=8), intent(in) :: valr(*)
    integer(kind=8), intent(in) :: nexcep
    character(len=*), intent(in) :: fname
!
    integer(kind=8), save :: recurs
    character(len=24) :: msgId
    character(len=16) :: compex
    character(len=8) :: nomres, k8b
    character(len=2) :: typm
    aster_logical :: lerror, lvalid, labort, suite, lstop, lerrm, ltrb
    integer(kind=8) :: lout, idf, i, lc, imaap, icode
    integer(kind=8) :: numex
!
    aster_logical, save :: isFirst = ASTER_TRUE
    type(Message), save :: firstMsg
    type(Message) :: excMsg
!
!
!     TYPES DE MESSAGES :
!     ERREURS :
!       F : ERREUR AVEC DESTRUCTION DU CONCEPT PRODUIT PAR LA COMMANDE
!       S : ERREUR AVEC VALIDATION DU CONCEPT, EXCEPTION
!       Z : LEVEE D'EXCEPTION PARTICULIERE, COMME 'S'
!       M : ERREUR SUIVIE DE MPI_ABORT, NE PAS LEVER D'EXCEPTION --> 'F'
!     MESSAGES :
!       E : SIMPLE MESSAGE D'ERREUR QUI SERA SUIVI D'UNE ERREUR 'F'
!       D : COMME 'E' MAIS AFFICHE AVEC 'F' POUR ASSURER UN 'D'IAGNOSTIC
!       I : INFORMATION
!       A : ALARME
!
!     LE TRACEBACK INTEL, SI DISPO, EST AFFICHE EN CAS D'ERREUR OU
!     EXCEPTION DVP_NNN, OU ERREUR 'D' CAR SUIVIE DE MPI_ABORT
    msgId = idmess
    typm = typ
    idf = index('EFIMASZD', typm(1:1))
    if (idf .eq. 0) then
        idf = 2
    end if
!
    if (idf .eq. 5) then
        icode = jdcget("WarningAsError")
        if (icode .ne. 0) then
            idf = 2
            typm(1:1) = 'F'
            call utprin('E', 0, 'CATAMESS_3', 1, [msgId], &
                        0, vali, 0, valr, fname)
        end if
    end if
!
!     --- COMPORTEMENT EN CAS D'ERREUR
    lstop = .false.
    call onerrf(' ', compex, lout)
!
    lerrm = idf .eq. 4
    if (lerrm) then
        idf = 2
        typm(1:1) = 'F'
!       L'EXCEPTION A-T-ELLE DEJA ETE LEVEE ?
        if (recurs .ne. 0) then
!         L'EXCEPTION A DEJA ETE LEVEE
            recurs = 0
        else
            lerrm = .false.
        end if
    end if
!
    lerror = idf .eq. 2 .or. idf .eq. 6 .or. idf .eq. 7
!     DOIT-ON VALIDER LE CONCEPT ?
    lvalid = (idf .eq. 6 .or. idf .eq. 7) .or. &
             (idf .eq. 2 .and. compex(1:lout) .eq. 'EXCEPTION+VALID')
!     DOIT-ON S'ARRETER BRUTALEMENT (POUR DEBUG) ?
    labort = idf .eq. 2 .and. compex(1:lout) .eq. 'ABORT'
!     AFFICHER LE TRACEBACK SI DISPONIBLE
    ltrb = labort .or. (lerror .and. msgId(1:4) .eq. 'DVP_') .or. idf .eq. 8
!
    numex = nexcep
    if (numex .eq. 0 .or. (lerror .and. idf .ne. 7)) then
!     SI EXCEPTION, NEXCEP EST FIXE PAR COMMON VIA UTEXCP
!     SINON ON LEVE L'EXCEPTION DE BASE ASTER.ERROR
        numex = 1
    end if
!
    suite = .false.
    if (len(typm) .gt. 1) then
        if (typm(2:2) .eq. '+') suite = .true.
    end if
!
!   Keep the first message in memory because this is one that will be used
!   to raise the exception
    if (isFirst) then
        call init_message(firstMsg, typm, msgId, &
                          nk=nk, valk=valk, &
                          ni=ni, vali=vali, &
                          nr=nr, valr=valr, &
                          num_except=numex)
        isFirst = ASTER_FALSE
    end if
! --- SE PROTEGER DES APPELS RECURSIFS POUR LES MESSAGES D'ERREUR
    if (lerror) then
        if (recurs .eq. 1234567891) then
            call jefini('ERREUR')
        end if
!
        if (recurs .eq. 1234567890) then
            recurs = 1234567891
!          ON EST DEJA PASSE PAR UTMESG... SANS EN ETRE SORTI
            call utprin('F', 0, 'CATAMESS_55', 0, valk, &
                        0, vali, 0, valr, fname)
!          ON NE FAIT PLUS RIEN ET ON SORT DE LA ROUTINE
            goto 999
        end if
        recurs = 1234567890
    end if
!
    call jevema(imaap)
    if (imaap .ge. 200) call jefini('ERREUR')
    if (isjvup() .eq. 1) then
        call jemarq()
    end if
!
    call utprin(typm, numex, msgId, nk, valk, &
                ni, vali, nr, valr, fname)
!
!     --- REMONTEE D'ERREUR SI DISPO
    if (ltrb) then
        call trabck('Traceback printed by Intel compiler', int(-1, 4))
    end if
! --- EN CAS DE MESSAGE AVEC SUITE, PAS D'ARRET, PAS D'EXCEPTION
    if (.not. suite) then
!
!     -- ABORT SUR ERREUR <F> "ORDINAIRE"
        if (labort) then
!           AVERTIR LE PROC #0 QU'ON A RENCONTRE UN PROBLEME !
            call asmpi_warn(0)
!
            call jefini('ERREUR')
!
!     -- LEVEE D'UNE EXCEPTION
        else if (lerror) then
!
!        -- QUELLE EXCEPTION ?
!           SI EXCEPTION, NEXCEP EST FIXE PAR COMMON VIA UTEXCP
!           IL A ETE COPIE DANS NUMEX POUR NE PAS ETRE MODIFIE SI
!           DES APPELS SONT IMBRIQUES
            if (idf .ne. 7) then
!           SINON ON LEVE L'EXCEPTION DE BASE ASTER.ERROR
                numex = 1
            end if
!
            if (isjvup() .eq. 1) then
                call superv_after(exception=.true.)
            end if

!           NOM DU CONCEPT COURANT
            call getres(nomres, k8b, k8b)
            if (nomres .ne. ' ') then
!             LE CONCEPT EST REPUTE VALIDE :
!               - SI ERREUR <S> OU EXCEPTION
!               - SI ERREUR <F> MAIS LA COMMANDE A DIT "EXCEPTION+VALID"
                if (lvalid) then
                    call utprin('I', 0, 'CATAMESS_70', 1, nomres, &
                                0, vali, 0, valr, fname)
!
!             SINON LE CONCEPT COURANT EST DETRUIT
                else
                    call utprin('I', 0, 'CATAMESS_69', 1, nomres, &
                                0, vali, 0, valr, fname)
                    lc = lxlgut(nomres)
                    if (lc .gt. 0) then
                        call jedetc(' ', nomres(1:lc), 1)
                    end if
                end if
            end if
!
            if (isjvup() .eq. 1) then
!
!             REMONTER LES N JEDEMA COURT-CIRCUITES
                call jevema(imaap)
                do i = imaap, 1, -1
                    call jedema()
                end do
!
            end if
!
!           AVERTIR LE PROC #0 QU'ON A RENCONTRE UN PROBLEME !
            excMsg = firstMsg
            call asmpi_warn(1)
!
!           ON REMONTE UNE EXCEPTION AU LIEU DE FERMER LES BASES
            if (lerror) recurs = 0
            lstop = .true.
            if (.not. lerrm) then
!               raise the exception with the first msg id & reinit id
                isFirst = ASTER_TRUE
                call ststat(ST_OK)
                call superv_after(exception=.true.)
                call uexcep(numex, excMsg%id, excMsg%nk, excMsg%valk, excMsg%ni, &
                            excMsg%vali, excMsg%nr, excMsg%valr)
                call free_message(excMsg)
                call free_message(firstMsg)
            end if
        else
!           info/warning, reinit id
            if (firstMsg%typ .ne. 'F' .and. recurs .eq. 0) then
                isFirst = ASTER_TRUE
                call free_message(firstMsg)
            end if
        end if
!
    end if
!
    if (lerror) recurs = 0
999 continue
    if (isjvup() .eq. 1 .and. .not. lstop) then
        call jedema()
    end if
end subroutine
