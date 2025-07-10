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

subroutine rsinch(nomsd, nomch, acces, rval, chextr, &
                  proldr, prolga, istop, base, prec, crit, ier)

    use searchlist_module, only: almostEqual
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/barych.h"
#include "asterfort/copisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsbary.h"
#include "asterfort/rsexch.h"
#include "asterfort/rslipa.h"
#include "asterfort/rsutro.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/rsinchpre.h"

!
    integer(kind=8), intent(in) :: istop
    real(kind=8), intent(in) :: rval
    character(len=*), intent(in) :: nomsd, nomch, acces, chextr, proldr, prolga
    character(len=*), intent(in) :: base
    real(kind=8), intent(in) :: prec
    character(len=8), intent(in) :: crit
    integer(kind=8), intent(out) :: ier

!      INTERPOLATION D'UN CHAMP_19 A PARTIR D'1 SD RESULTAT-COMPOSE
! ----------------------------------------------------------------------
! IN  : NOMSD  : NOM DE LA STRUCTURE "RESULTAT"
! IN  : NOMCH  : NOM SYMBOLIQUE DU CHAMP CHERCHE.
! IN  : ACCES  : NOM SYMBOLIQUE DE LA VARIABLE D'ACCES.
! IN  : RVAL   : VALEUR REEL DE LA VARIABLE D'ACCES.
! IN  : CHEXTR : NOM DU CHAMP A CREER. (S'IL EXISTE, ON LE DETRUIT).
! IN  : PROLDR : 'CONSTANT', 'LINEAIRE', OU 'EXCLU'
!                          (PROLONGEMENT VOULU A DROITE)
! IN  : PROLGA : 'CONSTANT', 'LINEAIRE', OU 'EXCLU'
!                          (PROLONGEMENT VOULU A GAUCHE)
! IN  : ISTOP  :  EN CAS D'ERREUR D'INTERPOLATION:
!                 0  --> N'ECRIT PAS DE MESSAGE , NE FAIT PAS STOP.
!                 1  --> ECRIT MESSAGES , NE FAIT PAS STOP.
!                 2  --> ECRIT MESSAGES , FAIT STOP.
! IN  : BASE   : BASE DU CHAMP CREE
!
! OUT : IER    : CODE_RETOUR :
!                LE CHAMP EST CALCULE:
!                00 --> LE CHAMP EST INTERPOLE ENTRE 2 VALEURS.
!                01 --> LE CHAMP EST PROLONGE A GAUCHE.
!                02 --> LE CHAMP EST PROLONGE A DROITE.
!
!                LE CHAMP N'EST PAS CALCULE:
!                10 --> IL N'EXISTE AUCUN CHAMP POUR L'INTERPOLATION.
!                11 --> LE PROLONGEMENT A GAUCHE INTERDIT.
!                12 --> SI PROLONGEMENT A DROITE INTERDIT.
!                20 --> LA VARIABLE D'ACCES EST ILLICITE.
! ----------------------------------------------------------------------
    real(kind=8) :: r1, r2, rbase
    real(kind=8) :: valr
    integer(kind=8) :: l1, l2
    character(len=1) :: stp, base2
    character(len=19) :: ch1, ch2
    character(len=8) :: prold2, prolg2
    character(len=19) :: noms2
    character(len=16) :: acce2, nomc2
    character(len=19) :: chext2
    character(len=24) :: valk(3)
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, i2, iaobj, iatach
    integer(kind=8) :: iadesc, ier1, ier2
    integer(kind=8) :: ip1, ip2, iposit, nbordr, nbvalid
    aster_logical, pointer :: lexi(:) => null()

!-----------------------------------------------------------------------
    call jemarq()
    acce2 = acces
    noms2 = nomsd
    nomc2 = nomch
    prold2 = proldr
    prolg2 = prolga
    chext2 = chextr
    base2 = base

    call rsinchpre(noms2, nomc2, acce2, ier)
    if (ier .eq. 0) then
!
        call rslipa(noms2, acces, '&&RSINCH.LIR8', iaobj, nbordr)
!
!     -- ON REPERE QUELS SONT LES CHAMPS EXISTANT REELLEMENT:
        AS_ALLOCATE(vl=lexi, size=nbordr)
        call jenonu(jexnom(noms2//'.DESC', nomc2), iadesc)
        call jeveuo(jexnum(noms2//'.TACH', iadesc), 'L', iatach)
        do i = 1, nbordr
            if (zk24(iatach-1+i) (1:1) .eq. ' ') then
                lexi(i) = .false.
            else
                lexi(i) = .true.
            end if
        end do
        nbvalid = count(lexi(1:nbordr))
!
        call rsbary(zr(iaobj), nbordr, ASTER_FALSE, lexi, rval, &
                    i1, i2, iposit, prec, crit)
        AS_DEALLOCATE(vl=lexi)
        ASSERT(iposit .ne. -2)
        call rsutro(nomsd, i1, ip1, ier1)
        call rsutro(nomsd, i2, ip2, ier2)
        ASSERT(ier1+ier2 .le. 0)
        rbase = zr(iaobj-1+i2)-zr(iaobj-1+i1)
!
        call rsexch(' ', nomsd, nomc2, ip1, ch1, l1)
        call rsexch(' ', nomsd, nomc2, ip2, ch2, l2)
        ASSERT(l1+l2 .le. 0)
!
        if (i1 .eq. i2) then
            ! --- RECOPIE DU CHAMP SI LES 2 POINTS IP1 ET IP2 ONT MEME ABSCISSE
            if (iposit .eq. 0) then
                ! CAS EXCEPTION OU ON DEMANDE UN INSTANT EXISTANT DANS LA LISTE
                call copisd('CHAMP_GD', base2, ch1(1:19), chext2(1:19))
                ier = 0
            else if ((iposit .eq. -1)) then
                ! CAS EXCEPTION DE PROLONGATION GAUCHE AVEC LISTE D'UN SEUL INSTANT
                ASSERT(nbvalid .eq. 1)
                if (prold2 .ne. "CONSTANT") then
                    ! SEULE L'OPTION CONSTANT A DU SENS
                    call utmess("A", "CALCULEL_28")
                end if
                call copisd('CHAMP_GD', base2, ch1(1:19), chext2(1:19))
                ier = 1
            else if ((iposit .eq. 1)) then
                ASSERT(nbvalid .eq. 1)
                ! CAS EXCEPTION DE PROLONGATION DROITE AVEC LISTE D'UN SEUL INSTANT
                if (prold2 .ne. "CONSTANT") then
                    ! SEULE L'OPTION CONSTANT A DU SENS
                    call utmess("A", "CALCULEL_28")
                end if
                call copisd('CHAMP_GD', base2, ch1(1:19), chext2(1:19))
                ier = 2
            else
                ASSERT(.false.)
            end if
        else
            !     -- INTERPOLATION VRAIE:
            !     -----------------------
            r1 = (zr(iaobj-1+i2)-rval)/rbase
            r2 = (rval-zr(iaobj-1+i1))/rbase

            if (iposit .eq. 0) then
                call barych(ch1, ch2, r1, r2, chext2, base2, nomsd)
                ier = 0
            else if (iposit .eq. -1) then
                !        -- PROLONGEMENT A GAUCHE:
                !        -------------------------
                if (prolg2(1:8) .eq. 'LINEAIRE') then
                    call barych(ch1, ch2, r1, r2, chext2, base2, nomsd)
                    ier = 1
                else if (prolg2(1:8) .eq. 'CONSTANT') then
                    call copisd('CHAMP_GD', base2, ch1(1:19), chext2(1:19))
                    ier = 1
                else
                    ier = 11
                end if

!        -- PROLONGEMENT A DROITE:
!        -------------------------
            else if (iposit .eq. 1) then
                if (prold2(1:8) .eq. 'LINEAIRE') then
                    call barych(ch1, ch2, r1, r2, chext2, base2, nomsd)
                    ier = 2
                else if (prold2(1:8) .eq. 'CONSTANT') then
                    call copisd('CHAMP_GD', base2, ch2(1:19), chext2(1:19))
                    ier = 2
                else
                    ier = 12
                end if
            else
                ASSERT(.false.)
            end if
        end if
        call jedetr('&&RSINCH.LIR8')
    end if
!
!     -- MESSAGES, ARRET?
!     -------------------
    if (istop .ne. 0) then
        if (istop .eq. 1) then
            stp = 'A'
        else if (istop .eq. 2) then
            stp = 'F'
        else
            ASSERT(.false.)
        end if
!
!
        if (ier .eq. 11) then
            call utmess(stp//'+', 'UTILITAI8_32')
        else if (ier .eq. 12) then
            call utmess(stp//'+', 'UTILITAI8_33')
        else if (ier .eq. 10) then
            valk(1) = nomc2
            call utmess(stp//'+', 'UTILITAI8_34', sk=valk(1))
        else if (ier .eq. 20) then
            valk(1) = acce2
            call utmess(stp//'+', 'UTILITAI8_35', sk=valk(1))
        else if (ier .eq. 21) then
            valk(1) = nomc2
            call utmess(stp//'+', 'UTILITAI8_36', sk=valk(1))
        end if
!
        if (ier .ge. 10) then
            valk(1) = nomsd
            valk(2) = nomch
            valk(3) = acces
            valr = rval
            call utmess(stp, 'UTILITAI8_37', nk=3, valk=valk, sr=valr)
        end if
    end if
!
    call jedema()
end subroutine
