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
subroutine lirtet(ifl, ilec, inom, cnl, nom, &
                  icl, iv, rv, cv, deblig)
    implicit none
!       BUT: SAUT D'UNE EVENTUELLE ENTETE AVEC LECTURE DE "NOM="
!       IN: IFL : FICHIER OU ON LIT
!          ILEC: 1=> PREMIERE LECTURE, 2=> SECONDE LECTURE
!          INOM: 1=> LECTURE DU NOM,   0=> NOM IGNORE
!          CNL : NUMERO DE LA LIGNE DANS UNE CHAINE
!          DEBLIG : IN ET OUT DANS LIRITM
!       OUT:
!          NOM : ITEM DERRIERE "NOM="
!          ICL, IV, RV, CV, DEBLIG : COMME LIRITM
!       EN SORTIE DE CETTE ROUTINE, ON A LU LA PREMIERE LIGNE
!       D'INFORMATIONS
!       ----------------------------------------------------------------
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/liritm.h"
#include "asterfort/lirlig.h"
#include "asterfort/utmess.h"
    character(len=14) :: cnl
    character(len=80) :: lig
    character(len=8) :: cvz
    character(len=24) :: nom
    character(len=*) :: cv
    integer(kind=8) :: deblig
    real(kind=8) :: rv
    aster_logical :: lnom, lent
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icl, ifl, ilec, inom, iv
    integer(kind=8) :: lcv, nbigno
!-----------------------------------------------------------------------
    lnom = .false.
    lent = .false.
    nbigno = 0
    nom = 'INDEFINI'
!
    ASSERT(inom .eq. 0 .or. inom .eq. 1)
!
    if (inom .eq. 0) then
1       continue
        cv = ' '
        iv = 0
        rv = 0.d0
        call liritm(ifl, icl, iv, rv, cv(1:8), &
                    cnl, deblig, ilec)
        if (deblig .eq. 0) then
! -     IL Y A UNE ENTETE
            if (icl .eq. 3) then
                if (cv(1:6) .eq. 'NBLIGE') then
                    call liritm(ifl, icl, iv, rv, cv(1:8), &
                                cnl, deblig, ilec)
                    if (icl .eq. 1) then
                        nbigno = iv
                    else
                        ASSERT(.false.)
                    end if
                    goto 9
                else
! -         L'IDENTIFICATEUR LU N'EST PAS "NBLIGE"
                    goto 1
                end if
            else
! -       L'ITEM LU N'EST PAS UN IDENTIFICATEUR
                goto 1
            end if
        else
! -     PAS D'ENTETE
            goto 9
        end if
!
    else if (inom .eq. 1) then
2       continue
        cv = ' '
        iv = 0
        rv = 0.d0
        call liritm(ifl, icl, iv, rv, cv(1:8), &
                    cnl, deblig, ilec)
        if (deblig .eq. 0) then
! -     IL Y A UNE ENTETE
            if (icl .eq. 3) then
! -       L'ITEM LU EST UN IDENTIFICATEUR
                if (cv(1:3) .eq. 'NOM') then
                    call liritm(ifl, icl, iv, rv, cv(1:24), &
                                cnl, deblig, ilec)
                    if (icl .eq. 3) then
                        if (iv .gt. 24) then
                            cvz = cv
                            call utmess('A', 'MODELISA4_97', sk=cvz(1:iv))
                        end if
                        lcv = min(iv, 24)
                        ASSERT(len(cv) .ge. lcv)
                        nom = cv(1:lcv)
                    else
                        call utmess('F', 'MODELISA4_98')
                    end if
                    lnom = .true.
                else if (cv(1:6) .eq. 'NBLIGE') then
                    cv = ' '
                    iv = 0
                    rv = 0.d0
                    call liritm(ifl, icl, iv, rv, cv(1:8), &
                                cnl, deblig, ilec)
                    if (icl .eq. 1) then
                        nbigno = iv
                    else
                        ASSERT(.false.)
                    end if
                    lent = .true.
                else
! -          L'IDENTIFICATEUR LU N'EST NI "NBLIGE", NI "NOM"
                    goto 2
                end if
                if (lnom .and. lent) then
                    goto 9
                else
                    goto 2
                end if
            else
! -       L'ITEM LU N'EST PAS UN IDENTIFICATEUR
                goto 2
            end if
        else
            if (lent) then
                ASSERT(.not. lnom)
                nbigno = nbigno-1
                goto 9
            else
! -       PAS D'ENTETE
                goto 9
            end if
        end if
    end if
!
9   continue
    if (nbigno .gt. 0) then
        do i = 1, nbigno-1
            call lirlig(ifl, cnl, lig, ilec)
        end do
        deblig = -1
        call liritm(ifl, icl, iv, rv, cv(1:8), &
                    cnl, deblig, ilec)
    end if
!
end subroutine
