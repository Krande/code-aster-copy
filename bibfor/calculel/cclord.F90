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
subroutine cclord(nuoplo, nbordr, lisord, nobase, optdem, &
                  minord, maxord, resuin, resuou, lisout)
    implicit none
!     --- ARGUMENTS ---
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
!
    aster_logical :: optdem
    integer(kind=8) :: nuoplo, nbordr, minord, maxord
    character(len=8) :: resuin, resuou, nobase
    character(len=19) :: lisord
    character(len=24) :: lisout
!  CALC_CHAMP - DETERMINATION LISTE DES NUMEROS D'ORDRE
!  -    -                     -                   ---
! ----------------------------------------------------------------------
!
!  CREATION D'UNE LISTE DE NUMEROS D'ORDRE POUR L'OPTION EN ARGUMENT
!  EN PRENANT EN COMPTE SA LISTE DE DEPENDANCE
!
!  IN  :
!   NUOPLO  I    INDICE DE L'OPTION POUR LAQUELLE ON SOUHAITE OBTENIR
!                LA LISTE DE NUMEROS D'ORDRE
!   NBORDR  I    NOMBRE DE NUMEROS D'ORDRE
!   LISORD  K19  LISTE DE NUMEROS D'ORDRE
!   NOBASE  K8   BASE DU NOM A PARTIR DE LAQUELLE LE NOM DES OBJETS DE
!                CCLIOP SERONT CONSTRUITS
!   OPTDEM  BOOL EST-CE QUE L'OPTION DEMANDEE EST UN OPTION DEMANDEE
!                PAR L'UTILISATEUR
!   MINORD  I    NUMERO D'ORDRE MIN
!   MAXORD  I    NUMERO D'ORDRE MAX
!   RESUIN  K8   NOM DE LA STRUCTURE DE DONNEES RESULTAT IN
!   RESUOU  K8   NOM DE LA STRUCTURE DE DONNEES RESULTAT OUT
!
!  OUT :
!   LISOUT  K24  NOM JEVEUX DE LA LISTE DE NUMERO D'ORDRE
! ----------------------------------------------------------------------
! person_in_charge: nicolas.sellenet at edf.fr
!
    integer(kind=8) :: jlisop, jliori, jlidep, jordop, ierd, inddeb, indfin
    integer(kind=8) :: nopous, iordr, curmax, curmin, iter, decal, numord, jlnoin
    integer(kind=8) :: jordo2, jlisde, jordr
!
    character(len=1) :: isodep
    character(len=5) :: numopt
    character(len=11) :: nobaop
    character(len=16) :: option
    character(len=19) :: nosyou
    character(len=24) :: noliop, nolori, noldep, noliin, nolisd
!
    call jemarq()
!
    call codent(nuoplo, 'D0', numopt)
    nobaop = nobase//'.OP'
    lisout = nobaop//numopt
    call jeveuo(lisord, 'L', jordr)
!
    isodep = ' '
    noliop = nobase//'.LISOPT'
    nolori = nobase//'.LISORI'
    noldep = nobase//'.LISDEP'
    noliin = nobase//'.LNOINS'
    nolisd = nobase//'.ISODEP'
!
    call jeveuo(noliop, 'L', jlisop)
    call jeveuo(nolori, 'L', jliori)
    call jeveuo(noldep, 'L', jlidep)
    call jeveuo(noliin, 'E', jlnoin)
    call jeveuo(nolisd, 'L', jlisde)
!
    option = zk24(jlisop+nuoplo-1)
    inddeb = zi(jliori+2*nuoplo-2)
    indfin = zi(jliori+2*nuoplo-1)
    isodep = zk8(jlisde+nuoplo-1)
!
    call wkvect(lisout, 'V V I ', nbordr+3, jordop)
!
    if (inddeb .ne. 0) then
!       CAS 1 : CETTE OPTION DEPEND D'AUTRES OPTIONS A CALCULER
!               AUQUEL CAS, IL FAUT REGARDER COMMENT ELLE EN DEPEND
!               ET LA LISTE DES NUMEROS D'ORDRE DE SES PARENTS
        call jeveuo(noliin, 'E', jlnoin)
        curmax = maxord
        curmin = minord
        do iter = inddeb, indfin
            call jeveuo(zk24(jlnoin+iter-1), 'L', jordo2)
!
            if (zk8(jlidep+iter-1) .eq. 'NP1') then
                decal = -1
            else if (zk8(jlidep+iter-1) .eq. 'NM1') then
                decal = +1
            else
                decal = 0
            end if
!
!         LA LISTE DE NUMEROS D'ORDRE PROVIENT DE OP0058
!         ELLE EST DONC CROISSANTE
            curmax = min(curmax, zi(jordo2+2)+decal)
            curmin = max(curmin, zi(jordo2+1)+decal)
        end do
!
        nopous = 0
        do iordr = 1, nbordr
            numord = zi(jordr-1+iordr)
            if ((isodep .eq. '-') .and. (numord .eq. minord)) then
                goto 30
            else if ((isodep .eq. '+') .and. (numord .eq. maxord)) then
                goto 30
            end if
            if (numord .ge. curmin) then
                if (numord .gt. curmax) goto 40
                nosyou = ' '
                ierd = 1
                if (.not. optdem) then
                    call rsexch(' ', resuin, option, numord, nosyou, &
                                ierd)
                end if
                if (ierd .ne. 0) then
                    call rsexch(' ', resuou, option, numord, nosyou, &
                                ierd)
                end if
!
                if (ierd .ne. 0) then
                    nopous = nopous+1
                    zi(jordop+nopous+2) = numord
                end if
            end if
30          continue
        end do
!
40      continue
        zi(jordop+1) = curmin
        zi(jordop+2) = curmax
        zi(jordop) = nopous
    else
!       CAS 2 : AUCUNE DEPENDANCE, ON PEUT DONC RECOPIER LA LISTE
!               DES NUMEROS D'ORDRE EN VERIFIANT QUE L'OPTION
!               N'EXISTE NI DANS RESUIN NI DANS RESUOU
        nopous = 0
        do iordr = 1, nbordr
            numord = zi(jordr-1+iordr)
            if ((isodep .eq. '-') .and. (numord .eq. minord)) then
                goto 20
            else if ((isodep .eq. '+') .and. (numord .eq. maxord)) then
                goto 20
            end if
            nosyou = ' '
            ierd = 1
            if (.not. optdem) then
                call rsexch(' ', resuin, option, numord, nosyou, &
                            ierd)
            end if
            if (ierd .ne. 0) then
                call rsexch(' ', resuou, option, numord, nosyou, &
                            ierd)
            end if
!
            if (ierd .ne. 0) then
                nopous = nopous+1
                zi(jordop+nopous+2) = numord
            end if
20          continue
        end do
        zi(jordop+1) = zi(jordr-1+1)
        zi(jordop+2) = zi(jordr-1+nbordr)
        zi(jordop) = nopous
    end if
!
    zk24(jlnoin+nuoplo-1) = lisout
!
    call jedema()
!
end subroutine
