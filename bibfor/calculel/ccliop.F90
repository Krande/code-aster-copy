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
subroutine ccliop(type, option, nobase, noliop, nopout)
    implicit none
!     --- ARGUMENTS ---
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nopout
    character(len=*) :: type
    character(len=8) :: nobase
    character(len=16) :: option
    character(len=24) :: noliop
!  CALC_CHAMP - DETERMINATION LISTE D'OPTIONS
!  -    -                     --      --
! ----------------------------------------------------------------------
!
!  CREATION D'UNE LISTE D'OPTIONS DONT OPTION DEPEND
!   EXEMPLE : LISTE DE TAILLE N
!       NOLIOP(N)   = OPTION (EN ARGUEMENT DE LA ROUTINE)
!       NOLIOP(N-1) = OPTIO2 (DONT DEPEND OPTION)
!       ...
!       NOLIOP(M)   = OPTIOM (DONT DEPEND OPTION)
!       ...
!       NOLIOP(P)   = OPTIOP (DONT DEPEND OPTIO2)
!       ...
!
!  D'OU NOLORI DE TAILLE 2*N
!       NOLORI(N)   = N-1
!       NOLORI(N-1) = M
!       NOLORI(3) = P
!       ...
!  LA LISTE NOLDEP RAPPELLE LA DEPENDANCE TEMPORELLE
!  EX : SI OPTION DEPEND DE L'INSTANT N+1 D'OPTIO2 ALORS :
!       NOLDEP(N-1) = 'NP1'
!
! IN  :
!   TYPE    K8   TYPE DE DEPENDANCE : OPTION OU CHAMP
!                - OPTION : OPTIONS A CALCULER
!                - CHAMP  : CHAMPS NECESSAIRES (POUR STANLEY)
!   OPTION  K24  NOM DE L'OPTION
!   NOBASE  K8   BASE DU NOM A PARTIR DE LAQUELLE LE NOM DES OBJETS DE
!                CCLIOP SERONT CONSTRUITS
! OUT :
!   NOLIOP  K24  NOM JEVEUX DE LA LISTE DES OPTIONS CREES
!   NOPOUT  I    TAILLE DE LA LISTE OUT
! ----------------------------------------------------------------------
! person_in_charge: nicolas.sellenet at edf.fr
    integer(kind=8) :: opt, iaopds, iaoplo, iopdeb, iop, nparin, ipara, opt2
    integer(kind=8) :: lopor1(100), lopor2(100), itmp, nopous, jlisop
    integer(kind=8) :: jliori, jlidep, jlnoin, jlisde
!
    aster_logical :: opajou
!
    character(len=1) :: isodep(100)
    character(len=4) :: lopdep(100)
    character(len=16) :: loptio(100), optio2, curopt
    character(len=24) :: nolori, noldep, noliin, nolisd
!
    call jemarq()
!
    nopout = 0
!
    noliop = nobase//'.LISOPT'
    nolori = nobase//'.LISORI'
    noldep = nobase//'.LISDEP'
    noliin = nobase//'.LNOINS'
    nolisd = nobase//'.ISODEP'
!
    if (option(6:9) .ne. 'NOEU') then
        call jenonu(jexnom('&CATA.OP.NOMOPT', option), opt)
        if (opt .eq. 0) then
            call utmess('F', 'CALCCHAMP_9', sk=option)
        end if
        call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iaopds)
        if (zi(iaopds-1+2) .eq. 0) goto 999
        call jeveuo(jexnum('&CATA.OP.LOCALIS', opt), 'L', iaoplo)
    end if
!
!     INITIALISATION DES ENTIERS
    iopdeb = 1
    nopout = 1
    nopous = nopout
!
!     REMPLISSAGE DU TABLEAU DES OPTIONS A CALCULER AVEC LA PREMIERE
!     OPTION DEMANDEE
    loptio(nopout) = option
    lopdep(nopout) = 'NSP'
!
40  continue
!
!     BOUCLE SUR LE TABLEAU DES OPTIONS QUI SERA ENRICHI A CHAQUE
!     PASSE
    do iop = iopdeb, nopous
        isodep(iop) = ' '
        curopt = loptio(iop)
        opajou = .false.
!
!       CAS D'UNE OPTION AUX NOEUDS
        if (curopt(6:9) .eq. 'NOEU') then
            optio2 = option(1:5)//'ELNO'
            call jenonu(jexnom('&CATA.OP.NOMOPT', optio2), opt)
            if (opt .ne. 0) then
                nopout = nopout+1
                loptio(nopout) = optio2
                if (.not. opajou) then
                    lopor1(iop) = nopout
                    lopor2(iop) = 0
                end if
                lopor1(nopout) = 0
                lopor2(nopout) = 0
                lopdep(nopout) = 'N'
                opajou = .true.
            end if
!
!       CAS D'UNE OPTION AUX ELEMENTS
        else
            if ((curopt .eq. 'DEPL') .or. (curopt .eq. 'VITE') .or. (curopt .eq. 'ACCE') .or. &
                (curopt .eq. 'TEMP') .or. (curopt .eq. 'VARI_ELGA')) then
!         POUR CES OPTIONS, IL N'Y A PAS DE DEPENDANCE
                nparin = 0
            else
                call jenonu(jexnom('&CATA.OP.NOMOPT', curopt), opt)
                call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iaopds)
                if (zi(iaopds-1+2) .eq. 0) goto 10
                call jeveuo(jexnum('&CATA.OP.LOCALIS', opt), 'L', iaoplo)
                nparin = zi(iaopds-1+2)
            end if
!
!         BOUCLE SUR LES PARAMETRES DE CETTE OPTION
            do ipara = 1, nparin
!           ON VERIFIE QUE L'OPTION CORRESPONDAND AU CHAMP EXISTE
                optio2 = zk24(iaoplo+3*ipara-2)
                call jenonu(jexnom('&CATA.OP.NOMOPT', optio2), opt2)
!
!           ON EVITE QU'UNE OPTION NE DEPENDE D'ELLE MEME
                if (loptio(iop) .eq. optio2) then
                    if (zk24(iaoplo+3*ipara-1) .eq. 'NP1') then
                        isodep(iop) = '+'
                    else if (zk24(iaoplo+3*ipara-1) .eq. 'NM1') then
                        isodep(iop) = '-'
!              ELSE
!                ASSERT(.FALSE.)
                    end if
                    goto 20
                end if
!
                if ((opt2 .ne. 0)) then
! --------- ON DEPEND D'UN CHAMP QUI EST UNE OPTION
                    nopout = nopout+1
                    loptio(nopout) = optio2
!             REMPLISSAGE DE LA LISTE DE DEPENDANCES DES OPTIONS
                    if (.not. opajou) then
                        lopor1(iop) = nopout
                        lopor2(iop) = 0
                    end if
                    lopor1(nopout) = 0
                    lopor2(nopout) = 0
                    lopdep(nopout) = zk24(iaoplo+3*ipara-1)
                    opajou = .true.
                    isodep(nopout) = ' '
                elseif (((optio2 .eq. 'DEPL') .or. (optio2 .eq. 'VITE') &
                         .or. (optio2 .eq. 'ACCE') .or. (optio2 .eq. 'TEMP') .or. &
                         (optio2 .eq. 'VARI_ELGA')) .and. (zk24(iaoplo+3*ipara-1) .eq. 'N') &
                        .and. (type(1:5) .eq. 'CHAMP')) then
! --------- ON DEPEND D'UN CHAMP QUI N'EST PAS UNE OPTION (POUR STANLEY)
                    nopout = nopout+1
                    loptio(nopout) = optio2
!             REMPLISSAGE DE LA LISTE DE DEPENDANCES DES OPTIONS
                    if (.not. opajou) then
                        lopor1(iop) = nopout
                        lopor2(iop) = 0
                    end if
                    lopor1(nopout) = 0
                    lopor2(nopout) = 0
                    lopdep(nopout) = zk24(iaoplo+3*ipara-1)
                    opajou = .true.
                    isodep(nopout) = ' '
                elseif (((optio2 .eq. 'DEPL') .or. (optio2 .eq. 'SIEF_ELGA') &
                         .or. (optio2 .eq. 'VARI_ELGA')) .and. &
                        ((zk24(iaoplo+3*ipara-1) .eq. 'NP1') &
                         .or. (zk24(iaoplo+3*ipara-1) (1:3) .eq. 'NM1'))) then
! --------- ON DEPEND DE SOIT-MEME A L'INSTANT NP1 OU NM1
                    nopout = nopout+1
                    loptio(nopout) = optio2
!             REMPLISSAGE DE LA LISTE DE DEPENDANCES DES OPTIONS
                    if (.not. opajou) then
                        lopor1(iop) = nopout
                        lopor2(iop) = 0
                    end if
                    lopor1(nopout) = 0
                    lopor2(nopout) = 0
                    lopdep(nopout) = zk24(iaoplo+3*ipara-1)
                    opajou = .true.
                    isodep(nopout) = ' '
                end if
20              continue
            end do
        end if
!
        if (.not. opajou) then
            lopor1(iop) = 0
            lopor2(iop) = 0
        else
            lopor2(iop) = nopout
        end if
10      continue
    end do
!
!     SI ON A AJOUTE UNE OPTION LORS DE LA DERNIERE PASSE, ON
!     DOIT CHERCHER SES DEPENDANCES
    if (opajou) then
        iopdeb = nopous+1
        nopous = nopout
        goto 40
    end if
!
!     TEMPORAIRE POUR EVITER LES DEPASSEMENTS DE TABLEAU
    if (nopout .gt. 100) then
        ASSERT(.false.)
    end if
!
    call wkvect(noliop, 'V V K24', nopout, jlisop)
    call wkvect(nolori, 'V V I', 2*nopout, jliori)
    call wkvect(noldep, 'V V K8', nopout, jlidep)
    call wkvect(noliin, 'V V K24', nopout, jlnoin)
    call wkvect(nolisd, 'V V K8', nopout, jlisde)
!
!     CONSTRUCTION DE LA LISTE DES PAS DE TEMPS
    do iop = 1, nopout
!       ON PARCOURT LA LISTE A L'ENVERS PUISQUE PAR CONSTRUCTION
!       LES OPTIONS 'D'EN HAUT' DEPENDENT DES OPTIONS 'D'EN BAS'
        itmp = nopout-iop+1
        optio2 = loptio(itmp)
        zk24(jlisop+iop-1) = optio2
!
!       ON REGARDE DE QUELLES OPTIONS DEPEND L'OPTION COURANTE
        if (lopor1(itmp) .ne. 0) then
            zi(jliori+2*iop-2) = nopout-lopor2(itmp)+1
            zi(jliori+2*iop-1) = nopout-lopor1(itmp)+1
        else
            zi(jliori+2*iop-2) = 0
            zi(jliori+2*iop-1) = 0
        end if
        zk8(jlidep+iop-1) = lopdep(itmp)
        zk8(jlisde+iop-1) = isodep(itmp)
    end do
!
999 continue
!
    call jedema()
!
end subroutine
