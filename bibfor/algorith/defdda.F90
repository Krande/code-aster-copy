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
subroutine defdda(nbec, nbcmp, numgd, ioc, motcle, &
                  iopt, icod)
!    P. RICHARD     DATE 18/02/91
!-----------------------------------------------------------------------
!  BUT:  DETERMINER LA LISTE DES TYPE DDL DEFINIS PAR L'UTILISATEUR
    implicit none
!       EN ARGUMENT D'UN MOT-CLE
!     INDEPENDAMENT DES DDL ACTIFS DANS LE MODELE
!             IL SORT UN ENTIER CODE
!  TRAITEMENT DU CAS DE L'ABSENCE DE MOT-CLE PAR IOPT
!-----------------------------------------------------------------------
!
! NBEC     /I/: NOMBRE D'ENTIER CODES GRANDEUR SOUS-JACENTE
! NBCMP    /I/: NOMBRE DE COMPOSANTE MAX DE LA GRANDEUR SOUS-JACENTE
! NUMGD    /I/: NUMERO DE LA GRANDEUR SOUS-JACENTE
! IOC      /I/: NUMERO OCCURENCE MOTFAC INTERFACE DEFINISSANT LES DDL
! MOTCLE   /I/: MOT CLE
! IOPT     /I/: CODE POUR ABSENCE MOT-CLE (1 TOUT DDL) (0 AUCUN DDL)
! ICOD     /O/: ENTIER CODE
!
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/getvtx.h"
#include "asterfort/iscode.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomcou
    character(len=*) :: motcle
    character(len=24) :: temddl, temidc
    character(len=24) :: valk
    integer(kind=8) :: nbec, icod(nbec)
    aster_logical :: ok, okg
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, iec, ioc, iopt, j, llncmp
    integer(kind=8) :: ltddl, ltidec, nbcmp, nbval, numgd
!-----------------------------------------------------------------------
    data okg/.false./
!-----------------------------------------------------------------------
!
    call jemarq()
!
    if (motcle(1:9) .eq. 'DDL_ACTIF') then
        nbval = 0
    else
        call getvtx('INTERFACE', motcle, iocc=ioc, nbval=0, nbret=nbval)
        nbval = -nbval
    end if
!
!----------ALLOCATION DU VECTEUR DES ENTIERS DE DECODAGE----------------
!
    temidc = '&&DEFDDA.IDEC'
    call wkvect(temidc, 'V V I', nbcmp, ltidec)
!
!--------------TRAITEMENT DES EXCEPTIONS: PAS DE MOT CLE----------------
!
    if (nbval .eq. 0 .and. iopt .eq. 1) then
        do i = 1, nbcmp
            zi(ltidec+i-1) = 1
        end do
        call iscode(zi(ltidec), icod, nbcmp)
        goto 999
    end if
!
    if (nbval .eq. 0 .and. iopt .eq. 0) then
        do iec = 1, nbec
            icod(iec) = 0
        end do
        goto 999
    end if
!
!---------RECUPERATION DU VECTEUR DES NOMS DE COMPOSANTES---------------
!
    call jeveuo(jexnum('&CATA.GD.NOMCMP', numgd), 'L', llncmp)
!
    temddl = '&&DEFDDA.DDL.DON'
    call wkvect(temddl, 'V V K80', nbval, ltddl)
!
    if (motcle(1:9) .eq. 'DDL_ACTIF') then
        ibid = 0
    else
        call getvtx('INTERFACE', motcle, iocc=ioc, nbval=nbval, vect=zk80(ltddl), &
                    nbret=ibid)
    end if
!
    do i = 1, nbval
        nomcou = zk80(ltddl+i-1)
        ok = .true.
        do j = 1, nbcmp
            if (nomcou .eq. zk8(llncmp+j-1)) then
                zi(ltidec+j-1) = 1
                ok = .false.
            end if
        end do
!
        if (ok) then
            okg = .true.
            valk = nomcou
            call utmess('E+', 'ALGORITH15_8', sk=valk)
            call utmess('E', 'VIDE_1')
        end if
!
    end do
!
    if (okg) then
        call utmess('F', 'ALGORITH15_10')
    end if
!
    call iscode(zi(ltidec), icod, nbcmp)
!
    call jedetr(temddl)
!
999 continue
    call jedetr(temidc)
!
    call jedema()
end subroutine
