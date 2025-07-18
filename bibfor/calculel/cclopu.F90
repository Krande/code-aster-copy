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
subroutine cclopu(resuin, resuou, lisord, nbordr, lisopt, &
                  nbropt)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbordr, nbropt
    character(len=8) :: resuin, resuou
    character(len=19) :: lisord, lisopt
! person_in_charge: nicolas.sellenet at edf.fr
! ----------------------------------------------------------------------
!  CALC_CHAMP - DETERMINATION DE LA LISTE D'OPTIONS DE L'UTILISATEUR
!  -    -                           -       --           -
! ----------------------------------------------------------------------
!
!  ROUTINE SERVANT A RECUPERER LA LISTE D'OPTIONS DE CALC_CHAMP
!    SOUHAITEES PAR L'UTILISATEUR
!
!  LES OPTIONS SONT FOURNIES DANS DES MOTS-CLES SIMPLES :
!    - CONTRAINTE
!    - DEFORMATION
!    - ENERGIE
!    - CRITERES
!    - VARI_INTERNE
!    - HYDRAULIQUE
!    - THERMIQUE
!    - ACOUSTIQUE
!    - FORCE
!    - PROPRIETES
!    - SOUS_POINT
!  ET LE MOT-CLE FACTEUR CHAM_UTIL.
!
! IN  :
!   RESUIN K8   NOM DE LA SD IN
!
! IN/OUT :
!   LISOPT K19  NOM DE LA LISTE D'OPTIONS
!
! OUT :
!   NBOPT  I    NOMBRE D'OPTIONS
! ----------------------------------------------------------------------
    integer(kind=8) :: ntymax
    parameter(ntymax=11)
!
    integer(kind=8) :: i, ityp, n1, jopt, postmp, nbopfa, ioc, ibid
    integer(kind=8) :: nuti, nsup, jord, iordr, iret
!
    character(len=9) :: mcfact
    character(len=12) :: typopt, tygrop(ntymax)
    character(len=16) :: option
    character(len=24) :: chn
    parameter(mcfact='CHAM_UTIL')
!
    aster_logical :: newcal, vu
    integer(kind=8), pointer :: nb_op_ty(:) => null()
    character(len=16), pointer :: oputil(:) => null()
!
    data tygrop/'CONTRAINTE  ', 'DEFORMATION ', 'ENERGIE     ', &
        'CRITERES    ', 'VARI_INTERNE', 'HYDRAULIQUE ', &
        'THERMIQUE   ', 'ACOUSTIQUE  ', 'FORCE       ', &
        'PROPRIETES  ', 'SOUS_POINT  '/
!
    call jemarq()
!
! --- PREMIERE BOUCLE POUR DETERMINER LE NOMBRE TOTAL D'OPTIONS
    AS_ALLOCATE(vi=nb_op_ty, size=ntymax)
    nbropt = 0
    do ityp = 1, ntymax
        typopt = tygrop(ityp)
        call getvtx(' ', typopt, nbval=0, nbret=n1)
        nb_op_ty(ityp) = -n1
        nbropt = nbropt-n1
    end do
!
    call wkvect(lisopt, 'V V K16', max(1, nbropt), jopt)
!
!     DEUXIEME BOUCLE POUR REMPLIR LE TABLEAU DES OPTIONS
    postmp = 0
    do ityp = 1, ntymax
        typopt = tygrop(ityp)
        nbopfa = nb_op_ty(ityp)
        if (nbopfa .eq. 0) goto 20
        call getvtx(' ', typopt, nbval=nbopfa, vect=zk16(jopt+postmp), nbret=n1)
        postmp = postmp+nbopfa
20      continue
    end do
!
! --- MOT-CLE FACTEUR CHAM_UTIL
!     POUR EVITER L'ALARME LIE AU RECALCUL D'UNE OPTION DEJA PRESENTE
!     ON REGARDE SI ELLE A DEJA ETE CALCULEE ET SI ELLE N'EST PAS DEJA
!     DANS LA LISTE DES OPTIONS DEMANDEES PAR L'UTILISATEUR
    call getfac(mcfact, nuti)
    if (nuti .eq. 0) then
        goto 999
    end if
!
    newcal = .false.
    call jeexin(resuou//'           .DESC', iret)
    if (iret .eq. 0) newcal = .true.
!
    AS_ALLOCATE(vk16=oputil, size=nuti)
    call jeveuo(lisord, 'L', jord)
    nsup = 0
    do ioc = 1, nuti
        call getvtx(mcfact, 'NOM_CHAM', iocc=ioc, scal=option, nbret=ibid)
        vu = .true.
!       OPTION PRESENTE DANS RESUIN A TOUS LES NUME_ORDRE A CALCULER ?
        do i = 1, nbordr
            iordr = zi(jord-1+i)
            call rsexch(' ', resuin, option, iordr, chn, &
                        iret)
            if (iret .ne. 0) then
                if (.not. newcal) call rsexch(' ', resuou, option, iordr, chn, &
                                              iret)
                if (iret .ne. 0) then
                    vu = .false.
                    goto 32
                end if
            end if
        end do
        goto 38
!
32      continue
!       OPTION DEJA DANS LA LISTE ?
        vu = .false.
        do i = 1, nbropt
            if (zk16(jopt-1+i) .eq. option) then
                vu = .true.
                goto 30
            end if
        end do
        do i = 1, nsup
            if (oputil(i) .eq. option) then
                vu = .true.
                goto 30
            end if
        end do
!
38      continue
!       ON AJOUTE L'OPTION A LA LISTE
        if (.not. vu) then
            nsup = nsup+1
            oputil(nsup) = option
        end if
30      continue
    end do
!
! --- REFAIRE OU AGRANDIR LISOPT
    if (nsup .gt. 0) then
        if (nbropt .eq. 0) then
            call jedetr(lisopt)
            call wkvect(lisopt, 'V V K16', nsup, jopt)
        else
            call juveca(lisopt, nbropt+nsup)
        end if
        do i = 1, nsup
            zk16(jopt-1+nbropt+i) = oputil(i)
        end do
        nbropt = nbropt+nsup
    end if
!
999 continue
    AS_DEALLOCATE(vi=nb_op_ty)
    AS_DEALLOCATE(vk16=oputil)
    call jedema()
!
end subroutine
