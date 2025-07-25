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

subroutine numcima(infcha, nu, ccid, base)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
!
#include "asterfort/asschc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/numchc.h"
#include "asterfort/wkvect.h"
    character(len=19) :: infcha
    character(len=*) :: nu, ccid
    character(len=1) :: base
! ----------------------------------------------------------------------
!  BUT : ON NOTE LES DDLS ELIMINES PAR LES CHARGES CINEMATIQUES
!
!  REMARQUE : LE RESTE DU TRAITEMENT DES CHARGES CINEMATIQUES EST FAIT
!             LORS DE LA RESOLUTION (ASMCHC+CSMBGG)
!             ON VERIFIE LA COHERENCE AVEC ASSCHC
!
! IN  K*  INFCHA : / SD_INFCHA (K19)
!                  / NOM D'UN OBJET JEVEUX (K24) CONTENANT
!                    LES NOMS DES CHARGES (K8)
! IN  NU      : NUME_DDL
! IN  BASE   : 'G' / 'V'
!
!----------------------------------------------------------------------
!     VARIABLES LOCALES
!----------------------------------------------------------------------
    character(len=19) :: infch2
    integer(kind=8) :: iret, iret1, iret2, iret3, ich, ncharg, jlchci
    integer(kind=8) :: nchci
    integer(kind=8), pointer :: infc(:) => null()
    character(len=24), pointer :: lcha(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    infch2 = infcha
    if (infch2 .eq. ' ') goto 999
!
    call jeexin(infch2//'.LCHA', iret1)
    call jeexin(infch2//'.INFC', iret2)
    call jeexin(infcha, iret3)
    if (iret1+iret2+iret3 .eq. 0) goto 999
!
!
!     -- CAS SD_INFCHA :
    if (iret1+iret2 .gt. 0) then
        call jeexin(infch2//'.LCHA', iret)
        if (iret .eq. 0) goto 999
!
        call jeveuo(infch2//'.LCHA', 'L', vk24=lcha)
        call jeveuo(infch2//'.INFC', 'L', vi=infc)
!
        ncharg = infc(1)
        if (ncharg .eq. 0) goto 999
        call wkvect('&&NUCIMA.LCHCI', 'V V K24', ncharg, jlchci)
!
        nchci = 0
        do ich = 1, ncharg
!
!         -- CAS DES SD_CHAR_CINE :
            if (infc(ich+1) .lt. 0) then
                nchci = nchci+1
                zk24(jlchci-1+nchci) = lcha(ich)
            end if
        end do
        if (nchci .eq. 0) goto 999
        call numchc(nu, ccid, nchci, zk24(jlchci), base)
        go to 100
!
!
!     -- CAS LISTE DE CHARGES CINEMATIQUES :
    else
        call jeveuo(infcha, 'L', jlchci)
        call jelira(infcha, 'LONMAX', nchci)
        call numchc(nu, ccid, nchci, zk24(jlchci), base)
        go to 100
    end if
!
!
999 continue
    call numchc(nu, ccid, 0, ["XXX"], base)
100 continue
!
    call jedetr('&&NUCIMA.LCHCI')
    call jedema()
end subroutine
