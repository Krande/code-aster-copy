! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine rsmode(resu)
    implicit none
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jelstc.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcopy.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: resu
! person_in_charge: jacques.pellet at edf.fr
!
!     FONCTION  :
!     SI LES CHAM_NO (DEPL_R) DE LA SD RESU NE SONT PAS NUMEROTES
!     COMME LE NUME_DDL DU .REFD, ON LES RENUMEROTE.
!
!    IN/JXVAR : RESU  NOM DU CONCEPT SD_RESULTAT
!-----------------------------------------------------------------------
    integer :: iret, neq, iordr, isymb, k, krang
    integer :: nbnosy, nbordr, iexi, nbval, jliprf, nbval2
    character(len=1) :: kbid, typ1
    character(len=8) :: resu8, nomgd, ma1, ma2
    character(len=19) :: resu19
    character(len=14) :: nu
    character(len=16) :: nomsym
    character(len=19) :: nume_equao, champt, nomcha, nume_equa1
    character(len=24) :: valk(5)
    integer, pointer :: ordr(:) => null()
!-----------------------------------------------------------------------
!
!
    call jemarq()
    resu19 = resu
    resu8 = resu
    call jeexin(resu19//'.REFD', iexi)
    if (iexi .eq. 0) goto 50
!
    nu = ' '
    call dismoi('NUME_DDL', resu, 'RESU_DYNA', repk=nu, arret='C', &
                ier=iret)
!
    if (nu .eq. ' ') goto 50
!
    nume_equa1 = nu(1:8)//'.NUME'
    call dismoi('NOM_MAILLA', nu, 'NUME_DDL', repk=ma1)
!
    champt = '&&RSMODE.CHAMPT'
!
    call jeveuo(resu19//'.ORDR', 'L', vi=ordr)
    call jelira(resu19//'.ORDR', 'LONUTI', nbordr)
    call jelira(resu19//'.DESC', 'NOMUTI', nbnosy)
!
    do isymb = 1, nbnosy
        call jenuno(jexnum(resu19//'.DESC', isymb), nomsym)
!
        do krang = 1, nbordr
            iordr = ordr(krang)
            call rsexch(' ', resu, nomsym, iordr, nomcha, &
                        iret)
            if (iret .ne. 0) goto 30
!
            call dismoi('NOM_GD', nomcha, 'CHAM_NO', repk=nomgd)
            if (nomgd(1:5) .ne. 'DEPL_') goto 30
!
            call dismoi('NUME_EQUA', nomcha, 'CHAM_NO', repk=nume_equao)
            if (nume_equao .eq. nume_equa1) goto 30
!
            call dismoi('NOM_MAILLA', nomcha, 'CHAM_NO', repk=ma2)
            if (ma1 .ne. ma2) then
                valk(1) = resu
                valk(2) = ma1
                valk(3) = nu
                valk(4) = ma2
                call utmess('F', 'UTILITAI_29', nk=4, valk=valk)
            end if
!
!        -- SI LE CHAMP NOMCHA N'A PAS LA BONNE NUMEROTATION,
!           IL FAUT LA MODIFIER :
            call jelira(nomcha//'.VALE', 'TYPE', cval=typ1)
            call vtcreb(champt, 'V', typ1, nume_ddlz=nu, nb_equa_outz=neq)
            call vtcopy(nomcha, champt, 'F', iret)
            call detrsd('CHAM_NO', nomcha)
            call copisd('CHAMP', 'G', champt, nomcha)
            call detrsd('CHAM_NO', champt)
30          continue
        end do
    end do
!
!
!     -- IL FAUT ENCORE DETRUIRE LES NUME_EQUA QUI ONT ETE CREES
!        INUTILEMENT SUE LA BASE GLOBALE (POUR SDVERI=OUI)
    if (resu .ne. nu(1:8)) then
        call jelstc('G', resu8//'.PRFCN', 1, 0, kbid, &
                    nbval)
        if (nbval .lt. 0) then
            nbval = -nbval
            call wkvect('&&RSMODES.LIPRFCN', 'V V K24', nbval, jliprf)
            call jelstc('G', resu8//'.PRFCN', 1, nbval, zk24(jliprf), &
                        nbval2)
            do k = 1, nbval2
                call detrsd('NUME_EQUA', zk24(jliprf-1+k) (1:19))
            end do
        end if
    end if
!
50  continue
    call jedema()
end subroutine
