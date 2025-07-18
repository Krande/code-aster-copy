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

subroutine rsmode(resultZ)
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
    character(len=*) :: resultZ
!
!     FONCTION  :
!     SI LES CHAM_NO (DEPL_R) DE LA SD RESU NE SONT PAS NUMEROTES
!     COMME LE NUME_DDL DU .REFD, ON LES RENUMEROTE.
!
!    IN/JXVAR : RESU  NOM DU CONCEPT SD_RESULTAT
!-----------------------------------------------------------------------
    integer(kind=8) :: iret, neq, numeStore, isymb, k, iStore
    integer(kind=8) :: nbnosy, nbStore, iexi, nbval, jliprf, nbval2
    character(len=1) :: kbid, typ1
    character(len=8) :: result8, nomgd, ma1, fieldMesh
    character(len=19) :: result19
    character(len=14) :: numeDof
    character(len=16) :: fieldType
    character(len=19) :: fieldNumeEqua, champt, fieldName, numeEqua
    integer(kind=8), pointer :: listStore(:) => null()
!-----------------------------------------------------------------------
!
!
    call jemarq()
    result19 = resultZ
    result8 = resultZ
    call jeexin(result19//'.REFD', iexi)
    if (iexi .eq. 0) goto 50
!
    numeDof = ' '
    call dismoi('NUME_DDL', resultZ, 'RESU_DYNA', repk=numeDof, arret='C', ier=iret)
!
    if (numeDof .eq. ' ') goto 50
!
    numeEqua = numeDof(1:8)//'.NUME'
    call dismoi('NOM_MAILLA', numeDof, 'NUME_DDL', repk=ma1)
!
    champt = '&&RSMODE.CHAMPT'
!
    call jeveuo(result19//'.ORDR', 'L', vi=listStore)
    call jelira(result19//'.ORDR', 'LONUTI', nbStore)
    call jelira(result19//'.DESC', 'NOMUTI', nbnosy)
!
    do isymb = 1, nbnosy
        call jenuno(jexnum(result19//'.DESC', isymb), fieldType)
!
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            call rsexch(' ', resultZ, fieldType, numeStore, fieldName, iret)
            if (iret .ne. 0) cycle
!
            call dismoi('NOM_GD', fieldName, 'CHAM_NO', repk=nomgd)
            if (nomgd(1:5) .ne. 'DEPL_') cycle
!
            call dismoi('NUME_EQUA', fieldName, 'CHAM_NO', repk=fieldNumeEqua)
            if (fieldNumeEqua .eq. numeEqua) cycle
!
            call dismoi('NOM_MAILLA', fieldName, 'CHAM_NO', repk=fieldMesh)
            if (ma1 .ne. fieldMesh) then
                call utmess('F', 'UTILITAI_29')
            end if
!
!        -- SI LE CHAMP NOMCHA N'A PAS LA BONNE NUMEROTATION,
!           IL FAUT LA MODIFIER :
            call jelira(fieldName//'.VALE', 'TYPE', cval=typ1)
            call vtcreb(champt, 'V', typ1, nume_ddlz=numeDof, nb_equa_outz=neq)
            call vtcopy(fieldName, champt, iret)
            if (iret .ne. 0) then
                call utmess("F", "FIELD0_14", sk=fieldType)
            end if
            call detrsd('CHAM_NO', fieldName)
            call copisd('CHAMP', 'G', champt, fieldName)
            call detrsd('CHAM_NO', champt)

        end do
    end do

!     -- IL FAUT ENCORE DETRUIRE LES NUME_EQUA QUI ONT ETE CREES
!        INUTILEMENT SUE LA BASE GLOBALE (POUR SDVERI=OUI)
    if (resultZ .ne. numeDof(1:8)) then
        call jelstc('G', result8//'.PRFCN', 1, 0, kbid, nbval)
        if (nbval .lt. 0) then
            nbval = -nbval
            call wkvect('&&RSMODES.LIPRFCN', 'V V K24', nbval, jliprf)
            call jelstc('G', result8//'.PRFCN', 1, nbval, zk24(jliprf), &
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
