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
subroutine mavegr(nomu)
    implicit none
#include "jeveux.h"
#include "asterfort/cpclma.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: nomu
!
!     SUPPRESSION DES GROUPES DE NOEUDS OU MAILLES DE NOM '      '
! ----------------------------------------------------------------------
!
    integer(kind=8) :: iret, i, j, nbgrma, nbgrmt, nbgrno, nbgrnt, nbma, nbno, jvg, jgg
    character(len=24) :: grpnoe, grpnov, grpmai, grpmav, gpptnn, gpptnm
    character(len=24) :: nomg, blanc
! ----------------------------------------------------------------------
    call jemarq()
!
    grpnoe = nomu//'.GROUPENO'
    grpnov = '&&MAVEGR.GROUPENO'
    gpptnn = nomu//'.PTRNOMNOE'
    grpmai = nomu//'.GROUPEMA'
    grpmav = '&&MAVEGR.GROUPEMA'
    gpptnm = nomu//'.PTRNOMMAI'
    blanc = ' '
!
! --- TRAITEMENT DES GROUP_MA
!
    call jeexin(grpmai, iret)
    if (iret .gt. 0) then
        call jelira(grpmai, 'NMAXOC', nbgrma)
        nbgrmt = nbgrma
        do i = 1, nbgrma
            call jeexin(jexnum(grpmai, i), iret)
            if (iret .eq. 0) goto 100
            call jenuno(jexnum(grpmai, i), nomg)
            if (nomg .eq. blanc) then
                nbgrmt = nbgrmt-1
                call utmess('A', 'MODELISA5_36')
            end if
100         continue
        end do
        if (nbgrmt .eq. 0) then
            call jedetr(grpmai)
        else if (nbgrmt .ne. nbgrma) then
            call cpclma(nomu, '&&MAVEGR', 'GROUPEMA', 'V')
            call jedetr(grpmai)
            call jedetr(gpptnm)
            call jecreo(gpptnm, 'G N K24')
            call jeecra(gpptnm, 'NOMMAX', nbgrmt)
            call jecrec(grpmai, 'G V I', 'NO '//gpptnm, 'DISPERSE', 'VARIABLE', &
                        nbgrmt)
            do i = 1, nbgrma
                call jeexin(jexnum(grpmav, i), iret)
                if (iret .eq. 0) goto 110
                call jenuno(jexnum(grpmav, i), nomg)
                if (nomg .eq. blanc) goto 110
                call jecroc(jexnom(grpmai, nomg))
                call jeveuo(jexnum(grpmav, i), 'L', jvg)
                call jelira(jexnum(grpmav, i), 'LONUTI', nbma)
                call jeecra(jexnom(grpmai, nomg), 'LONMAX', max(1, nbma))
                call jeecra(jexnom(grpmai, nomg), 'LONUTI', nbma)
                call jeveuo(jexnom(grpmai, nomg), 'E', jgg)
                do j = 0, nbma-1
                    zi(jgg+j) = zi(jvg+j)
                end do
110             continue
            end do
            call jedetr(grpmav)
        end if
    end if
!
! --- TRAITEMENT DES GROUP_NO
!
    call jeexin(grpnoe, iret)
    if (iret .gt. 0) then
        call jelira(grpnoe, 'NMAXOC', nbgrno)
        nbgrnt = nbgrno
        do i = 1, nbgrno
            call jeexin(jexnum(grpnoe, i), iret)
            if (iret .eq. 0) goto 200
            call jenuno(jexnum(grpnoe, i), nomg)
            if (nomg .eq. blanc) then
                nbgrnt = nbgrnt-1
                call utmess('A', 'MODELISA5_37')
            end if
200         continue
        end do
        if (nbgrnt .eq. 0) then
            call jedetr(grpnoe)
        else if (nbgrnt .ne. nbgrno) then
            call cpclma(nomu, '&&MAVEGR', 'GROUPENO', 'V')
            call jedetr(grpnoe)
            call jedetr(gpptnn)
            call jecreo(gpptnn, 'G N K24')
            call jeecra(gpptnn, 'NOMMAX', nbgrnt)
            call jecrec(grpnoe, 'G V I', 'NO '//gpptnn, 'DISPERSE', 'VARIABLE', &
                        nbgrnt)
            do i = 1, nbgrno
                call jeexin(jexnum(grpnov, i), iret)
                if (iret .eq. 0) goto 210
                call jenuno(jexnum(grpnov, i), nomg)
                if (nomg .eq. blanc) goto 210
                call jecroc(jexnom(grpnoe, nomg))
                call jeveuo(jexnum(grpnov, i), 'L', jvg)
                call jelira(jexnum(grpnov, i), 'LONUTI', nbno)
                call jeecra(jexnom(grpnoe, nomg), 'LONMAX', max(1, nbno))
                call jeecra(jexnom(grpnoe, nomg), 'LONUTI', nbno)
                call jeveuo(jexnom(grpnoe, nomg), 'E', jgg)
                do j = 0, nbno-1
                    zi(jgg+j) = zi(jvg+j)
                end do
210             continue
            end do
            call jedetr(grpnov)
        end if
    end if
!
    call jedema()
!
end subroutine
