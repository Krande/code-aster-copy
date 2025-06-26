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
function idenob(obj1, obj2)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
!
    aster_logical :: idenob
    character(len=*) :: obj1, obj2
!
! but : determiner l'identite de 2 objets jeveux
!  in  k24    obj1   : nom du 1er objet jeveux
!  in  k24    obj2   : nom du 2eme objet jeveux
!
! resultat:
!   idenob : .true.    si obj1(*) == obj2(*)
!                      (ou si obj1 et obj2 sont inexistants tous les 2)
!            .false.   sinon
!----------------------------------------------------------------------
    character(len=24) :: ob1, ob2, k241, k242
    character(len=1) :: typ1, typ2
    character(len=4) :: genr1, genr2, genr, xous1, xous2, xous, type
    integer(kind=8) :: iret1, iret2, ltyp1, ltyp2, ltyp, l, l1, l2, k, iad1, iad2
    integer(kind=8) :: iobj, nbobj
!----------------------------------------------------------------------
    call jemarq()
    ob1 = obj1
    ob2 = obj2
    idenob = .true.
!
    if (ob1 .eq. ob2) goto 220
!
    call jeexin(ob1, iret1)
    call jeexin(ob2, iret2)
    if (iret1 .eq. 0) then
        if (iret2 .gt. 0) goto 210
        goto 220
    else
        if (iret2 .eq. 0) goto 210
    end if
!
!
!   -- determination de type: r/c/k8,...
!   -------------------------------------
    call jelira(ob1, 'TYPE', cval=typ1)
    call jelira(ob2, 'TYPE', cval=typ2)
    if (typ1 .ne. typ2) then
        goto 210
    end if
!
    if (typ1 .eq. 'K') then
        call jelira(ob1, 'LTYP', ltyp1)
        call jelira(ob2, 'LTYP', ltyp2)
        if (ltyp1 .ne. ltyp2) then
            goto 210
        else
            ltyp = ltyp1
        end if
!
!
        if (ltyp .eq. 8) then
            type = 'K8'
        else if (ltyp .eq. 16) then
            type = 'K16'
        else if (ltyp .eq. 24) then
            type = 'K24'
        else if (ltyp .eq. 32) then
            type = 'K32'
        else if (ltyp .eq. 80) then
            type = 'K80'
        end if
    else
        type = typ1
    end if
!
!
!
!   -- determination de xous et genr
!   -------------------------------------
    call jelira(ob1, 'XOUS', cval=xous1)
    call jelira(ob2, 'XOUS', cval=xous2)
    if (xous1 .ne. xous2) then
        goto 210
    else
        xous = xous1
    end if
!
    call jelira(ob1, 'GENR', cval=genr1)
    call jelira(ob2, 'GENR', cval=genr2)
    if (genr1 .ne. genr2) then
        goto 210
    else
        genr = genr1
    end if
!
!
!
!   3- on compare le contenu des objets
!   -------------------------------------
!
!   3.1 : cas des objets simples :
!   ------------------------------
    if (xous .eq. 'S') then
        ASSERT((genr .eq. 'V') .or. (genr .eq. 'N'))
        if (genr .eq. 'V') then
!
            call jelira(ob1, 'LONMAX', l1)
            call jelira(ob2, 'LONMAX', l2)
            if (l1 .ne. l2) then
                goto 210
            end if
!
            call jelira(ob1, 'LONUTI', l1)
            call jelira(ob2, 'LONUTI', l2)
            if (l1 .ne. l2) then
                goto 210
            else
                l = l1
            end if
!
            call jeveuo(ob1, 'L', iad1)
            call jeveuo(ob2, 'L', iad2)
!
            if (type .eq. 'R') then
                do k = 1, l
                    if (zr(iad1-1+k) .ne. zr(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'I') then
                do k = 1, l
                    if (zi(iad1-1+k) .ne. zi(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'C') then
                do k = 1, l
                    if (zc(iad1-1+k) .ne. zc(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'L') then
                do k = 1, l
                    if (.not. (zl(iad1-1+k) .or. (.not. zl(iad2-1+k)))) goto 210
                end do
!
            else if (type .eq. 'K8') then
                do k = 1, l
                    if (zk8(iad1-1+k) .ne. zk8(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'K16') then
                do k = 1, l
                    if (zk16(iad1-1+k) .ne. zk16(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'K24') then
                do k = 1, l
                    if (zk24(iad1-1+k) .ne. zk24(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'K32') then
                do k = 1, l
                    if (zk32(iad1-1+k) .ne. zk32(iad2-1+k)) goto 210
                end do
!
            else if (type .eq. 'K80') then
                do k = 1, l
                    if (zk80(iad1-1+k) .ne. zk80(iad2-1+k)) goto 210
                end do
!
            else
                ASSERT(.false.)
            end if
!
!
        else if (genr .eq. 'N') then
!       ------------------------------
            call jelira(ob1, 'NOMMAX', l1)
            call jelira(ob2, 'NOMMAX', l2)
            if (l1 .ne. l2) then
                goto 210
            end if
!
            call jelira(ob1, 'NOMUTI', l1)
            call jelira(ob2, 'NOMUTI', l2)
            if (l1 .ne. l2) then
                goto 210
            else
                l = l1
            end if
!
            do k = 1, l
                call jenuno(jexnum(ob1, k), k241)
                call jenuno(jexnum(ob2, k), k242)
                if (k241 .ne. k242) goto 210
            end do
!
!
!
        end if
!
!
!   3.2 : cas des collections :
!   ------------------------------
    else
        if (genr .eq. 'V') then
!
            call jelira(ob1, 'NMAXOC', l1)
            call jelira(ob2, 'NMAXOC', l2)
            if (l1 .ne. l2) then
                goto 210
            end if
!
            call jelira(ob1, 'NUTIOC', l1)
            call jelira(ob2, 'NUTIOC', l2)
            if (l1 .ne. l2) then
                goto 210
            end if
            nbobj = l1
!
            do iobj = 1, nbobj
!
                call jelira(jexnum(ob1, iobj), 'LONMAX', l1)
                call jelira(jexnum(ob2, iobj), 'LONMAX', l2)
                if (l1 .ne. l2) then
                    goto 210
                end if
!
                call jelira(jexnum(ob1, iobj), 'LONUTI', l1)
                call jelira(jexnum(ob2, iobj), 'LONUTI', l2)
                if (l1 .ne. l2) then
                    goto 210
                else
                    l = l1
                end if
!
                call jeveuo(jexnum(ob1, iobj), 'L', iad1)
                call jeveuo(jexnum(ob2, iobj), 'L', iad2)
!
                if (type .eq. 'R') then
                    do k = 1, l
                        if (zr(iad1-1+k) .ne. zr(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'I') then
                    do k = 1, l
                        if (zi(iad1-1+k) .ne. zi(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'C') then
                    do k = 1, l
                        if (zc(iad1-1+k) .ne. zc(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'L') then
                    do k = 1, l
                        if (.not. (zl(iad1-1+k) .or. (.not. zl(iad2-1+k)))) goto 210
                    end do
!
                else if (type .eq. 'K8') then
                    do k = 1, l
                        if (zk8(iad1-1+k) .ne. zk8(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'K16') then
                    do k = 1, l
                        if (zk16(iad1-1+k) .ne. zk16(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'K24') then
                    do k = 1, l
                        if (zk24(iad1-1+k) .ne. zk24(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'K32') then
                    do k = 1, l
                        if (zk32(iad1-1+k) .ne. zk32(iad2-1+k)) goto 210
                    end do
!
                else if (type .eq. 'K80') then
                    do k = 1, l
                        if (zk80(iad1-1+k) .ne. zk80(iad2-1+k)) goto 210
                    end do
!
                else
                    ASSERT(.false.)
                end if
            end do
!
        else
            ASSERT(.false.)
        end if
    end if
!
    goto 220
!
210 continue
    idenob = .false.
!
220 continue
!
!
    call jedema()
end function
