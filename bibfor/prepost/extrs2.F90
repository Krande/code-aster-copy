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

subroutine extrs2(resu0, resu1, typcon, lrest, mailla, &
                  modele, cara, chmat, nbordr, nuordr, nbacc, nomacc, &
                  nbarch, nuarch, nbexcl, chexcl, nbnosy)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/rdtchp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: nbordr, nuordr(*), nbarch, nbacc, nuarch(*), nbexcl, nbnosy
    character(len=16) :: nomacc(*), chexcl(*)
    character(len=*) :: resu0, resu1
    character(len=16) :: typcon
    character(len=8), intent(in) :: mailla
    character(len=8), intent(in) :: modele
    character(len=8), intent(in) :: cara
    character(len=8), intent(in) :: chmat
    aster_logical :: lrest
! person_in_charge: nicolas.sellenet at edf.fr
!     OPERATEUR D'EXTRACTION
!     ------------------------------------------------------------------
!
!
! 0.3. ==> VARIABLES LOCALES
!
    integer(kind=8) :: vali(2)
!
    integer(kind=8) :: i, j, ire1, ire2, iadin, iadou, iret
    integer(kind=8) :: cret
    character(len=3) :: type
    character(len=4) :: tych
    character(len=8) :: noma1, noma2, noma3, nomavr
    character(len=16) :: nomsym
    character(len=16) :: nopara
    character(len=19) :: resuin, resuou, ligrel
    character(len=24) :: chamin, chamou, corrn, corrm
    character(len=24) :: valk(3)
    character(len=8), pointer :: lgrf(:) => null()
    character(len=8), pointer :: maor(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    resuin = ' '
    resuou = '  '
    i = lxlgut(resu0)
    resuin(1:i) = resu0(1:i)
    i = lxlgut(resu1)
    resuou(1:i) = resu1(1:i)
!
!
    call jeexin(resuou//'.DESC', iret)
    if (iret .eq. 0) then
        call rscrsd('G', resuou, typcon, nbarch)
    end if
!
    if (lrest) then
        if (modele .ne. ' ') then
            call jeveuo(modele//'.MODELE    .LGRF', 'L', vk8=lgrf)
            noma2 = lgrf(1)
            ligrel = modele//'.MODELE'
        else
            noma2 = mailla
            ligrel = ' '
        end if
        call jeveuo(noma2//'.MAOR', 'L', vk8=maor)
        noma1 = maor(1)
        corrn = noma2//'.CRNO'
        corrm = noma2//'.CRMA'

!       -- si cara et chmat sont fournis, on verifie qu'ils s'appuient sur le
!          "bon" maillage :
        if (cara .ne. ' ') then
            call dismoi('NOM_MAILLA', cara, 'CARA_ELEM', repk=noma3)
            if (noma3 .ne. noma2) then
                valk(1) = cara
                valk(2) = noma2
                valk(3) = noma3
                call utmess('F', 'CALCULEL3_48', nk=3, valk=valk)
            end if
        end if
        if (chmat .ne. ' ') then
            call dismoi('NOM_MAILLA', chmat, 'CHAM_MATER', repk=noma3)
            if (noma3 .ne. noma2) then
                valk(1) = chmat
                valk(2) = noma2
                valk(3) = noma3
                call utmess('F', 'CALCULEL3_49', nk=3, valk=valk)
            end if
        end if
    end if
!
    do i = 1, nbnosy
!
        call jenuno(jexnum(resuin//'.DESC', i), nomsym)
        if (nomsym .eq. 'COMPORTEMENT') then
        else
            do j = 1, nbexcl
                if (chexcl(j) .eq. nomsym) goto 30
            end do
        end if
        do j = 1, nbordr
            if (nuarch(j) .eq. 0) goto 20
            call rsexch(' ', resuin, nomsym, nuordr(j), chamin, &
                        ire1)
            if (ire1 .gt. 0) goto 20
!
            call rsexch(' ', resuou, nomsym, nuordr(j), chamou, &
                        ire2)
            if (ire2 .eq. 0) then
            else if (ire2 .eq. 100) then
            else
                vali(1) = nuordr(j)
                vali(2) = ire2
                valk(1) = chamou
                call utmess('F', 'PREPOST5_16', sk=valk(1), ni=2, vali=vali)
            end if
            if (lrest) then
                call dismoi('NOM_MAILLA', chamin, 'CHAMP', repk=nomavr)
                if (noma1 .ne. nomavr) then
                    valk(1) = resuin
                    valk(2) = nomavr
                    valk(3) = noma1
                    call utmess('F', 'PREPOST5_49', nk=3, valk=valk)
                end if
                call dismoi('TYPE_CHAMP', chamin, 'CHAMP', repk=tych)
                if (tych(1:2) .eq. 'EL') then
                    if (ligrel .eq. ' ') call utmess('F', 'PREPOST5_48')
                end if
                call rdtchp(corrn, corrm, chamin(1:19), chamou(1:19), 'G', &
                            noma1, noma2, ligrel, cret)
                if (cret .ne. 0) goto 20
            else
                call copisd('CHAMP_GD', 'G', chamin, chamou)
            end if
            call rsnoch(resuou, nomsym, nuordr(j))
20          continue
        end do
30      continue
    end do
!
!
    do i = 1, nbordr
        if (nuarch(i) .eq. 0) goto 50
        do j = 1, nbacc
            nopara = nomacc(j)
            call rsadpa(resuin, 'L', 1, nopara, nuordr(i), &
                        1, sjv=iadin, styp=type, istop=0)
            call rsadpa(resuou, 'E', 1, nopara, nuordr(i), &
                        1, sjv=iadou, styp=type)
            if (type(1:1) .eq. 'I') then
                zi(iadou) = zi(iadin)
            else if (type(1:1) .eq. 'R') then
                zr(iadou) = zr(iadin)
            else if (type(1:1) .eq. 'C') then
                zc(iadou) = zc(iadin)
            else if (type(1:3) .eq. 'K80') then
                zk80(iadou) = zk80(iadin)
            else if (type(1:3) .eq. 'K32') then
                zk32(iadou) = zk32(iadin)
            else if (type(1:3) .eq. 'K24') then
                zk24(iadou) = zk24(iadin)
                if (nopara(1:5) .eq. 'EXCIT' .and. zk24(iadin) (1:2) .ne. '  ') then
                    if (.not. lrest) then
                        zk24(iadou) = resuou(1:8)//zk24(iadin) (9:)
                        call copisd(' ', 'G', zk24(iadin) (1:19), zk24(iadou) (1:19))
                    end if
                end if
            else if (type(1:3) .eq. 'K16') then
                zk16(iadou) = zk16(iadin)
            else if (type(1:2) .eq. 'K8') then
                zk8(iadou) = zk8(iadin)
            end if
        end do

!       -- pour RESTREINT, on surcharge les parametres : MODELE, CARA_ELEM, CHAM_MATER et EXCIT:
        if (lrest) then
            call rsadpa(resuou, 'E', 1, 'MODELE', nuordr(i), 1, sjv=iadou, styp=type)
            zk8(iadou) = modele
            call rsadpa(resuou, 'E', 1, 'CHAMPMAT', nuordr(i), 1, sjv=iadou, styp=type)
            zk8(iadou) = chmat
            call rsadpa(resuou, 'E', 1, 'CARAELEM', nuordr(i), 1, sjv=iadou, styp=type)
            zk8(iadou) = cara
            call rsadpa(resuou, 'E', 1, 'EXCIT', nuordr(i), 1, sjv=iadou, styp=type)
            zk24(iadou) = ' '
        end if

50      continue
    end do

    call jedema()

end subroutine
