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
subroutine prexno(champ, ioc, nomax, cmpmax, valmax, &
                  nomin, cmpmin, valmin, noamax, cmamax, &
                  vaamax, noamin, cmamin, vaamin)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8vide.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: ioc
    real(kind=8) :: valmin, valmax, vaamin, vaamax
    character(len=8) :: nomax, nomin, cmpmax, cmpmin
    character(len=8) :: noamax, noamin, cmamax, cmamin
    character(len=*) :: champ
!
!     COMMANDE : POST_RELEVE_T
!                DETERMINE LES EXTREMA POUR UN CHAM_NO
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: jcnsl, nbno, ncmp, nbn
    integer(kind=8) :: ibid, nbnoeu, idnoeu, nbc, nbcmp
    integer(kind=8) :: i100, i110, icp, ino, inomax, inomin, inamax, inamin
    real(kind=8) :: x
    character(len=8) :: nocmp, ma
    character(len=16) :: motcle(4), typmcl(4)
    character(len=19) :: chams1
    character(len=24) :: mesnoe
    character(len=8), pointer :: nom_cmp(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
! ---------------------------------------------------------------------
!
    motcle(1) = 'GROUP_NO'
    motcle(2) = 'NOEUD'
    motcle(3) = 'GROUP_MA'
    motcle(4) = 'MAILLE'
    typmcl(1) = 'GROUP_NO'
    typmcl(2) = 'NOEUD'
    typmcl(3) = 'GROUP_MA'
    typmcl(4) = 'MAILLE'
    mesnoe = '&&PREXNO.MES_NOEUDS'
!
    chams1 = '&&PREXNO.CHAMS1'
    call cnocns(champ, 'V', chams1)
!
    call jeveuo(chams1//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(chams1//'.CNSD', 'L', vi=cnsd)
    call jeveuo(chams1//'.CNSC', 'L', vk8=cnsc)
    call jeveuo(chams1//'.CNSV', 'L', vr=cnsv)
    call jeveuo(chams1//'.CNSL', 'L', jcnsl)
    ma = cnsk(1)
    nbno = cnsd(1)
    ncmp = cnsd(2)
!
    call reliem(' ', ma, 'NU_NOEUD', 'ACTION', ioc, &
                4, motcle, typmcl, mesnoe, nbn)
    if (nbn .gt. 0) then
        nbnoeu = nbn
        call jeveuo(mesnoe, 'L', idnoeu)
    else
        nbnoeu = nbno
    end if
!
    call getvtx('ACTION', 'NOM_CMP', iocc=ioc, nbval=0, nbret=nbc)
    if (nbc .ne. 0) then
        nbcmp = -nbc
        AS_ALLOCATE(vk8=nom_cmp, size=nbcmp)
        call getvtx('ACTION', 'NOM_CMP', iocc=ioc, nbval=nbcmp, vect=nom_cmp, &
                    nbret=ibid)
    else
        nbcmp = ncmp
    end if
!
    inomax = 0
    valmax = -r8vide()
    inomin = 0
    valmin = r8vide()
!
    inamax = 0
    vaamax = -r8vide()
    inamin = 0
    vaamin = r8vide()
!
    do i100 = 1, nbcmp
        if (nbc .ne. 0) then
            nocmp = nom_cmp(i100)
            icp = indik8(cnsc, nocmp, 1, ncmp)
            if (icp .eq. 0) goto 100
        else
            icp = i100
            nocmp = cnsc(i100)
        end if
!
        do i110 = 1, nbnoeu
            if (nbn .gt. 0) then
                ino = zi(idnoeu+i110-1)
            else
                ino = i110
            end if
!
            if (zl(jcnsl-1+(ino-1)*ncmp+icp)) then
!
                x = cnsv((ino-1)*ncmp+icp)
!
                if (x .gt. valmax) then
                    inomax = ino
                    valmax = x
                    cmpmax = nocmp
                end if
!
                if (abs(x) .gt. vaamax) then
                    inamax = ino
                    vaamax = abs(x)
                    cmamax = nocmp
                end if
!
                if (x .lt. valmin) then
                    inomin = ino
                    valmin = x
                    cmpmin = nocmp
                end if
!
                if (abs(x) .lt. vaamin) then
                    inamin = ino
                    vaamin = abs(x)
                    cmamin = nocmp
                end if
!
            end if
!
        end do
!
100     continue
    end do
!
    if (inomax .eq. 0) call utmess('F', 'POSTRELE_18')
    nomax = int_to_char8(inomax)
    nomin = int_to_char8(inomin)
    noamax = int_to_char8(inamax)
    noamin = int_to_char8(inamin)
!
! --- MENAGE
    call detrsd('CHAM_NO_S', chams1)
    call jedetr(mesnoe)
    AS_DEALLOCATE(vk8=nom_cmp)
!
end subroutine
