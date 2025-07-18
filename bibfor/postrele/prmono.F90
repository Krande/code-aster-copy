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
subroutine prmono(champ, ioc, som, nbcmp, nocmp)
!     COMMANDE : POST_RELEVE_T
!                DETERMINE LA MOYENNE SUR DES ENTITES POUR UN CHAM_NO
!
! ----------------------------------------------------------------------
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/cnocns.h"
#include "asterfort/detrsd.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
    integer(kind=8) :: ioc, nbcmp
    real(kind=8) :: som(1)
    character(len=8) :: nocmp(1)
    character(len=*) :: champ
!
    integer(kind=8) :: jcnsl, nbno, ncmp, nbn
    integer(kind=8) :: ibid, nbnoeu, idnoeu, nbc
    integer(kind=8) :: i100, i110, icp, ino
    real(kind=8) :: x
    character(len=8) :: ma
    character(len=16) :: motcle(4), typmcl(4)
    character(len=19) :: chams1
    character(len=24) :: mesnoe
    character(len=8), pointer :: nom_cmp(:) => null()
    character(len=8), pointer :: cnsc(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
    real(kind=8), pointer :: cnsv(:) => null()
    integer(kind=8), pointer :: cnsd(:) => null()
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
    mesnoe = '&&PRMONO.MES_NOEUDS'
!
    chams1 = '&&PRMONO.CHAMS1'
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
    do i100 = 1, nbcmp
        if (nbc .ne. 0) then
            nocmp(i100) = nom_cmp(i100)
            icp = indik8(cnsc, nocmp(i100), 1, ncmp)
            if (icp .eq. 0) goto 100
        else
            icp = i100
            nocmp(i100) = cnsc(i100)
        end if
        som(i100) = 0.d0
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
                som(i100) = som(i100)+x
!
            end if
!
        end do
        som(i100) = som(i100)/nbnoeu
!
100     continue
    end do
!
! --- MENAGE
    call detrsd('CHAM_NO_S', chams1)
    call jedetr(mesnoe)
    AS_DEALLOCATE(vk8=nom_cmp)
!
end subroutine
