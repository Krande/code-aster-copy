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

subroutine rdtcns(ma2, corrn, cns1, base, cns2)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cnscre.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=8) :: ma2
    character(len=19) :: cns1, cns2
    character(len=*) :: corrn
    character(len=1) :: base
! ---------------------------------------------------------------------
! BUT: REDUIRE UN CHAM_NO_S SUR UN MAILLAGE REDUIT
! ---------------------------------------------------------------------
! ARGUMENTS:
! MA2    IN       K8  : MAILLAGE REDUIT
! CNS1   IN/JXIN  K19 : CHAM_NO_S A REDUIRE
! CNS2   IN/JXOUT K19 : CHAM_NO_S REDUIT
! BASE   IN       K1  : BASE DE CREATION POUR CNS2Z : G/V/L
! CORRN  IN       K*  : NOM DE L'OBJET CONTENANT LA CORRESPONDANCE
!                       INO_RE -> INO
!
!-----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) :: nbno1, nbno2, jcorrn
    integer(kind=8) ::   jcn1v, jcn1l
    integer(kind=8) ::  jcn2v, jcn2l, jcn2c
    integer(kind=8) :: ncmp, icmp
    integer(kind=8) :: ino2, ino1
    character(len=8) :: nomgd
    character(len=3) :: tsca
    integer(kind=8), pointer :: cn1d(:) => null()
    integer(kind=8), pointer :: cn2d(:) => null()
    character(len=8), pointer :: cn1c(:) => null()
    character(len=8), pointer :: cnsk(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    ASSERT(cns2 .ne. ' ')
    ASSERT(cns1 .ne. cns2)
!
    call jeveuo(cns1//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(cns1//'.CNSD', 'L', vi=cn1d)
    call jeveuo(cns1//'.CNSC', 'L', vk8=cn1c)
    call jeveuo(cns1//'.CNSV', 'L', jcn1v)
    call jeveuo(cns1//'.CNSL', 'L', jcn1l)
!
    nomgd = cnsk(2)
    nbno1 = cn1d(1)
    ncmp = cn1d(2)
    ASSERT(ncmp .gt. 0)
!
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
!
!
!
!     1- CREATION DE CNS2 :
!     ---------------------------------------
    call cnscre(ma2, nomgd, ncmp, cn1c, base, &
                cns2)
    call jeveuo(cns2//'.CNSD', 'L', vi=cn2d)
    call jeveuo(cns2//'.CNSC', 'L', jcn2c)
    call jeveuo(cns2//'.CNSV', 'E', jcn2v)
    call jeveuo(cns2//'.CNSL', 'E', jcn2l)
    nbno2 = cn2d(1)
    ASSERT(nbno2 .gt. 0)
    ASSERT(nbno2 .le. nbno1)
!
!
!     3- REMPLISSAGE DES OBJETS .CNSL ET .CNSV :
!     ------------------------------------------
    call jeveuo(corrn, 'L', jcorrn)
!
    do icmp = 1, ncmp
!
        do ino2 = 1, nbno2
            ino1 = zi(jcorrn-1+ino2)
            if (zl(jcn1l-1+(ino1-1)*ncmp+icmp)) then
                zl(jcn2l-1+(ino2-1)*ncmp+icmp) = .true.
!
                if (tsca .eq. 'R') then
                    zr(jcn2v-1+(ino2-1)*ncmp+icmp) = zr(jcn1v-1+(ino1-1) &
                                                        *ncmp+icmp)
                else if (tsca .eq. 'C') then
                    zc(jcn2v-1+(ino2-1)*ncmp+icmp) = zc(jcn1v-1+(ino1-1) &
                                                        *ncmp+icmp)
                else if (tsca .eq. 'I') then
                    zi(jcn2v-1+(ino2-1)*ncmp+icmp) = zi(jcn1v-1+(ino1-1) &
                                                        *ncmp+icmp)
                else if (tsca .eq. 'L') then
                    zl(jcn2v-1+(ino2-1)*ncmp+icmp) = zl(jcn1v-1+(ino1-1) &
                                                        *ncmp+icmp)
                else if (tsca .eq. 'K8') then
                    zk8(jcn2v-1+(ino2-1)*ncmp+icmp) = zk8(jcn1v-1+(ino1- &
                                                                   1)*ncmp+icmp)
                else
                    ASSERT(.false.)
                end if
!
            end if
!
        end do
    end do
!
!
    call jedema()
end subroutine
