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

subroutine cnsdot(cns1z, cns2z, pscal, ier)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: cns1z, cns2z
    real(kind=8) :: pscal
    integer(kind=8) :: ier
! ---------------------------------------------------------------------
! BUT: CALCULER LE "PRODUIT SCALAIRE" DE 2 CHAM_NO_S
! ---------------------------------------------------------------------
!     ARGUMENTS:
! CNS1Z  IN/JXIN  K19 : SD CHAM_NO_S
! CNS2Z  IN/JXIN  K19 : SD CHAM_NO_S
! PSCAL  OUT      R   : "PRODUIT SCALAIRE"
! IER    OUT      I   : CODE_RETOUR :
!                       0 -> OK
!                       1 -> LES 2 CHAM_NO_S N'ONT PAS LES MEMES CMPS
!
! REMARQUE :
! CETTE ROUTINE BOUCLE SUR TOUS LES DDLS DE CNS1 ET CNS2PSCAL=0
!  POUR IN IN CNS1 :
!      CMP1=CNS2
!-----------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    integer(kind=8) ::    jcnsl1
    integer(kind=8) ::    jcnsl2
    integer(kind=8) :: nbno, k, ino, ncmp, nbno1, nbno2, ncmp1, ncmp2
    character(len=8) :: ma1, nomgd1, ma2, nomgd2
    character(len=3) :: tsca1
    character(len=19) :: cns1, cns2
    real(kind=8) :: x1, x2
    character(len=8), pointer :: cnsk1(:) => null()
    character(len=8), pointer :: cnsk2(:) => null()
    character(len=8), pointer :: cnsc1(:) => null()
    character(len=8), pointer :: cnsc2(:) => null()
    real(kind=8), pointer :: cnsv1(:) => null()
    real(kind=8), pointer :: cnsv2(:) => null()
    integer(kind=8), pointer :: cnsd1(:) => null()
    integer(kind=8), pointer :: cnsd2(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
!
    cns1 = cns1z
    cns2 = cns2z
!
    call jeveuo(cns1//'.CNSK', 'L', vk8=cnsk1)
    call jeveuo(cns1//'.CNSD', 'L', vi=cnsd1)
    call jeveuo(cns1//'.CNSC', 'L', vk8=cnsc1)
    call jeveuo(cns1//'.CNSV', 'L', vr=cnsv1)
    call jeveuo(cns1//'.CNSL', 'L', jcnsl1)
!
    ma1 = cnsk1(1)
    nomgd1 = cnsk1(2)
    nbno1 = cnsd1(1)
    ncmp1 = cnsd1(2)
!
!
    call jeveuo(cns2//'.CNSK', 'L', vk8=cnsk2)
    call jeveuo(cns2//'.CNSD', 'L', vi=cnsd2)
    call jeveuo(cns2//'.CNSC', 'L', vk8=cnsc2)
    call jeveuo(cns2//'.CNSV', 'L', vr=cnsv2)
    call jeveuo(cns2//'.CNSL', 'L', jcnsl2)
!
    ma2 = cnsk2(1)
    nomgd2 = cnsk2(2)
    nbno2 = cnsd2(1)
    ncmp2 = cnsd2(2)
!
!     -- COHERENCE DES 2 CHAMPS :
    ASSERT(ma1 .eq. ma2)
    ASSERT(nomgd1 .eq. nomgd2)
    ASSERT(nbno1 .eq. nbno2)
    ASSERT(ncmp1 .eq. ncmp2)
    nbno = nbno1
    ncmp = ncmp1
!
!
!
    call dismoi('TYPE_SCA', nomgd1, 'GRANDEUR', repk=tsca1)
    ASSERT(tsca1 .eq. 'R')
!
    do k = 1, ncmp
        ASSERT(cnsc1(k) .eq. cnsc2(k))
    end do
!
!
!     CALCUL DE LA SOMME DES PRODUITS DES CMPS :
!     -------------------------------------------
    pscal = 0.d0
    ier = 0
    do ino = 1, nbno
        do k = 1, ncmp
            if (zl(jcnsl1-1+(ino-1)*ncmp+k)) then
                if (.not. zl(jcnsl2-1+(ino-1)*ncmp+k)) then
                    ier = 1
                    goto 40
!
                end if
                x1 = cnsv1((ino-1)*ncmp+k)
                x2 = cnsv2((ino-1)*ncmp+k)
                pscal = pscal+x1*x2
            else
                if (zl(jcnsl2-1+(ino-1)*ncmp+k)) then
                    ier = 1
                    goto 40
!
                end if
            end if
        end do
    end do
!
40  continue
    call jedema()
end subroutine
