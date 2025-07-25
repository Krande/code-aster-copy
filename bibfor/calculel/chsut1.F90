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

subroutine chsut1(chs1, nomgd2, ncmp, lcmp1, lcmp2, &
                  base, chs2)
! person_in_charge: jacques.pellet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/verigd.h"
    integer(kind=8) :: ncmp
    character(len=*) :: chs1, nomgd2, base, chs2
    character(len=8) :: lcmp1(ncmp), lcmp2(ncmp)
! ---------------------------------------------------------------------
! BUT: CHANGER LA GRANDEUR ET LE NOM DES CMPS D'UN CHAMP_S
! ---------------------------------------------------------------------
!     ARGUMENTS:
! CHS1   IN/JXIN  K19 : SD CHAMP_S A MODIFIER
! CHS2   IN/JXOUT K19 : SD CHAMP_S MODIFIEE
! BASE   IN       K1  : /G/V/L
! NCMP   IN       I   : NOMBRE DE CMPS DE LCMP1 ET LCMP2
!                       IL FAUT QUE LCMP1 CONTIENNE TOUTES LES CMPS
!                       DE CHS1
! NOMGD2  IN      K8   : NOM DE LA GRANDEUR "APRES"
! LCMP1   IN      L_K8 : LISTE DES CMPS "AVANT"
! LCMP2   IN      L_K8 : LISTE DES CMPS "APRES"
!
! REMARQUE : CHS2 PEUT ETRE IDENTIQUE A CHS1 (CHAMP_S MODIFIE)
!-----------------------------------------------------------------------
    character(len=24) :: valk(3)
!     ------------------------------------------------------------------
    character(len=19) :: chsa, chsb, chsp
    integer(kind=8) :: i1, i2, jcs1k, jcs1d, jcs1c, jcs2k, jcs2c, k, kk
    integer(kind=8) :: iret, ncmpch
!
    character(len=8) :: nocmp, nomgd1, tsca1, tsca2
!
    chsa = chs1
    chsb = chs2
    chsp = '&&CHUT1.CHAMP_S_IN'
!
    call exisd('CHAM_NO_S', chsa, i1)
    call exisd('CHAM_ELEM_S', chsa, i2)
    if (i1*i2 .ne. 0) then
        call utmess('A', 'CALCULEL2_2', sk=chsa)
    end if
    if (i1+i2 .eq. 0) then
        call utmess('A', 'CALCULEL2_3', sk=chsa)
    end if
!
!
!     1.  ON RECOPIE LE CHAMP "IN" ET ON RECUPERE LES ADRESSES JEVEUX :
!     -----------------------------------------------------------------
    if (i1 .gt. 0) then
!      -- CAS D'UN CHAM_NO_S :
        call copisd('CHAM_NO_S', 'V', chsa, chsp)
        call copisd('CHAM_NO_S', base, chsp, chsb)
        call jeveuo(chsp//'.CNSK', 'L', jcs1k)
        call jeveuo(chsp//'.CNSD', 'L', jcs1d)
        call jeveuo(chsp//'.CNSC', 'L', jcs1c)
        call jeveuo(chsb//'.CNSK', 'E', jcs2k)
        call jeveuo(chsb//'.CNSC', 'E', jcs2c)
!
    else
!      -- CAS D'UN CHAM_ELEM_S :
        call copisd('CHAM_ELEM_S', 'V', chsa, chsp)
        call copisd('CHAM_ELEM_S', base, chsp, chsb)
        call jeveuo(chsp//'.CESK', 'L', jcs1k)
        call jeveuo(chsp//'.CESD', 'L', jcs1d)
        call jeveuo(chsp//'.CESC', 'L', jcs1c)
        call jeveuo(chsb//'.CESK', 'E', jcs2k)
        call jeveuo(chsb//'.CESC', 'E', jcs2c)
    end if
!
!
!     2. QUELQUES VERIFICATIONS :
!     ----------------------------
!
!     2.1 : LES TYPES SCALAIRES DE NOMGD1 ET NOMGD2 SONT LES MEMES:
    nomgd1 = zk8(jcs1k-1+2)
    call dismoi('TYPE_SCA', nomgd1, 'GRANDEUR', repk=tsca1)
    call dismoi('TYPE_SCA', nomgd2, 'GRANDEUR', repk=tsca2)
    if (tsca1 .ne. tsca2) then
        valk(1) = tsca1
        valk(2) = tsca2
        call utmess('F', 'CALCULEL4_4', nk=2, valk=valk)
    end if
!
!     2.2 : NOMGD1 ET LCMP1 SONT COHERENTS :
    call verigd(nomgd1, lcmp1, ncmp, iret)
    ASSERT(iret .le. 0)
!
!     2.3 : NOMGD2 ET LCMP2 SONT COHERENTS :
    call verigd(nomgd2, lcmp2, ncmp, iret)
    ASSERT(iret .le. 0)
!
!
!      3. MODIFICATION DE CHS2 :
!      -------------------------
    zk8(jcs2k-1+2) = nomgd2
    ncmpch = zi(jcs1d-1+2)
    do k = 1, ncmpch
        nocmp = zk8(jcs1c-1+k)
        kk = indik8(lcmp1, nocmp, 1, ncmp)
!       SI KK.EQ.0 : ON NE SAIT PAS RENOMMER LA CMP
        ASSERT(kk .ne. 0)
        zk8(jcs2c-1+k) = lcmp2(kk)
    end do
!
!
!
!     5. MENAGE :
!     -----------
    if (i1 .gt. 0) then
        call detrsd('CHAM_NO_S', chsp)
    else
        call detrsd('CHAM_ELEM_S', chsp)
    end if
!
end subroutine
