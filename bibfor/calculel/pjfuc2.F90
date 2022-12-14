! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine pjfuc2(c1, c2, base, c3)
! person_in_charge: nicolas.greffet at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedupo.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=16) :: c1, c2, c3
    character(len=1) :: base
! ======================================================================
!     BUT : VARIANTE DE PJFUCO POUR LE COUPLAGE IFS VIA YACS
!       FUSIONNER 2 SD CORRESP_2_MAILLA C1 ET C2 POUR FABRIQUER C3
!
!  IN/JXIN   C1  : SD CORRESP_2_MAILLA
!  IN/JXIN   C2  : SD CORRESP_2_MAILLA
!  IN/JXOUT  C3  : SD CORRESP_2_MAILLA RESULTAT DE LA FUSION
!  IN        BASE: NOM DE LA BASE POUR CREER C3
!
!----------------------------------------------------------------------
!
!
    integer :: ii, ino, ilengt, ideca1, ideca2
    integer :: jno1, nbno1
    integer :: jno2, nbno2
    integer :: jnb3, jm13, jcf3, jnu3
    character(len=8) :: ma1, ma2
    character(len=24) :: valk(2)
    integer, pointer :: nu1(:) => null()
    integer, pointer :: nu2(:) => null()
    integer, pointer :: nb1(:) => null()
    integer, pointer :: nb2(:) => null()
    integer, pointer :: m11(:) => null()
    integer, pointer :: m12(:) => null()
    real(kind=8), pointer :: cf1(:) => null()
    real(kind=8), pointer :: cf2(:) => null()
!
! DEB ------------------------------------------------------------------
    call jemarq()
!
!
!     1- OBJET '.PJEF_NO' :
!     ----------------------
!
! PJEF_NO : DEVIENT PJXX_K1 DEPUIS 10.1.9
!      CALL JEVEUO(C1//'.PJEF_NO','L',JNO1)
!      CALL JEVEUO(C2//'.PJEF_NO','L',JNO2)
    call jeveuo(c1//'.PJXX_K1', 'L', jno1)
    call jeveuo(c2//'.PJXX_K1', 'L', jno2)
    do ii = 1, 2
        ma1 = zk24(jno1-1+ii)(1:8)
        ma2 = zk24(jno2-1+ii)(1:8)
        if (ma1 .ne. ma2) then
            valk(1) = ma1
            valk(2) = ma2
            call utmess('F', 'CALCULEL4_65', nk=2, valk=valk)
        endif
    end do
    call jedupo(c1//'.PJXX_K1', base, c3//'.PJXX_K1', .false._1)
!
!     2- RECUPERATION DES POINTEURS
!     -----------------------------
    call jelira(c1//'.PJEF_NB', 'LONMAX', nbno1)
    call jelira(c2//'.PJEF_NB', 'LONMAX', nbno2)
    call jeveuo(c1//'.PJEF_NB', 'L', vi=nb1)
    call jeveuo(c2//'.PJEF_NB', 'L', vi=nb2)
    call jeveuo(c1//'.PJEF_M1', 'L', vi=m11)
    call jeveuo(c2//'.PJEF_M1', 'L', vi=m12)
    call jeveuo(c1//'.PJEF_CF', 'L', vr=cf1)
    call jeveuo(c2//'.PJEF_CF', 'L', vr=cf2)
    call jeveuo(c1//'.PJEF_NU', 'L', vi=nu1)
    call jeveuo(c2//'.PJEF_NU', 'L', vi=nu2)
!
!     3- AFFECTATION DE PJEF_NB ET PJEF_M1
!     ------------------------------------
    call wkvect(c3//'.PJEF_NB', base//' V I', nbno1+nbno2, jnb3)
    call wkvect(c3//'.PJEF_M1', base//' V I', nbno1+nbno2, jm13)
!
    ilengt = 0
    do ino = 1, nbno1
        zi(jnb3-1+ino) = nb1(ino)
        zi(jm13-1+ino) = m11(ino)
        ilengt = ilengt + nb1(ino)
    end do
    do ino = 1, nbno2
        zi(jnb3-1+nbno1+ino) = nb2(ino)
        zi(jm13-1+nbno1+ino) = m12(ino)
        ilengt = ilengt + nb2(ino)
    end do
!
!     4 - AFFECTATION DE PJEF_CF ET PJEF_NU
!     -------------------------------------
    call wkvect(c3//'.PJEF_CF', base//' V R', ilengt, jcf3)
    call wkvect(c3//'.PJEF_NU', base//' V I', ilengt, jnu3)
!
    ideca1 = 0
    do ino = 1, nbno1
        do ii = 1, nb1(ino)
            zr(jcf3-1+ideca1+ii) = cf1(ideca1+ii)
            zi(jnu3-1+ideca1+ii) = nu1(ideca1+ii)
        end do
        ideca1 = ideca1 + nb1(ino)
    end do
    ideca2 = 0
    do ino = 1, nbno2
        do ii = 1, nb2(ino)
            zr(jcf3-1+ideca1+ii) = cf2(ideca2+ii)
            zi(jnu3-1+ideca1+ii) = nu2(ideca2+ii)
        end do
        ideca1 = ideca1 + nb2(ino)
        ideca2 = ideca2 + nb2(ino)
    end do
!
!     5 - LIBERATION DE LA MEMOIRE
!     ----------------------------
    call jedema()
!
! FIN ------------------------------------------------------------------
end subroutine
