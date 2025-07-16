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
subroutine rc36rs(nomres, noma, nbma, listma, chindi, &
                  chresu)
    implicit none
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nbma, listma(*)
    character(len=8) :: nomres, noma
    character(len=24) :: chindi, chresu
!
!     OPERATEUR POST_RCCM, TRAITEMENT DE FATIGUE_B3600
!
!     TRANSFERT DE CHRESU DANS LA TABLE
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ibid, im, ima, decal, nbcmp, nbpt, ipt, ino, jconx1
    integer(kind=8) :: jconx2, npara, nbcin, decin, icmp, iad
    parameter(npara=8)
    real(kind=8) :: valer(5), type
    complex(kind=8) :: c16b
    character(len=8) :: valek(3), typara(npara)
    character(len=16) :: nopara(npara)
    character(len=24) :: connex
    integer(kind=8), pointer :: cesd(:) => null()
    integer(kind=8), pointer :: cind(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: cinv(:) => null()
!     ------------------------------------------------------------------
    data nopara/'MAILLE', 'TYPE_MAILLE', 'NOEUD', 'SM',&
     &              'SN_MAX', 'SN/3SM', 'SALT_MAX', 'FACT_USAGE_CUMU'/
    data typara/'K8', 'K8', 'K8', 'R',&
     &              'R', 'R', 'R', 'R'/
! DEB ------------------------------------------------------------------
!
    ibid = 0
    c16b = (0.d0, 0.d0)
    call tbcrsd(nomres, 'G')
    call tbajpa(nomres, npara, nopara, typara)
!
    connex = noma//'.CONNEX         '
    call jeveuo(connex, 'L', jconx1)
    call jeveuo(jexatr(connex, 'LONCUM'), 'L', jconx2)
!
! --- LE CHAMP INDICE DE CONTRAINTES
!
    call jeveuo(chindi(1:19)//'.CESV', 'L', vr=cinv)
    call jeveuo(chindi(1:19)//'.CESD', 'L', vi=cind)
    nbcin = cind(2)
!
! --- LE CHAM_ELEM RESULTAT
!
    call jeveuo(chresu(1:19)//'.CESD', 'L', vi=cesd)
    call jeveuo(chresu(1:19)//'.CESV', 'L', vr=cesv)
    nbcmp = cesd(2)
!
    do im = 1, nbma
!
        ima = listma(im)
        valek(1) = int_to_char8(ima)
!
        nbpt = cesd(5+4*(ima-1)+1)
        decal = cesd(5+4*(ima-1)+4)
        decin = cind(5+4*(ima-1)+4)
!
        do ipt = 1, nbpt
!
            icmp = 7
            iad = decin+(ipt-1)*nbcin+icmp
            type = cinv(iad)
            if (type .eq. 0.d0) then
                valek(2) = '???'
            else if (type .eq. 10.d0) then
                valek(2) = 'DRO'
            else if (type .eq. 20.d0) then
                valek(2) = 'COU'
            else if (type .eq. 30.d0) then
                valek(2) = 'TRN'
            else if (type .eq. 40.d0) then
                valek(2) = 'TEE'
            end if
!
            ino = zi(jconx1-1+zi(jconx2+ima-1)+ipt-1)
            valek(3) = int_to_char8(ino)
!
            do icmp = 1, nbcmp
!
                iad = decal+(ipt-1)*nbcmp+icmp
                valer(icmp) = cesv(iad)
!
            end do
!
            call tbajli(nomres, npara, nopara, [ibid], valer, &
                        [c16b], valek, 0)
!
        end do
!
    end do
!
end subroutine
