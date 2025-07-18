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
subroutine tbcopi(base, sd1, sd2)
    implicit none
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: base, sd1, sd2
!
!   BUT:
!   DUPLIQUER UNE STRUCTURE DE DONNEES TABLE.
!
!     IN:
!     BASE     : 'G' , 'V' , ... : BASE DE CREATION DE SD2
!     SD1 (K*) : NOM DE LA SD A DUPPLIQUER
!     SD2 (K*) : NOM DE LA SD A CREER
!
!     OUT:
!     SD2 EST CREEE ET A LE MEME CONTENU QUE SD1
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: i, j, nbpm, nbpu, jnjv, knjv, kvale, jvale, nbpara, jtbba, ktbnp
    integer(kind=8) :: nblign, jtbnp, ndim, jtblp, ktblp
    character(len=1) :: bas2
    character(len=4) :: type, knume
    character(len=19) :: tab1, tab2
    character(len=24) :: nomjv
!
! DEB-------------------------------------------------------------------
!
    call jemarq()
    bas2 = base
!
    tab1 = sd1
    tab2 = sd2
!
    call wkvect(tab2//'.TBBA', bas2//' V K8', 1, jtbba)
    zk8(jtbba) = bas2
    call jeveuo(tab1//'.TBNP', 'L', ktbnp)
    nbpara = zi(ktbnp)
    nblign = zi(ktbnp+1)
!
    call wkvect(tab2//'.TBNP', bas2//' V I', 2, jtbnp)
    zi(jtbnp) = nbpara
    zi(jtbnp+1) = nblign
    ndim = 4*nbpara
    call jecreo(tab2//'.TBLP', bas2//' V K24')
    call jeecra(tab2//'.TBLP', 'LONMAX', ndim)
    call jeecra(tab2//'.TBLP', 'LONUTI', ndim)
    call jeveuo(tab2//'.TBLP', 'E', jtblp)
    call jeveuo(tab1//'.TBLP', 'L', ktblp)
    do i = 1, nbpara
        zk24(jtblp+4*(i-1)) = zk24(ktblp+4*(i-1))
        type = zk24(ktblp+4*(i-1)+1)
        zk24(jtblp+4*(i-1)+1) = type
        nomjv = zk24(ktblp+4*(i-1)+2)
        call jelira(nomjv, 'LONMAX', nbpm)
        call jelira(nomjv, 'LONUTI', nbpu)
        call codent(i, 'D0', knume)
        nomjv = tab2(1:17)//'LG.'//knume
        zk24(jtblp+4*(i-1)+3) = nomjv
        call jecreo(nomjv, bas2//' V I')
        call jeecra(nomjv, 'LONMAX', nbpm)
        call jeecra(nomjv, 'LONUTI', nbpu)
        call jeveuo(nomjv, 'E', jnjv)
        nomjv = tab1(1:17)//'LG.'//knume
        call jeveuo(nomjv, 'L', knjv)
        do j = 1, nbpm
            zi(jnjv+j-1) = zi(knjv+j-1)
        end do
        nomjv = tab1//'.'//knume
        call jeveuo(nomjv, 'L', kvale)
        nomjv = tab2//'.'//knume
        zk24(jtblp+4*(i-1)+2) = nomjv
        call jecreo(nomjv, bas2//' V '//type)
        call jeecra(nomjv, 'LONMAX', nbpm)
        call jeecra(nomjv, 'LONUTI', nbpu)
        call jeveuo(nomjv, 'E', jvale)
        if (type(1:1) .eq. 'I') then
            do j = 1, nbpm
                zi(jvale+j-1) = zi(kvale+j-1)
            end do
        else if (type(1:1) .eq. 'R') then
            do j = 1, nbpm
                zr(jvale+j-1) = zr(kvale+j-1)
            end do
        else if (type(1:1) .eq. 'C') then
            do j = 1, nbpm
                zc(jvale+j-1) = zc(kvale+j-1)
            end do
        else if (type(1:3) .eq. 'K80') then
            do j = 1, nbpm
                zk80(jvale+j-1) = zk80(kvale+j-1)
            end do
        else if (type(1:3) .eq. 'K32') then
            do j = 1, nbpm
                zk32(jvale+j-1) = zk32(kvale+j-1)
            end do
        else if (type(1:3) .eq. 'K24') then
            do j = 1, nbpm
                zk24(jvale+j-1) = zk24(kvale+j-1)
            end do
        else if (type(1:3) .eq. 'K16') then
            do j = 1, nbpm
                zk16(jvale+j-1) = zk16(kvale+j-1)
            end do
        else if (type(1:2) .eq. 'K8') then
            do j = 1, nbpm
                zk8(jvale+j-1) = zk8(kvale+j-1)
            end do
        end if
        call jeecra(nomjv, 'LONUTI', nbpu)
    end do
!
    call jedema()
end subroutine
