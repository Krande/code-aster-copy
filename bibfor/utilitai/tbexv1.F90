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
subroutine tbexv1(nomta, para, nomobj, basobj, nbval, &
                  typval)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nbval
    character(len=*) :: nomta, para, nomobj, basobj, typval
!      LECTURE DES VALEURS D'UNE COLONNE D'UNE TABLE
!              EN ELIMINANT LES DOUBLONS.
! ----------------------------------------------------------------------
! IN  : NOMTA  : NOM DE LA STRUCTURE "TABLE".
! IN  : PARA   : PARAMETRE DESIGNANT LA COLONNE A EXTRAIRE
! IN  : NOMOBJ : NOM DE L'OBJET JEVEUX CONTENANT LES VALEURS
! IN  : BASOBJ : BASE SUR LAQUELLE ON CREE LE VECTEUR
! OUT : NBVAL  : NOMBRE DE VALEURS EXTRAITES
! OUT : TYPVAL : TYPE JEVEUX DES VALEURS EXTRAITES
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, ipar
    integer(kind=8) :: i, j, iv, jvale, jvall, kvale
    character(len=1) :: base
    character(len=4) :: type
    character(len=19) :: nomtab
    character(len=24) :: nomjv, nomjvl, inpar, jnpar
    character(len=24) :: valk
    character(len=24), pointer :: tblp(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
! DEB------------------------------------------------------------------
!
    call jemarq()
!
    nomtab = nomta
    base = basobj(1:1)
    inpar = para
!
!     --- VERIFICATION DE LA BASE ---
!
    ASSERT(base .eq. 'V' .or. base .eq. 'G')
!
!     --- VERIFICATION DE LA TABLE ---
!
    call jeexin(nomtab//'.TBBA', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
!
    call jeveuo(nomtab//'.TBNP', 'L', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    if (nbpara .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nblign .eq. 0) then
        call utmess('F', 'UTILITAI4_76')
    end if
!
!     --- VERIFICATION QUE LE PARAMETRE EXISTE DANS LA TABLE ---
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
    do ipar = 1, nbpara
        jnpar = tblp(1+4*(ipar-1))
        if (inpar .eq. jnpar) goto 12
    end do
    valk = inpar
    call utmess('F', 'UTILITAI6_89', sk=valk)
12  continue
!
    type = tblp(1+4*(ipar-1)+1)
    nomjv = tblp(1+4*(ipar-1)+2)
    nomjvl = tblp(1+4*(ipar-1)+3)
!
    call jeveuo(nomjv, 'L', jvale)
    call jeveuo(nomjvl, 'L', jvall)
    nbval = 0
    do i = 1, nblign
        if (zi(jvall+i-1) .eq. 1) nbval = nbval+1
    end do
!
    iv = 0
    if (type(1:1) .eq. 'I') then
        call wkvect(nomobj, base//' V I', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zi(kvale+j-1) .eq. zi(jvale+i-1)) goto 100
                end do
                iv = iv+1
                zi(kvale+iv-1) = zi(jvale+i-1)
            end if
100         continue
        end do
!
    else if (type(1:1) .eq. 'R') then
        call wkvect(nomobj, base//' V R', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zr(kvale+j-1) .eq. zr(jvale+i-1)) goto 200
                end do
                iv = iv+1
                zr(kvale+iv-1) = zr(jvale+i-1)
            end if
200         continue
        end do
!
    else if (type(1:1) .eq. 'C') then
        call wkvect(nomobj, base//' V C', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zc(kvale+j-1) .eq. zc(jvale+i-1)) goto 300
                end do
                iv = iv+1
                zc(kvale+iv-1) = zc(jvale+i-1)
            end if
300         continue
        end do
!
    else if (type(1:3) .eq. 'K80') then
        call wkvect(nomobj, base//' V K80', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zk80(kvale+j-1) .eq. zk80(jvale+i-1)) goto 400
                end do
                iv = iv+1
                zk80(kvale+iv-1) = zk80(jvale+i-1)
            end if
400         continue
        end do
!
    else if (type(1:3) .eq. 'K32') then
        call wkvect(nomobj, base//' V K32', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zk32(kvale+j-1) .eq. zk32(jvale+i-1)) goto 500
                end do
                iv = iv+1
                zk32(kvale+iv-1) = zk32(jvale+i-1)
            end if
500         continue
        end do
!
    else if (type(1:3) .eq. 'K24') then
        call wkvect(nomobj, base//' V K24', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zk24(kvale+j-1) .eq. zk24(jvale+i-1)) goto 600
                end do
                iv = iv+1
                zk24(kvale+iv-1) = zk24(jvale+i-1)
            end if
600         continue
        end do
!
    else if (type(1:3) .eq. 'K16') then
        call wkvect(nomobj, base//' V K16', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zk16(kvale+j-1) .eq. zk16(jvale+i-1)) goto 700
                end do
                iv = iv+1
                zk16(kvale+iv-1) = zk16(jvale+i-1)
            end if
700         continue
        end do
!
    else if (type(1:2) .eq. 'K8') then
        call wkvect(nomobj, base//' V K8', nbval, kvale)
        do i = 1, nblign
            if (zi(jvall+i-1) .eq. 1) then
                do j = 1, iv
                    if (zk8(kvale+j-1) .eq. zk8(jvale+i-1)) goto 800
                end do
                iv = iv+1
                zk8(kvale+iv-1) = zk8(jvale+i-1)
            end if
800         continue
        end do
    end if
!
    typval = type
    nbval = iv
    call jeecra(nomobj, 'LONUTI', nbval)
!
    call jedema()
end subroutine
