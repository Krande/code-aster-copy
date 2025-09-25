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
subroutine juveca(nom, long)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: nom
    integer(kind=8) :: long
!     REDIMENSIONNEMENT D'UN OBJET SIMPLE JEVEUX DEJA EXISTANT
!     ------------------------------------------------------------------
! IN  NOM  : K24 : NOM DE L'OBJET A REDIMENSIONNER
! IN  LONG : I   : NOUVELLE LONGUEUR DU VECTEUR
!     ------------------------------------------------------------------
!     REMARQUE: LES VALEURS SONT RECOPIEES
!      SI LA NOUVELLE LONGUEUR EST INFERIEURE A L'ANCIENNE, DES VALEURS
!      SONT PERDUES
!     ------------------------------------------------------------------
!
!
    character(len=8) :: base, type
    character(len=32) :: valk(2)
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ldec, ll, lonma2, lonmax, lonuti
    integer(kind=8) :: ltamp, ltyp
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(nom, 'L', ldec)
!
!     --- TYPE, LONGUEUR ET BASE DE L'OBJET A REDIMENSIONNER
    call jelira(nom, 'TYPE  ', cval=type)
    call jelira(nom, 'LONMAX', lonmax)
    call jelira(nom, 'LONUTI', lonuti)
    call jelira(nom, 'CLAS', cval=base)
!
!     -- LONMA2 : LONGUEUR DE RECOPIE :
    ASSERT(lonmax .gt. 0)
    ASSERT(long .gt. 0)
    lonma2 = min(long, lonmax)
!
!     --- ALLOCATION D'UN TAMPON ---
    if (type(1:1) .ne. 'K') then
        call wkvect('&&JUVECA.TAMPON', 'V V '//type, lonma2, ltamp)
    else
        call jelira(nom, 'LTYP', ltyp)
        call codent(ltyp, 'G', type(2:))
        call wkvect('&&JUVECA.TAMPON', 'V V '//type, lonma2, ltamp)
    end if
!
!     --- RECOPIE L'OBJET DANS LE TAMPON ---
    if (type .eq. 'I') then
        do i = 1, lonma2
            zi(ltamp+i-1) = zi(ldec+i-1)
        end do
    else if (type .eq. 'S') then
        do i = 1, lonma2
            zi4(ltamp+i-1) = zi4(ldec+i-1)
        end do
    else if (type .eq. 'R') then
        do i = 1, lonma2
            zr(ltamp+i-1) = zr(ldec+i-1)
        end do
    else if (type .eq. 'C') then
        do i = 1, lonma2
            zc(ltamp+i-1) = zc(ldec+i-1)
        end do
    else if (type .eq. 'L') then
        do i = 1, lonma2
            zl(ltamp+i-1) = zl(ldec+i-1)
        end do
    else if (type(1:1) .eq. 'K') then
        if (ltyp .eq. 8) then
            do i = 1, lonma2
                zk8(ltamp+i-1) = zk8(ldec+i-1)
            end do
        else if (ltyp .eq. 16) then
            do i = 1, lonma2
                zk16(ltamp+i-1) = zk16(ldec+i-1)
            end do
        else if (ltyp .eq. 24) then
            do i = 1, lonma2
                zk24(ltamp+i-1) = zk24(ldec+i-1)
            end do
        else if (ltyp .eq. 32) then
            do i = 1, lonma2
                zk32(ltamp+i-1) = zk32(ldec+i-1)
            end do
        else if (ltyp .eq. 80) then
            do i = 1, lonma2
                zk80(ltamp+i-1) = zk80(ldec+i-1)
            end do
        else
            valk(1) = nom
            valk(2) = type
            call utmess('F', 'JEVEUX_31', nk=2, valk=valk)
        end if
    else
        valk(1) = nom
        valk(2) = type
        call utmess('F', 'JEVEUX_31', nk=2, valk=valk)
    end if
!
!     --- DESTRUCTION DU VIEUX ET CREATION DU NEUF ---
    call jedetr(nom)
    call wkvect(nom, base//' V '//type, long, ldec)
!
!     --- RECOPIE DU TAMPON DANS L'OBJET DEFINITIF ---
    if (type .eq. 'I') then
        do i = 1, lonma2
            zi(ldec+i-1) = zi(ltamp+i-1)
        end do
    else if (type .eq. 'S') then
        do i = 1, lonma2
            zi4(ldec+i-1) = zi4(ltamp+i-1)
        end do
    else if (type .eq. 'R') then
        do i = 1, lonma2
            zr(ldec+i-1) = zr(ltamp+i-1)
        end do
    else if (type .eq. 'C') then
        do i = 1, lonma2
            zc(ldec+i-1) = zc(ltamp+i-1)
        end do
    else if (type .eq. 'L') then
        do i = 1, lonma2
            zl(ldec+i-1) = zl(ltamp+i-1)
        end do
        do i = lonma2+1, long
            zl(ldec+i-1) = .false.
        end do
    else if (type(1:1) .eq. 'K') then
        if (ltyp .eq. 8) then
            do i = 1, lonma2
                zk8(ldec+i-1) = zk8(ltamp+i-1)
            end do
        else if (ltyp .eq. 16) then
            do i = 1, lonma2
                zk16(ldec+i-1) = zk16(ltamp+i-1)
            end do
        else if (ltyp .eq. 24) then
            do i = 1, lonma2
                zk24(ldec+i-1) = zk24(ltamp+i-1)
            end do
        else if (ltyp .eq. 32) then
            do i = 1, lonma2
                zk32(ldec+i-1) = zk32(ltamp+i-1)
            end do
        else if (ltyp .eq. 80) then
            do i = 1, lonma2
                zk80(ldec+i-1) = zk80(ltamp+i-1)
            end do
        end if
    end if
    ll = min(lonuti, long)
    if (lonuti .gt. 0) call jeecra(nom, 'LONUTI', ll)
!
!     --- DESTRUCTION DU TAMPON ---
    call jedetr('&&JUVECA.TAMPON')
    call jedema()
end subroutine
