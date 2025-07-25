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
subroutine jenuno(nomlu, nomo)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjvern.h"
#include "asterfort/jxveuo.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: nomlu
    character(len=*), intent(out) :: nomo
! ----------------------------------------------------------------------
! RENVOIE LE NOM ASSOCIE A UN IDENTIFICATEUR
!
! IN  NOMLU  : NOM DE LA COLLECTION OU DU REPERTOIRE
!              L' APPEL DOIT S'EFFECTUER PAR L'INTERMEDIAIRE DE JEXNUM
! IN  NOMO   : NOM DANS LE REPERTOIRE
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmex, iadmi, ibacol, ideco, idenom, ipgcex, ixnom
    integer(kind=8) :: jcara, jctab, jdate, jdocu, jgenr, jhcod, jiadd
    integer(kind=8) :: jiadm, jlong, jlono, jltyp, jluti, jmarq, jorig
    integer(kind=8) :: jrnom, jtype, k, kadm, lnom, lutii, n
    integer(kind=8) :: nk
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
! ----------------------------------------------------------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: numec
    common/inumje/numec
! ----------------------------------------------------------------------
    integer(kind=8) :: ideno, ilnom
    parameter(ideno=2, ilnom=3)
! ----------------------------------------------------------------------
    integer(kind=8) :: idnom
    parameter(idnom=5)
! ----------------------------------------------------------------------
    character(len=32) :: noml32
    character(len=1) :: genri
    integer(kind=8) :: icre, iret, itab(1), vali(2)
! DEB ------------------------------------------------------------------
    ipgcex = ipgc
    ipgc = -2
!
    icre = 0
    noml32 = nomlu
    call jjvern(noml32, icre, iret)
!
    if (iret .eq. 0) then
        call utmess('F', 'JEVEUX_25', sk=noml32(1:24))
    else
        if (iret .eq. 1) then
!
! ------- OBJET DE TYPE REPERTOIRE
!
            genri = genr(jgenr(iclaos)+idatos)
            if (genri .ne. 'N') then
                call utmess('F', 'JEVEUX1_12', sk=noml32)
            end if
            lutii = luti(jluti(iclaos)+idatos)
            if (lutii .lt. numec .or. numec .le. 0) then
                vali(1) = lutii
                vali(2) = numec
                call utmess('F', 'JEVEUX1_13', sk=noml32, ni=2, vali=vali)
            end if
            iadmi = iadm(jiadm(iclaos)+2*idatos-1)
            iadmex = iadmi
            if (iadmex .eq. 0) then
                call jxveuo('L', itab, iret, jctab)
                iadmi = iadm(jiadm(iclaos)+2*idatos-1)
            end if
            kadm = iadmi
            idenom = iszon(jiszon+kadm-1+ideno)
            lnom = iszon(jiszon+kadm-1+ilnom)
            ideco = (kadm-1)*lois+idenom+lnom*(numec-1)
            nk = min(len(nomo), lnom)
            nomo = ' '
            do k = 1, nk
                nomo(k:k) = k1zon(jk1zon+ideco+k)
            end do
            if (iadmex .eq. 0) then
                call jjlide('JENUNO', noml32, iret)
            end if
        else if (iret .eq. 2) then
!
! ------- REPERTOIRE DE COLLECTION
!
            call jjallc(iclaco, idatco, 'L', ibacol)
            ixnom = iszon(jiszon+ibacol+idnom)
            if (ixnom .eq. 0) then
                call utmess('F', 'JEVEUX1_14', sk=noml32)
            end if
            lutii = luti(jluti(iclaco)+ixnom)
            if (lutii .lt. numec .or. numec .le. 0) then
                vali(1) = lutii
                vali(2) = numec
                call utmess('F', 'JEVEUX1_13', sk=noml32, ni=2, vali=vali)
            end if
            iadmi = iadm(jiadm(iclaco)+2*ixnom-1)
            kadm = iadmi
            idenom = iszon(jiszon+kadm-1+ideno)
            lnom = iszon(jiszon+kadm-1+ilnom)
            ideco = (kadm-1)*lois+idenom+lnom*(numec-1)
            nomo = ' '
            do k = 1, min(len(nomo), lnom)
                nomo(k:k) = k1zon(jk1zon+ideco+k)
            end do
            call jjlide('JENUNO', nomlu(1:24), 2)
        else
            ASSERT(.false.)
        end if
    end if
    ipgc = ipgcex
! FIN ------------------------------------------------------------------
end subroutine
