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

subroutine jerazo(nomlu, ni, i1)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "jeveux.h"
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjalty.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjvern.h"
#include "asterfort/jxlocs.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: nomlu
    integer(kind=8), intent(in) :: ni, i1
! ----------------------------------------------------------------------
!     REMISE A "ZERO" DU SEGMENT DE VALEURS ASSOCIE A UN OBJET JEVEUX
! IN  NI    : NOMBRE DE VALEURS A REINITIALISER
! IN  I1    : INDICE DE LA PREMIERE VALEUR
! IN  NOMLU : NOM DE L'OBJET JEVEUX
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibacol, iblono, icre, inat, inatb, iret
    integer(kind=8) :: ixdeso, ixiadd, ixlono, j1, j2, jcara, jctab
    integer(kind=8) :: jdate, jdocu, jgenr, jhcod, jiadd, jiadm, jini
    integer(kind=8) :: jlong, jlono, jltyp, jluti, jmarq, jorig, jrnom
    integer(kind=8) :: jtype, lonoi, ltypi, n
!-----------------------------------------------------------------------
    parameter(n=5)
! ----------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
! ----------------------------------------------------------------------
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: numatr
    common/idatje/numatr
! -------------------------------------------------
    character(len=32) :: noml32
    character(len=8) :: noml8
    character(len=1) :: typei, genri
! ----------------------------------------------------------------------
    integer(kind=8) :: iddeso, idiadd, idlono
    parameter(iddeso=1, idiadd=2,&
     &               idlono=8)
! ----------------------------------------------------------------------
    noml32 = nomlu
    noml8 = noml32(25:32)

!
    icre = 0
    call jjvern(noml32, icre, iret)
    inat = iret
    inatb = iret
    if (iret .eq. 0) then
        call utmess('F', 'JEVEUX_26', sk=noml32(1:24))
        goto 100
    else if (iret .eq. 1) then
        genri = genr(jgenr(iclaos)+idatos)
        typei = type(jtype(iclaos)+idatos)
        ltypi = ltyp(jltyp(iclaos)+idatos)
        if (genri .eq. 'N') then
            call utmess('F', 'JEVEUX1_20', sk=noml32)
        end if
        goto 100
    else if (iret .eq. 2) then
        call jjallc(iclaco, idatco, 'E', ibacol)
        ixiadd = iszon(jiszon+ibacol+idiadd)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        if (noml8 .eq. '$$XATR  ') then
            ixlono = numatr
            iblono = iadm(jiadm(iclaco)+2*ixlono-1)
            genri = genr(jgenr(iclaco)+ixlono)
            ltypi = ltyp(jltyp(iclaco)+ixlono)
            lonoi = lono(jlono(iclaco)+ixlono)*ltypi
            call jxlocs(zi, genri, ltypi, lonoi, iblono, &
                        .false._1, jctab)
            goto 1000
        else
            if (noml8 .ne. '        ') then
                inat = 3
                call jjcroc(noml8, icre)
!            ------ CAS D'UN OBJET DE COLLECTION  ------
                if (ixiadd .ne. 0) inatb = 3
            else
                if (ixiadd .ne. 0) then
!            ----------- COLLECTION DISPERSEE
                    call utmess('F', 'JEVEUX1_21', sk=noml32)
                end if
            end if
            genri = genr(jgenr(iclaco)+ixdeso)
            typei = type(jtype(iclaco)+ixdeso)
            ltypi = ltyp(jltyp(iclaco)+ixdeso)
        end if
    else
        ASSERT(.false.)
    end if
100 continue
    call jjalty(typei, ltypi, 'E', inatb, jctab)
    if (inat .eq. 3 .and. ixiadd .eq. 0) then
        ixlono = iszon(jiszon+ibacol+idlono)
        if (ixlono .gt. 0) then
            iblono = iadm(jiadm(iclaco)+2*ixlono-1)
            lonoi = iszon(jiszon+iblono-1+idatoc+1)-iszon(jiszon+iblono-1+idatoc)
            if (lonoi .gt. 0) then
                jctab = jctab+(iszon(jiszon+iblono-1+idatoc)-1)
            else
                call utmess('F', 'JEVEUX1_22', sk=noml32)
            end if
        else
            jctab = jctab+long(jlong(iclaco)+ixdeso)*(idatoc-1)
        end if
    end if
1000 continue
!
    jini = jctab+i1-1
    j1 = 0
    j2 = ni-1
    if (typei .eq. 'I') then
        do i = j1, j2
            zi(jini+i) = 0
        end do
    else if (typei .eq. 'S') then
        do i = j1, j2
            zi4(jini+i) = 0
        end do
    else if (typei .eq. 'R') then
        do i = j1, j2
            zr(jini+i) = 0.d0
        end do
    else if (typei .eq. 'C') then
        do i = j1, j2
            zc(jini+i) = (0.d0, 0.d0)
        end do
    else if (typei .eq. 'L') then
        do i = j1, j2
            zl(jini+i) = .false.
        end do
    else if (typei .eq. 'K') then
        if (ltypi .eq. 8) then
            do i = j1, j2
                zk8(jini+i) = ' '
            end do
        else if (ltypi .eq. 16) then
            do i = j1, j2
                zk16(jini+i) = ' '
            end do
        else if (ltypi .eq. 24) then
            do i = j1, j2
                zk24(jini+i) = ' '
            end do
        else if (ltypi .eq. 32) then
            do i = j1, j2
                zk32(jini+i) = ' '
            end do
        else if (ltypi .eq. 80) then
            do i = j1, j2
                zk80(jini+i) = ' '
            end do
        end if
    end if
!
end subroutine
