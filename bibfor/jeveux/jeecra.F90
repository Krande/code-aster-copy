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

subroutine jeecra(nomlu, catr, ival, cval)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjalls.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjecrs.h"
#include "asterfort/jjprem.h"
#include "asterfort/jjvern.h"
#include "asterfort/utmess.h"
#include "asterfort/assert.h"
    character(len=*), intent(in) :: nomlu
    character(len=*), intent(in) :: catr
    integer(kind=8), intent(in), optional :: ival
    character(len=*), intent(in), optional :: cval
! ----------------------------------------------------------------------
! ROUTINE UTILISATEUR D'AFFECTATION D'UN ATTRIBUT
!
! IN  NOMLU  : NOM DE L'OBJET JEVEUX
! IN  CATR   : NOM DE L'ATTRIBUT
! IN  IVAL   : VALEUR EN ENTIER DE L'ATTRIBUT
! IN  CVAL   : VALEUR EN CARACTERE DE L'ATTRIBUT
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!     ------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmi, iadyn, iblong, iblono, ibluti, ic, id, ibnum
    integer(kind=8) :: il0, il1, ixlono, ixluti, jcara, jdate, jdocu, ixnum
    integer(kind=8) :: jgenr, jhcod, jiadd, jiadm, jitab, jlong
    integer(kind=8) :: jlono, jltyp, jluti, jmarq, jorig, jrnom, jtype
    integer(kind=8) :: longi, longj, lonoi, lonoj, lonok, lont, lonti
    integer(kind=8) :: ltypi, n, nbl, nhc, nmaxi, nbv
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
!     ------------------------------------------------------------------
    character(len=32) :: noml32
    character(len=1) :: genri, typei
    character(len=8) :: catrlu
    character(len=32), dimension(2) :: lvalk
    aster_logical :: lconst, lconti, llong, lluti
    integer(kind=8) :: icre, iret, itab(1), jtab, irt
    integer(kind=8) :: ibacol, ixiadd, ixdeso, ixlong, ii
!     ------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idlong, idlono, idluti, idnum
    parameter(ivnmax=0, iddeso=1, idiadd=2,&
     &               idlong=7,&
     &               idlono=8, idluti=9, idnum=10)
    integer(kind=8) :: ilorep, ideno, ilnom, ilmax, iluti, idehc
    parameter(ilorep=1, ideno=2, ilnom=3, ilmax=4, iluti=5, idehc=6)
! DEB ------------------------------------------------------------------
!
    catrlu = catr
    noml32 = nomlu
    irt = 0
!
! --- CAS GENERAL
!
    icre = 0
    call jjvern(noml32, icre, iret)
!
    if (iret .eq. 0) then
        call utmess('F', 'JEVEUX_26', sk=noml32(1:24))
    else if (iret .eq. 1) then
        ic = iclaos
        id = idatos
        lconst = .true.
        lconti = .false.
        ixlong = id
        ixlono = id
        ixluti = id
    else
        ic = iclaco
        id = idatco
        call jjallc(ic, id, 'E', ibacol)
        if (noml32(25:32) .ne. '        ') then
            iret = 3
            call jjcroc(noml32(25:32), icre)
        end if
        ixdeso = iszon(jiszon+ibacol+iddeso)
        id = ixdeso
        ixiadd = iszon(jiszon+ibacol+idiadd)
        lconti = (ixiadd .eq. 0)
        ixlong = iszon(jiszon+ibacol+idlong)
        ixlono = iszon(jiszon+ibacol+idlono)
        ixluti = iszon(jiszon+ibacol+idluti)
        ixnum = iszon(jiszon+ibacol+idnum)
        lconst = (ixlong .eq. 0)
        nmaxi = iszon(jiszon+ibacol+ivnmax)
    end if
!
    genri = genr(jgenr(ic)+id)
    typei = type(jtype(ic)+id)
    if (catrlu .eq. 'LONT    ') then
        if (.not. lconti) then
            call utmess('F', 'JEVEUX_98', sk=catrlu)
        else
            llong = .false.
            lluti = .false.
        end if
    else if (catrlu .eq. 'LONMAX' .or. catrlu .eq. 'NOMMAX' .or. catrlu .eq. 'LONUTI' .or. &
             catrlu .eq. 'DOCU' .or. catrlu .eq. 'DATE' .or. catrlu .eq. 'NUTIOC') then
        llong = (catrlu(4:6) .eq. 'MAX')
        lluti = (catrlu(4:6) .eq. 'UTI')
        if ((genri .ne. 'N' .and. catrlu(1:3) .eq. 'NOM') .or. &
            (genri .ne. 'V' .and. catrlu(1:4) .eq. 'LONM') .or. &
            (genri .ne. 'V' .and. catrlu(1:4) .eq. 'LONU')) then
            call utmess('F', 'JEVEUX_99', sk=genri)
        end if
    else
        call utmess('F', 'JEVEUX1_23', sk=catrlu)
    end if
!
    if (catrlu .eq. 'LONT    ' .and. lconti) then
        lono(jlono(ic)+id) = ival
        if (lconst) long(jlong(ic)+id) = ival/nmaxi
    else if (catrlu .eq. 'DATE    ') then
        date(jdate(ic)+id) = ival
    else if (catrlu .eq. 'DOCU    ') then
        docu(jdocu(ic)+id) = cval
    else if (catrlu .eq. 'NUTIOC') then
        if (iret .eq. 2 .and. lconti) then
            ibnum = iadm(jiadm(ic)+2*ixnum-1)
            nbv = iszon(jiszon+ibnum)
            if (ival > nbv) then
                print *, "ival = ", ival, "nbv = ", nbv
            end if
            ASSERT(ival .le. nbv)
            iszon(jiszon+ibnum-1+2) = ival
            if (.not. lconst) then
                iblong = iadm(jiadm(ic)+2*ixlong-1)
                iblono = iadm(jiadm(ic)+2*ixlono-1)
                do ii = 1, ival
                    iszon(jiszon+iblong-1+ii) = ( &
                                                iszon(jiszon+iblono+ii)-iszon(jiszon+iblono-1+ii) &
                                                )
                end do
            end if
        else
            ASSERT(.false.)
        end if
    else if (lconst) then
        if (llong) then
            lonoi = lono(jlono(ic)+id)
            if (catrlu .eq. 'LONMAX  ' .or. catrlu .eq. 'NOMMAX  ') then
                longi = long(jlong(ic)+id)
                longj = ival
                if (ival .le. 0) then
                    lvalk(1) = catrlu
                    lvalk(2) = noml32(1:24)
                    call utmess('F', 'JEVEUX1_67', nk=2, valk=lvalk, si=ival)
                end if
            end if
            if (longi .ne. 0) then
                call utmess('F', 'JEVEUX1_01', sk=catrlu)
            else
                long(jlong(ic)+id) = longj
                if (lonoi .ne. 0 .and. iret .eq. 1) then
                    call utmess('F', 'JEVEUX1_02', sk=catrlu)
                else
                    if (genri .eq. 'V') then
                        lono(jlono(ic)+id) = longj
                    else if (genri .eq. 'N') then
                        ltypi = ltyp(jltyp(ic)+id)
                        lonok = (idehc+jjprem(longj, irt))*lois+(longj+1)*ltypi
                        if (mod(lonok, ltypi) .gt. 0) then
                            lonok = (lonok/ltypi+1)
                        else
                            lonok = lonok/ltypi
                        end if
                        lono(jlono(ic)+id) = lonok
                        luti(jluti(ic)+id) = 0
                        if (iadm(jiadm(ic)+2*id-1) .eq. 0) then
                            nbl = lonok*ltypi
                            call jjalls(nbl, ic, genri, typei, ltypi, &
                                        'INIT', itab, jtab, iadmi, iadyn)
                            iadm(jiadm(ic)+2*id-1) = iadmi
                            iadm(jiadm(ic)+2*id) = iadyn
                            call jjecrs(iadmi, ic, id, 0, 'E', &
                                        imarq(jmarq(ic)+2*id-1))
                            nhc = jjprem(ival, irt)
                            jitab = jiszon+iadmi-1
                            iszon(jitab+ilorep) = nhc
                            iszon(jitab+ideno) = (idehc+nhc)*lois
                            iszon(jitab+ilnom) = ltypi
                            iszon(jitab+ilmax) = ival
                            iszon(jitab+iluti) = 0
                            iszon(jitab+idehc) = idehc
                        end if
                    end if
                    if (lconti) then
                        if (lonoi .ne. 0 .and. lonoi .lt. nmaxi*lono(jlono(ic)+id)) then
                            call utmess('F', 'JEVEUX1_03', sk=catrlu)
                        else
                            lono(jlono(ic)+id) = nmaxi*lono( &
                                                 jlono(ic)+id)
                        end if
                    end if
                end if
            end if
        else if (lluti) then
            if (catrlu .eq. 'LONUTI  ') then
                luti(jluti(ic)+id) = ival
            end if
        else
            call utmess('F', 'JEVEUX1_04', sk=catrlu)
        end if
    else if (iret .eq. 3) then
        if (llong .and. .not. lconst) then
            iblong = iadm(jiadm(ic)+2*ixlong-1)
            iblono = iadm(jiadm(ic)+2*ixlono-1)
            if (lconti) then
                if (idatoc .eq. 1) then
                    if (iszon(jiszon+iblono-1+idatoc) .eq. 0) then
                        iszon(jiszon+iblono-1+idatoc) = 1
                    end if
                end if
                il1 = jiszon+iblono-1+idatoc+1
                il0 = jiszon+iblono-1+idatoc
                if (iszon(il0) .eq. 0) then
                    call utmess('F', 'JEVEUX1_05', sk=catrlu)
                else
                    lonti = iszon(il0)
                    lonoi = 0
                    if (iszon(il1) .ne. 0) then
                        lonti = max(iszon(il1), iszon(il0))
                        lonoi = max(iszon(il1)-iszon(il0), 0)
                    end if
                end if
            else
                lonoi = iszon(jiszon+iblono-1+idatoc)
            end if
            if (lonoi .ne. 0) then
                call utmess('F', 'JEVEUX1_06', sk=catrlu)
            end if
            if (catrlu .eq. 'LONMAX  ') then
                longi = iszon(jiszon+iblong-1+idatoc)
                longj = ival
                lonoj = longj
            end if
            if (longi .ne. 0) then
                call utmess('F', 'JEVEUX1_01', sk=catrlu)
            else
                if (lconti) then
                    lont = lono(jlono(ic)+id)
                    if (lont .ne. 0 .and. lonti-1+lonoj .gt. lont) then
                        call utmess('F', 'JEVEUX1_07', sk=catrlu)
                    else
                        iszon(jiszon+iblono-1+idatoc+1) = lonti+lonoj
                    end if
                else
                    iszon(jiszon+iblono-1+idatoc) = lonoj
                end if
                iszon(jiszon+iblong-1+idatoc) = longj
                luti(jluti(ic)+ixlono) = 1+luti(jluti(ic)+ &
                                                ixlono)
            end if
        else if (lluti) then
            ibluti = iadm(jiadm(ic)+2*ixluti-1)
            if (catrlu .eq. 'LONUTI  ') then
                iszon(jiszon+ibluti-1+idatoc) = ival
            end if
        end if
    else
        call utmess('F', 'JEVEUX1_04', sk=catrlu)
    end if
! FIN ------------------------------------------------------------------
end subroutine
