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
subroutine jedetr(nomlu)
    implicit none
#include "asterf_config.h"
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjcren.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jjmzat.h"
#include "asterfort/jjvern.h"
#include "asterfort/jxlibd.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: nomlu
! ----------------------------------------------------------------------
! DESTRUCTION D'UN OBJET JEVEUX
!
! IN  NOMLU  : NOM D'OBJET JEVEUX
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmar, iadmi, iadmoc, iadyn, iadyoc, ibacol, ibiadd
    integer(kind=8) :: ibiadm, iblong, iblono, ibmarq, ic, ixdeso, ixiadd
    integer(kind=8) :: ixiadm, ixlong, ixlono, ixmarq, ixnom, jcara, jdate
    integer(kind=8) :: jdocu, jgenr, jhcod, jiadd, jiadm, jlong, jlono
    integer(kind=8) :: jltyp, jluti, jmarq, jorig, jrnom, jtype, k
    integer(kind=8) :: lonoi, n, nmax
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: nivo
    common/jvnivo/nivo
!     ------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idmarq, idnom, idlong, idlono
    integer(kind=8) :: idnum
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &               idmarq=4, idnom=5, idlong=7,&
     &               idlono=8, idnum=10)
!     ------------------------------------------------------------------
    character(len=32) :: noml32, nom32
    integer(kind=8) :: icre, iret, id(idnum), iaddi(2)
! DEB ------------------------------------------------------------------
    noml32 = nomlu
    icre = 0
! #ifdef ASTER_DEBUG_CXX_OBJECTS
!     nivo = 2
! #endif
    call jjvern(noml32, icre, iret)
!
    if (iret .eq. 0) then
        goto 999
    else if (iret .eq. 1) then
        ic = iclaos
        iadmi = iadm(jiadm(ic)+2*idatos-1)
        iadyn = iadm(jiadm(ic)+2*idatos)
        call jjlidy(iadyn, iadmi)
        iaddi(1) = iadd(jiadd(ic)+2*idatos-1)
        iaddi(2) = iadd(jiadd(ic)+2*idatos)
        if (iaddi(1) .gt. 0) then
            lonoi = lono(jlono(ic)+idatos)*ltyp(jltyp(ic)+idatos)
            call jxlibd(0, idatos, ic, iaddi, lonoi)
        end if
        if (nivo .ge. 2) then
            call utmess('I', 'JEVEUX_07', sk=noml32)
        end if
        call jjcren(noml32, -1, iret)
        nomos = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        call jjmzat(iclaos, idatos)
    else
        ic = iclaco
        call jjallc(ic, idatco, 'E', ibacol)
        if (noml32(25:32) .ne. '        ') then
            call jjcroc(noml32(25:32), icre)
            iret = 3
        end if
        if (iret .eq. 2) then
            ixiadm = iszon(jiszon+ibacol+idiadm)
            ixiadd = iszon(jiszon+ibacol+idiadd)
            ixlono = iszon(jiszon+ibacol+idlono)
            ixdeso = iszon(jiszon+ibacol+iddeso)
            ixmarq = iszon(jiszon+ibacol+idmarq)
            if (ixiadm .gt. 0) then
                ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
                ibmarq = iadm(jiadm(ic)+2*ixmarq-1)
                nmax = iszon(jiszon+ibacol+ivnmax)
                do k = 1, nmax
                    iadmar = iszon(jiszon+ibmarq-1+2*k)
                    if (iadmar .ne. 0) then
                        iszon(jiszon+kdesma(1)+iadmar-1) = 0
                    end if
                    iadmoc = iszon(jiszon+ibiadm-1+2*k-1)
                    iadyoc = iszon(jiszon+ibiadm-1+2*k)
                    call jjlidy(iadyoc, iadmoc)
                    ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
                    iaddi(1) = iszon(jiszon+ibiadd-1+2*k-1)
                    iaddi(2) = iszon(jiszon+ibiadd-1+2*k)
                    if (iaddi(1) .gt. 0) then
                        if (ixlono .gt. 0) then
                            iblono = iadm(jiadm(ic)+2*ixlono-1)
                            lonoi = iszon(jiszon+iblono+k-1)*ltyp( &
                                    jltyp(ic)+ixdeso)
                        else
                            lonoi = lono(jlono(ic)+ixdeso)*ltyp(jltyp(ic)+ixdeso)
                        end if
                        call jxlibd(idatco, k, ic, iaddi, lonoi)
                    end if
                end do
            end if
            do k = 1, idnum
                id(k) = iszon(jiszon+ibacol+k)
                if (id(k) .gt. 0) then
                    nom32 = rnom(jrnom(ic)+id(k))
                    if (nom32(1:24) .eq. noml32(1:24) .or. nom32(25:26) .eq. '&&') then
                        iadmi = iadm(jiadm(ic)+2*id(k)-1)
                        iadyn = iadm(jiadm(ic)+2*id(k))
                        call jjlidy(iadyn, iadmi)
                        iaddi(1) = iadd(jiadd(ic)+2*id(k)-1)
                        iaddi(2) = iadd(jiadd(ic)+2*id(k))
                        if (iaddi(1) .gt. 0) then
                            lonoi = lono(jlono(ic)+id(k))*ltyp(jltyp(ic)+id(k))
                            call jxlibd(0, id(k), ic, iaddi, lonoi)
                        end if
                    else
                        id(k) = 0
                    end if
                end if
            end do
            do k = 1, idnum
                if (id(k) .gt. 0) then
                    nom32 = rnom(jrnom(ic)+id(k))
                    if (nivo .ge. 2) then
                        call utmess('I', 'JEVEUX_07', sk=noml32(1:24))
                    end if
                    call jjcren(nom32, -2, iret)
                    call jjmzat(ic, id(k))
                end if
            end do
            iadyn = iadm(jiadm(ic)+2*idatco)
            call jjlidy(iadyn, ibacol)
            iaddi(1) = iadd(jiadd(ic)+2*idatco-1)
            iaddi(2) = iadd(jiadd(ic)+2*idatco)
            if (iaddi(1) .gt. 0) then
                lonoi = lono(jlono(ic)+idatco)*ltyp(jltyp(ic)+idatco)
                call jxlibd(0, idatco, ic, iaddi, lonoi)
            end if
            if (nivo .ge. 2) then
                call utmess('I', 'JEVEUX_07', sk=noml32(1:24))
            end if
            call jjcren(noml32(1:24), -2, iret)
            call jjmzat(ic, idatco)
            nomco = '$$$$$$$$$$$$$$$$$$$$$$$$'
        else if (iret .eq. 3) then
            ixiadm = iszon(jiszon+ibacol+idiadm)
            ixiadd = iszon(jiszon+ibacol+idiadd)
            ixlong = iszon(jiszon+ibacol+idlong)
            ixnom = iszon(jiszon+ibacol+idnom)
            ixdeso = iszon(jiszon+ibacol+iddeso)
            ixlono = iszon(jiszon+ibacol+idlono)
            ixmarq = iszon(jiszon+ibacol+idmarq)
!
!         DESTRUCTION D''UN OBJET DE COLLECTION CONTIGUE REFUSEE
            ASSERT(ixiadd .gt. 0)
!
!         DESTRUCTION DANS UNE COLLECTION NON NOMMEE REFUSEE
            ASSERT(ixnom .gt. 0)
!
            ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
            iaddi(1) = iszon(jiszon+ibiadd-1+2*idatoc-1)
            iaddi(2) = iszon(jiszon+ibiadd-1+2*idatoc)
            if (iaddi(1) .gt. 0) then
                if (ixlono .gt. 0) then
                    iblono = iadm(jiadm(ic)+2*ixlono-1)
                    lonoi = iszon(jiszon+iblono+idatoc-1)*ltyp(jltyp(ic)+ixdeso)
                else
                    lonoi = lono(jlono(ic)+ixdeso)*ltyp(jltyp(ic)+ixdeso)
                end if
                call jxlibd(idatco, idatoc, ic, iaddi, lonoi)
            end if
            iszon(jiszon+ibiadd+idatoc-1) = 0
            ibmarq = iadm(jiadm(ic)+2*ixmarq-1)
            iadmar = iszon(jiszon+ibmarq-1+2*idatoc)
            if (iadmar .ne. 0) then
                iszon(jiszon+kdesma(1)+iadmar-1) = 0
            end if
            ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
            iadmi = iszon(jiszon+ibiadm-1+2*idatoc-1)
            iadyn = iszon(jiszon+ibiadm-1+2*idatoc)
            call jjlidy(iadyn, iadmi)
            iszon(jiszon+ibiadm-1+2*idatoc-1) = 0
            iszon(jiszon+ibiadm-1+2*idatoc) = 0
            if (ixlong .gt. 0) then
                iblong = iadm(jiadm(ic)+2*ixlong-1)
                iszon(jiszon+iblong+idatoc-1) = 0
            end if
            if (nivo .ge. 2) then
                call utmess('I', 'JEVEUX_07', sk=noml32)
            end if
            call jjcroc(nomlu(25:32), -3)
            nomoc = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
        end if
    end if
999 continue
! FIN ------------------------------------------------------------------
end subroutine
