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

subroutine jedetc(clas, souch, ipos)
    implicit none
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjcren.h"
#include "asterfort/jjlidy.h"
#include "asterfort/jjmzat.h"
#include "asterfort/jxlibd.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: clas, souch
    integer(kind=8), intent(in) :: ipos
! ----------------------------------------------------------------------
! DESTRUCTION D'UN ENSEMBLE D'OBJETS JEVEUX
!
! IN  CLAS   : CLASSE DES OBJETS ( ' ' TOUTES LES BASES OUVERTES)
! IN  SOUCH  : SOUS-CHAINE RECHERCHEE
! IN  IPOS   : POSITION DE LA SOUS-CHAINE RECHERCHEE
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmar, iadmi, iadmoc, iadyn, iadyoc, ibacol, ibiadd
    integer(kind=8) :: ibiadm, iblono, ibmarq, ixdeso, ixiadd, ixiadm, ixlono
    integer(kind=8) :: ixmarq, jcara, jdate, jdocu, jgenr, jhcod, jiadd
    integer(kind=8) :: jiadm, jlong, jlono, jltyp, jluti, jmarq, jorig
    integer(kind=8) :: jrnom, jtype, k, kk, l, lonoi, n
!
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    integer(kind=8) :: nivo
    common/jvnivo/nivo
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
!     ------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idmarq, idlono, idnum
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &               idmarq=4,&
     &               idlono=8, idnum=10)
!     ------------------------------------------------------------------
    integer(kind=8) :: ncla1, ncla2, ic, j, iret, id(idnum), nmax, iaddi(2)
    character(len=32) :: crnom, nom32
    character(len=1) :: kclas
! DEB ------------------------------------------------------------------
    l = len(souch)
    ASSERT(ipos+l .le. 25 .and. ipos .ge. 0 .and. l .ne. 0)
    kclas = clas(1:min(1, len(clas)))
    if (kclas .eq. ' ') then
        ncla1 = 1
        ncla2 = index(classe, '$')-1
        if (ncla2 .lt. 0) ncla2 = n
    else
        ncla1 = index(classe, kclas)
        ncla2 = ncla1
    end if
    do ic = ncla1, ncla2
        do kk = 1, 2
            do j = 1, nremax(ic)
                crnom = rnom(jrnom(ic)+j)
                if (crnom(1:1) .eq. '?' .or. crnom(25:32) .ne. '        ') then
                    cycle
                end if
                if (souch .eq. crnom(ipos:ipos+l-1)) then
                    call jjcren(crnom(1:24), 0, iret)
                    if (iret .eq. 1 .and. kk .eq. 2) then
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
                            call utmess('I', 'JEVEUX_7', sk=crnom(1:24))
                        end if
                        call jjcren(crnom(1:24), -1, iret)
                        nomos = '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
                        call jjmzat(iclaos, idatos)
                    else if (iret .gt. 1) then
                        call jjallc(ic, idatco, 'E', ibacol)
                        ixiadm = iszon(jiszon+ibacol+idiadm)
                        ixiadd = iszon(jiszon+ibacol+idiadd)
                        ixlono = iszon(jiszon+ibacol+idlono)
                        ixdeso = iszon(jiszon+ibacol+iddeso)
                        ixmarq = iszon(jiszon+ibacol+idmarq)
                        if (ixiadm .gt. 0) then
                            ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
                            ibiadd = iadm(jiadm(ic)+2*ixiadd-1)
                            ibmarq = iadm(jiadm(ic)+2*ixmarq-1)
                            nmax = iszon(jiszon+ibacol+ivnmax)
                            do k = 1, nmax
                                iadmar = iszon(jiszon+ibmarq-1+2*k)
                                if (iadmar .ne. 0) then
                                    iszon(jiszon+kdesma(1)+iadmar-1) = &
                                        0
                                end if
                                iadmoc = iszon(jiszon+ibiadm-1+2*k-1)
                                iadyoc = iszon(jiszon+ibiadm-1+2*k)
                                call jjlidy(iadyoc, iadmoc)
                                iaddi(1) = iszon(jiszon+ibiadd-1+2*k-1)
                                iaddi(2) = iszon(jiszon+ibiadd-1+2*k)
                                if (iaddi(1) .gt. 0) then
                                    if (ixlono .gt. 0) then
                                        iblono = iadm(jiadm(ic)+2* &
                                                      ixlono-1)
                                        lonoi = iszon(jiszon+iblono+k- &
                                                      1)*ltyp(jltyp(ic)+ixdeso)
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
                                if (nom32(1:24) .eq. crnom(1:24) .or. nom32(25:26) .eq. &
                                    '&&') then
                                    iadmi = iadm(jiadm(ic)+2*id(k)-1)
                                    iadyn = iadm(jiadm(ic)+2*id(k))
                                    call jjlidy(iadyn, iadmi)
                                    iaddi(1) = iadd(jiadd(ic)+2*id(k)-1)
                                    iaddi(2) = iadd(jiadd(ic)+2*id(k))
                                    if (iaddi(1) .gt. 0) then
                                        lonoi = lono(jlono(ic)+id(k)) &
                                                *ltyp(jltyp(ic)+id(k))
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
                                    call utmess('I', 'JEVEUX_7', sk=nom32)
                                end if
                                call jjcren(nom32, -2, iret)
                                call jjmzat(ic, id(k))
                            end if
                        end do
                        crnom = rnom(jrnom(ic)+idatco)
                        iadyn = iadm(jiadm(ic)+2*idatco)
                        call jjlidy(iadyn, ibacol)
                        iaddi(1) = iadd(jiadd(ic)+2*idatco-1)
                        iaddi(2) = iadd(jiadd(ic)+2*idatco)
                        if (iaddi(1) .gt. 0) then
                            lonoi = lono(jlono(ic)+idatco)*ltyp(jltyp(ic)+idatco)
                            call jxlibd(0, idatco, ic, iaddi, lonoi)
                        end if
                        if (nivo .ge. 2) then
                            call utmess('I', 'JEVEUX_7', sk=crnom(1:24))
                        end if
                        call jjcren(crnom(1:24), -2, iret)
                        call jjmzat(ic, idatco)
                        nomco = '$$$$$$$$$$$$$$$$$$$$$$$$'
                    end if
                end if
            end do
        end do
    end do
! FIN ------------------------------------------------------------------
end subroutine
