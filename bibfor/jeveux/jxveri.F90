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
subroutine jxveri()
    implicit none
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjvern.h"
! ----------------------------------------------------------------------
! VERIFIE L'INTEGRITE DU CHAINAGE AVANT DES SEGMENTS DE VALEURS ET DE LA
! ZONE MEMOIRE UTILISEE
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
    integer(kind=8) :: isstat
    common/iconje/isstat
    integer(kind=8) :: istat
    common/istaje/istat(4)
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmi, iadmoc, iadyn, iadyoc, ibacol, ibiadm, ic
    integer(kind=8) :: idm, il, iret, isd, isdc, isf
    integer(kind=8) :: ixiadm, j, jcara, jdate, jdocu
    integer(kind=8) :: jgenr, jhcod, jiadd, jiadm, jlong, jlono, jltyp
    integer(kind=8) :: jluti, jmarq, jorig, jrnom, jtype, k, n
    integer(kind=8) :: ncla1, ncla2, nmax
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
    integer(kind=8) :: ivnmax, idiadm
    parameter(ivnmax=0, idiadm=3)
! ----------------------------------------------------------------------
    character(len=32) :: nom32
    character(len=1) :: cgenr
! DEB ------------------------------------------------------------------
!
    nom32 = '??'
!
!     ON TRAITE LES OBJETS ALLOUES EN MEMOIRE DYNAMIQUE
!
    if (ldyn .ne. 1 .and. ldyn .ne. 2) goto 300
    ncla1 = 1
    ncla2 = index(classe, '$')-1
    if (ncla2 .lt. 0) ncla2 = n
    do ic = ncla2, ncla1, -1
        do j = 1, nremax(ic)
            iadmi = iadm(jiadm(ic)+2*j-1)
            iadyn = iadm(jiadm(ic)+2*j)
            if (iadmi .eq. 0 .or. iadyn .eq. 0) goto 205
            cgenr = genr(jgenr(ic)+j)
            nom32 = rnom(jrnom(ic)+j)
!
            isdc = iszon(jiszon+iadmi-1)/isstat
            ASSERT(isdc .eq. 1 .or. isdc .eq. 2)
            if (cgenr .eq. 'X' .and. isdc .eq. 2) then
                call jjvern(nom32, 0, iret)
                call jjallc(ic, j, 'L', ibacol)
                ixiadm = iszon(jiszon+ibacol+idiadm)
                nmax = iszon(jiszon+ibacol+ivnmax)
                if (ixiadm .gt. 0) then
                    ibiadm = iadm(jiadm(ic)+2*ixiadm-1)
                    do k = 1, nmax
                        iadmoc = iszon(jiszon+ibiadm-1+2*k-1)
                        iadyoc = iszon(jiszon+ibiadm-1+2*k)
                        if (iadyoc .ne. 0) then
                            idm = iadmoc-4
                            isd = iszon(jiszon+idm+3)/isstat
                            ASSERT(isd .eq. 1 .or. isd .eq. 2)
                            isf = iszon(jiszon+iszon(jiszon+idm)-4)/isstat
                            ASSERT(isf .eq. 3 .or. isf .eq. 4)
                            il = iszon(jiszon+idm)-8-idm
                            ASSERT(il .gt. 0)
                        end if
                    end do
                end if
                call jjlide('JEIMPO', nom32(1:24), 2)
                goto 205
            else
                idm = iadmi-4
                isd = iszon(jiszon+idm+3)/isstat
                ASSERT(isd .eq. 1 .or. isd .eq. 2)
                isf = iszon(jiszon+iszon(jiszon+idm)-4)/isstat
                ASSERT(isf .eq. 3 .or. isf .eq. 4)
                il = iszon(jiszon+idm)-8-idm
                ASSERT(il .gt. 0)
            end if
205         continue
        end do
    end do
!
300 continue
! FIN ------------------------------------------------------------------
end subroutine
