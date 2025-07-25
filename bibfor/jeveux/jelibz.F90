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

subroutine jelibz(clas)
    implicit none
#include "jeveux_private.h"
#include "asterfort/jjlide.h"
    character(len=*) :: clas
! ----------------------------------------------------------------------
! LIBERATION DE L'ENSEMBLE DES OBJETS MARQUES PAR -1
!
! IN  CLAS   : CLASSE DES OBJETS A LIBERER
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: jcara, jdate, jdocu, jgenr, jhcod, jiadd, jiadm
    integer(kind=8) :: jlong, jlono, jltyp, jluti, jmarq, jorig, jrnom
    integer(kind=8) :: jtype, n, nmax
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
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
    integer(kind=8) :: ivnmax, idiadm, idmarq, idnum
    parameter(ivnmax=0, idiadm=3,&
     &               idmarq=4,&
     &                 idnum=10)
! ----------------------------------------------------------------------
    integer(kind=8) :: ncla1, ncla2, ibacol, ibmarq, ic, id, ix
    integer(kind=8) :: j, k, marqi, iclasi
    character(len=32) :: crnom, d32
    character(len=1) :: kclas
    data d32/'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'/
! DEB ------------------------------------------------------------------
    kclas = clas
    iclasi = iclas
    if (kclas .eq. ' ') then
        ncla1 = 1
        ncla2 = index(classe, '$')-1
        if (ncla2 .lt. 0) ncla2 = n
    else
        ncla1 = index(classe, kclas)
        ncla2 = ncla1
    end if
    do ic = ncla1, ncla2
        do j = 1, nremax(ic)
            crnom = rnom(jrnom(ic)+j)
            if (crnom(1:1) .eq. '?' .or. crnom(25:26) .eq. '$$') goto 150
!          CALL JJCREN ( CRNOM , 0 , IRET )
            if (genr(jgenr(ic)+j) .eq. 'X') then
                iclas = ic
                iclaco = ic
                idatco = j
                nomco = crnom
                nomoc = d32
                if (iclasi .ne. iclaco) then
                    nomos = d32
                end if
                ibacol = iadm(jiadm(ic)+2*j-1)
                if (ibacol .eq. 0) goto 150
                id = iszon(jiszon+ibacol+idiadm)
                if (id .gt. 0) then
!
! ------------- COLLECTION DISPERSEE ( OBJETS DE COLLECTION )
!
                    ix = iszon(jiszon+ibacol+idmarq)
                    ibmarq = iadm(jiadm(ic)+2*ix-1)
                    nmax = iszon(jiszon+ibacol+ivnmax)
                    do k = 1, nmax
                        marqi = iszon(jiszon+ibmarq-1+2*k-1)
                        if (marqi .eq. -1) then
                            call jjlide('JELIBZ', crnom, 2)
                            goto 171
                        end if
                    end do
                end if
!
! ---------- COLLECTION CONTIGUE OU DISPERSEE ( OBJETS ATTRIBUTS )
!
                do k = idnum, 1, -1
                    id = iszon(jiszon+ibacol+k)
                    if (id .gt. 0) then
                        marqi = imarq(jmarq(ic)+2*id-1)
                        if (marqi .eq. -1) then
                            call jjlide('JELIBZ', crnom, 2)
                            goto 171
                        end if
                    end if
                end do
171             continue
            else
!
! --------- OBJET SIMPLE
!
                iclas = ic
                iclaos = ic
                idatos = j
                nomos = crnom
                if (iclasi .ne. iclaos) then
                    nomco = d32
                    nomoc = d32
                end if
                marqi = imarq(jmarq(ic)+2*j-1)
                if (marqi .eq. -1) then
                    call jjlide('JELIBZ', crnom, 1)
                end if
            end if
150         continue
        end do
    end do
! FIN ------------------------------------------------------------------
end subroutine
