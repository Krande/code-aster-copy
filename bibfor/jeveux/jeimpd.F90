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
subroutine jeimpd(unit, clas, cmess)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjlide.h"
    integer(kind=8) :: unit
    character(len=*) :: clas, cmess
! ---------------------------------------------------------------------
! ROUTINE UTILISATEUR D'IMPRESSION DE LA LISTE DES OBJETS PRESENTS SUR
! LE FICHIER D'ACCES DIRECT ASSOCIE A UNE BASE
!
! IN  UNIT  : NUMERO D'UNITE LOGIQUE ASSOCIE AU FICHIER D'IMPRESSION
! IN  CLAS   : CLASSE ASSOCIEE A LA BASE ( ' ' : TOUTES LES CLASSES )
! IN  CMESS  : MESSAGE D'INFORMATION
! ---------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ---------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacc, ibacol, ibiadd, iblono, iiadd
    integer(kind=8) :: ilono, iltyp, ixdeso, ixiadd, ixlono, j, jcara
    integer(kind=8) :: jdate, jdocu, jgenr, jhcod, jiacce, jiadd, jiadm
    integer(kind=8) :: jlong, jlono, jltyp, jluti, jmarq, jorig, jrnom
    integer(kind=8) :: jtype, k, kiadd, kj, koc, liadd, n
    integer(kind=8) :: nbacce, ncla1, ncla2, nmax
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
!
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    character(len=8) :: nombas
    common/kbasje/nombas(n)
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
! ---------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/jiacce/jiacce(n), nbacce(2*n)
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
! ---------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idlono
    parameter(ivnmax=0, iddeso=1, idiadd=2,&
     &               idlono=8)
! ---------------------------------------------------------------------
    character(len=1) :: kclas, cgenr, ctype, clasi, cgen2
    character(len=32) :: crnom
    aster_logical :: lcol, lente
    integer(kind=8) :: ipgcex, lgbl
! DEB -----------------------------------------------------------------
    ipgcex = ipgc
    ipgc = -2
    lente = .true.
    kclas = clas(1:min(1, len(clas)))
    if (unit .le. 0) goto 999
    if (kclas .eq. ' ') then
        ncla1 = 1
        ncla2 = index(classe, '$')-1
        if (ncla2 .lt. 0) ncla2 = n
    else
        ncla1 = index(classe, kclas)
        ncla2 = ncla1
    end if
    do i = ncla1, ncla2
        clasi = classe(i:i)
        if (clasi .ne. ' ') then
            write (unit, '(''1'',4A)') ('--------------------', k=1, 4)
            write (unit, *) '                                  '
            write (unit, '(1X,2A)')&
     &          '       CONTENU DE LA BASE ', clasi,&
     &          '        ', cmess(1:min(72, len(cmess)))
            write (unit, *) ' NOM DE LA BASE               : ', nombas(i)
            write (unit, *) ' NB D''ENREGISTREMENTS MAXIMUM : ', nblmax(i)
            lgbl = 1024*longbl(i)*lois
            write (unit, *) ' LONGUEUR D''ENREGISTREMENT (OCTETS): ', &
                lgbl
            write (unit, *) '                                  '
            write (unit, '(    1X,4A)') ('--------------------', k=1, 4)
            kj = 1
            do j = 1, nremax(i)
                crnom = rnom(jrnom(i)+j)
                if (crnom(1:1) .eq. '?') goto 5
                if (mod(kj, 25) .eq. 1 .and. lente) then
                    write (unit, '(/,A,A/)')&
     &     '---- NUM ------------- NOM ---------------- G T -L-'&
     &     , ' -LOTY- -IADD- --LIADD- NB AC'
                    lente = .false.
                end if
                cgenr = genr(jgenr(i)+j)
                ctype = type(jtype(i)+j)
                iltyp = ltyp(jltyp(i)+j)
                ilono = lono(jlono(i)+j)
                iiadd = iadd(jiadd(i)+2*j-1)
                if (iiadd .eq. 0) goto 6
                kj = kj+1
                lente = .true.
                lcol = .false.
                liadd = iadd(jiadd(i)+2*j)
                iacc = iacce(jiacce(i)+iiadd)
                write (unit, 1001) j, crnom, cgenr, ctype, iltyp, ilono, &
                    iiadd, liadd, iacc
6               continue
                if (cgenr .eq. 'X') then
                    idatco = j
                    iclaco = i
                    lcol = .true.
                    call jjallc(i, j, 'L', ibacol)
                    ixiadd = iszon(jiszon+ibacol+idiadd)
                    ixdeso = iszon(jiszon+ibacol+iddeso)
                    if (ixiadd .eq. 0) goto 51
                    ixlono = iszon(jiszon+ibacol+idlono)
                    nmax = iszon(jiszon+ibacol+ivnmax)
                    cgen2 = genr(jgenr(i)+ixdeso)
                    ctype = type(jtype(i)+ixdeso)
                    iltyp = ltyp(jltyp(i)+ixdeso)
                    do koc = 1, nmax
                        ibiadd = iadm(jiadm(i)+2*ixiadd-1)
                        kiadd = iszon(jiszon+ibiadd-1+2*koc-1)
                        if (kiadd .eq. 0) goto 50
                        if (mod(kj, 25) .eq. 1 .and. lente) then
                            write (unit, '(/,A,A/)')&
     &            '---- NUM ------------- NOM -------------- E G T -L-'&
     &           , ' -LOTY- -IADD- --LIADD- NB AC'
                            lente = .false.
                        end if
                        iiadd = iszon(jiszon+ibiadd-1+2*koc-1)
                        liadd = iszon(jiszon+ibiadd-1+2*koc)
                        iacc = iacce(jiacce(i)+iiadd)
                        if (ixlono .eq. 0) then
                            ilono = lono(jlono(i)+ixdeso)
                        else
                            iblono = iadm(jiadm(i)+2*ixlono-1)
                            ilono = iszon(jiszon+iblono-1+koc)
                        end if
                        kj = kj+1
                        lente = .true.
                        write (crnom(25:32), '(I8)') koc
                        write (unit, 1001) j, crnom, cgen2, ctype, iltyp, &
                            ilono, iiadd, liadd, iacc
50                      continue
                    end do
51                  continue
                    if (lcol) then
                        call jjlide('JEIMPO', crnom(1:24), 2)
                    end if
                end if
5               continue
            end do
            write (unit, '(/)')
        end if
    end do
999 continue
    ipgc = ipgcex
1001 format(i8, 1x, a, '  -', 2(a, '-'), i3, i7, i7, i9, i6)
! FIN -----------------------------------------------------------------
end subroutine
