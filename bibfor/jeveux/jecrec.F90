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

subroutine jecrec(nomlu, listat, accelu, stoclu, longlu, &
                  nmax)
    implicit none
#include "jeveux_private.h"
#include "asterfort/jjanal.h"
#include "asterfort/jjcrec.h"
#include "asterfort/jjcren.h"
#include "asterfort/jjvern.h"
#include "asterfort/jxveuo.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: nomlu, listat, accelu, stoclu, longlu
    integer(kind=8), intent(in) :: nmax
!     ------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: iadcol, iadnum, icl, ilongu(1), ipgcex, jcara, jdate
    integer(kind=8) :: jdocu, jgenr, jhcod, jiadd, jiadm, jlong, jlongu
    integer(kind=8) :: jlono, jltyp, jluti, jmarq, jnom, jorig, jrnom
    integer(kind=8) :: jtype, la, lenk, n, nb, nbl
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!     ------------------------------------------------------------------
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    character(len=24) :: nomco
    character(len=32) :: nomuti, nomos, nomoc, bl32
    common/nomcje/nomuti, nomos, nomco, nomoc, bl32
!     ------------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
!     ------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idmarq, idnom, idlong, idlono
    integer(kind=8) :: idluti, idnum
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &               idmarq=4, idnom=5, idlong=7,&
     &               idlono=8, idluti=9, idnum=10)
!     ------------------------------------------------------------------
    character(len=8) :: stocka, cval(3)
    integer(kind=8) :: ibid, iadzon, knom(1), lval(3), icre, iret, iretc, icrei
    character(len=2) :: ta
    character(len=32) :: noml32, nom32
    character(len=24) :: nom24l
    character(len=26) :: nompar
! DEB ------------------------------------------------------------------
    ipgcex = ipgc
    nom24l = nomlu
    if (nmax .le. 0) then
        call utmess('F', 'JEVEUX_67', sk=nom24l)
    end if
!
    call jjanal(listat, 3, 3, lval, cval)
    iclas = index(classe, cval(1) (1:1))
    if (iclas .eq. 0) then
        call utmess('F', 'JEVEUX_68', sk=cval(1) (1:1))
    end if
!
    icre = 2
    call jjvern(nom24l//'        ', icre, iretc)
    icre = 1
!
    if (iretc .eq. 1) then
        call utmess('F', 'JEVEUX_69', sk=nom24l)
    else
        stocka = stoclu
        if (lval(3) .eq. 2) read (cval(3) (2:lval(3)), '(I1)') lenk
        if (lval(3) .eq. 3) read (cval(3) (2:lval(3)), '(I2)') lenk
        if (lval(3) .gt. 3) lenk = 512
        if (stocka .ne. 'CONTIG  ' .and. stocka .ne. 'DISPERSE') then
            call utmess('F', 'JEVEUX_70', sk=stocka)
        else if (longlu .ne. 'CONSTANT' .and. cval(2) (1:1) .eq. 'E') &
            then
            call utmess('F', 'JEVEUX_71', sk=nom24l)
        else if (stocka .eq. 'CONTIG  ' .and. longlu .ne. 'CONSTANT'&
 &            .and. cval(3) (1:1) .eq. 'K' .and. lenk .ne. 8&
 &            .and. lenk .ne. 16 .and. lenk .ne. 24) then
            call utmess('F', 'JEVEUX_72', sk=nom24l)
        end if
!
        call jjcrec(iclaco, idatco, 'X', 'I', idnum+1, &
                    iadcol)
        iszon(jiszon+iadcol+ivnmax) = nmax
!     ------------------------------------------------------------------
        nom32 = nom24l//'$$DESO  '
        call jjcren(nom32, icre, iret)
        call jjcrec(iclaos, idatos, cval(2) (1:1), cval(3) (1:lval(3)), 0, &
                    iadzon)
        if (cval(2) (1:1) .eq. 'E') then
            if (stocka .eq. 'CONTIG  ') then
                lono(jlono(iclaos)+idatos) = nmax
            else if (longlu .eq. 'CONSTANT') then
                lono(jlono(iclaos)+idatos) = 1
            end if
        end if
        iszon(jiszon+iadcol+iddeso) = idatos
!     ------------------------------------------------------------------
        nb = nmax+1
        if (stocka .eq. 'DISPERSE') then
            nb = nmax
            nom32 = nom24l//'$$IADD  '
            call jjcren(nom32, icre, iret)
            call jjcrec(iclaos, idatos, 'V', 'I', 2*nmax, &
                        ibid)
            iszon(jiszon+iadcol+idiadd) = idatos
!
            nom32 = nom24l//'$$IADM  '
            call jjcren(nom32, icre, iret)
            call jjcrec(iclaos, idatos, 'V', 'I', 2*nmax, &
                        ibid)
            iszon(jiszon+iadcol+idiadm) = idatos
!
            nom32 = nom24l//'$$MARQ  '
            call jjcren(nom32, icre, iret)
            call jjcrec(iclaos, idatos, 'V', 'I', 2*nmax, &
                        ibid)
            iszon(jiszon+iadcol+idmarq) = idatos
        end if
!     ------------------------------------------------------------------
        nompar = nom24l//'$$'
        if ((longlu .ne. 'CONSTANT' .and. longlu .ne. 'VARIABLE') .or. len(longlu) .ne. 8) then
            call utmess('F', 'JEVEUX_02')
            icl = iclaco
            if (len(longlu) .gt. 24) then
                call utmess('F', 'JEVEUX_73', sk=longlu)
            end if
            noml32 = longlu
            nompar = noml32(1:24)//'&&'
            icrei = 0
            call jjvern(noml32, icrei, iret)
            if (iret .eq. 0) then
                call jjvern(noml32, icre, iret)
                call jjcrec(iclaos, idatos, 'V', 'I', nmax, &
                            ibid)
            else if (iret .ne. 1) then
                call utmess('F', 'JEVEUX_73', sk=noml32)
            else
                if (icl .ne. iclaos) then
                    call utmess('F', 'JEVEUX_74', sk=noml32)
                end if
                nbl = long(jlong(iclaos)+idatos)
                if (nbl .lt. nmax) then
                    call utmess('F', 'JEVEUX_75', sk=noml32)
                else if (type(jtype(iclaos)+idatos) .ne. 'I') then
                    call utmess('F', 'JEVEUX_76', sk=noml32)
                end if
                ipgc = -1
                call jxveuo('E', ilongu, 1, jlongu)
                ipgc = ipgcex
            end if
!     ------------------------------------------------------------------
        else if (longlu .eq. 'VARIABLE' .and. len(longlu) .eq. 8) then
            nom32 = nom24l//'$$LONG  '
            call jjcren(nom32, icre, iret)
            call jjcrec(iclaos, idatos, 'V', 'I', nmax, &
                        ibid)
        end if
        if (longlu .ne. 'CONSTANT') then
            iszon(jiszon+iadcol+idlong) = idatos
!     ------------------------------------------------------------------
            nom32 = nompar//'LONO  '
            if (nom32(25:26) .eq. '&&') then
                icrei = 0
                call jjcren(nom32, icrei, iret)
                if (iret .eq. 0) then
                    call jjcren(nom32, icre, iret)
                    call jjcrec(iclaos, idatos, 'V', 'I', nb, &
                                ibid)
                end if
                ipgc = -1
                call jxveuo('E', ilongu, 1, jlongu)
                ipgc = ipgcex
            else
                call jjcren(nom32, icre, iret)
                call jjcrec(iclaos, idatos, 'V', 'I', nb, &
                            ibid)
            end if
            iszon(jiszon+iadcol+idlono) = idatos
!     ------------------------------------------------------------------
            nom32 = nompar//'LUTI  '
            if (nom32(25:26) .eq. '&&') then
                icrei = 0
                call jjcren(nom32, icrei, iret)
                if (iret .eq. 0) then
                    call jjcren(nom32, icre, iret)
                    call jjcrec(iclaos, idatos, 'V', 'I', nmax, &
                                ibid)
                end if
                ipgc = -1
                call jxveuo('E', ilongu, 1, jlongu)
                ipgc = ipgcex
            else
                call jjcren(nom32, icre, iret)
                call jjcrec(iclaos, idatos, 'V', 'I', nmax, &
                            ibid)
            end if
            iszon(jiszon+iadcol+idluti) = idatos
        end if
!     ------------------------------------------------------------------
        ta = accelu(1:2)
        if (index('NO $NU $', ta//' $') .eq. 0) then
            call utmess('F', 'JEVEUX_81', sk=ta)
        else
            la = len(accelu)
            if (la .gt. 3) then
                if (accelu(3:3) .ne. ' ') then
                    call utmess('F', 'JEVEUX_82', sk=accelu)
                end if
                if (la .gt. 28) then
                    call utmess('F', 'JEVEUX_83', sk=accelu)
                end if
                noml32 = accelu(4:min(la, len(noml32)))
            else
                noml32 = ' '
            end if
        end if
        if (ta .eq. 'NO' .and. noml32 .ne. bl32) then
            icl = iclaco
            icrei = 0
            call jjvern(noml32, icrei, iret)
            if (iret .eq. 0) then
                call jjvern(noml32, icre, iret)
                call jjcrec(iclaos, idatos, 'N', 'K8', nmax, &
                            ibid)
            else if (iret .ne. 1) then
                call utmess('F', 'JEVEUX_77', sk=noml32)
            else
                if (icl .ne. iclaos) then
                    call utmess('F', 'JEVEUX_78', sk=noml32)
                end if
                nbl = long(jlong(iclaos)+idatos)
                if (nbl .lt. nmax) then
                    call utmess('F', 'JEVEUX_79', sk=noml32)
                else if (genr(jgenr(iclaos)+idatos) .ne. 'N') then
                    call utmess('F', 'JEVEUX_80', sk=noml32)
                end if
                ipgc = -1
                call jxveuo('E', knom, 1, jnom)
                if (iadm(jiadm(iclaos)+2*idatos-1) .eq. 0) then
                    knom(jnom+4-1) = 0
                    knom(jnom+5-1) = 0
                end if
                ipgc = ipgcex
            end if
            iszon(jiszon+iadcol+idnom) = idatos
        else if (ta .eq. 'NO') then
            nom32 = nom24l//'$$NOM   '
            call jjcren(nom32, icre, iret)
            call jjcrec(iclaos, idatos, 'N', 'K8', nmax, &
                        ibid)
            iszon(jiszon+iadcol+idnom) = idatos
        else if (ta .eq. 'NU') then
            nom32 = nom24l//'$$NUM   '
            call jjcren(nom32, icre, iret)
            call jjcrec(iclaos, idatos, 'V', 'I', 2, &
                        iadnum)
            iszon(jiszon+iadnum) = nmax
            iszon(jiszon+iadcol+idnum) = idatos
        end if
    end if
! FIN ------------------------------------------------------------------
end subroutine
