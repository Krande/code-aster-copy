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
subroutine jeimpo(unit, nomlu, mess)
! person_in_charge: j-pierre.lefebvre at edf.fr
! IMPRIME LE CONTENU D'UN OBJET JEVEUX
!
! IN  UNIT  : NUMERO D'UNITE LOGIQUE ASSOCIE AU FICHIER D'IMPRESSION
! IN  NOMLU : NOM DE L'OBJET JEVEUX A IMPRIMER
! IN  MESS  : MESSAGE D'INFORMATION IMPRIME
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjalty.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjimpo.h"
#include "asterfort/jjlide.h"
#include "asterfort/jjvern.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: unit
    character(len=*) :: nomlu, mess
!     ------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!     -----------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmex, iadmi, ibdeso, ibiadd, ibiadm, iblono, ideci
    integer(kind=8) :: inat, ipgcex, ixiadm, ixlono, jcara, jdate, jdocu
    integer(kind=8) :: jgenr, jhcod, jiadd, jiadm, jlong, jlono, jltyp
    integer(kind=8) :: jluti, jmarq, jorig, jrnom, jtype, k, n
    integer(kind=8) :: nbmax
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: ipgc, kdesma(2), lgd, lgduti, kposma(2), lgp, lgputi
    common/iadmje/ipgc, kdesma, lgd, lgduti, kposma, lgp, lgputi
    integer(kind=8) :: numec
    common/inumje/numec
!     ------------------------------------------------------------------
    character(len=32) :: noml32
    character(len=1) :: genri, typei
    integer(kind=8) :: icre, iret, jctab, ltypi, lonoi, iaddi(2)
    integer(kind=8) :: ibacol, ixiadd, ixdeso
    aster_logical :: lconst, lcol
!     ------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idlong, idlono
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &                              idlong=7,&
     &               idlono=8)
! DEB ------------------------------------------------------------------
    ipgcex = ipgc
    ipgc = -2
    noml32 = nomlu
!
    lcol = .false.
    icre = 0
    call jjvern(noml32, icre, iret)
!
    if (iret .eq. 0) then
        call utmess('A', 'JEVEUX_26', sk=noml32(1:24))
        goto 999
    else if (iret .eq. 1) then
!
! ----  CAS D'UN OBJET SIMPLE
!
        inat = 1
        iadmi = iadm(jiadm(iclaos)+2*idatos-1)
        iadmex = iadmi
        genri = genr(jgenr(iclaos)+idatos)
        typei = type(jtype(iclaos)+idatos)
        ltypi = ltyp(jltyp(iclaos)+idatos)
        lonoi = lono(jlono(iclaos)+idatos)*ltypi
        if (iadmex .eq. 0) then
            iaddi(1) = iadd(jiadd(iclaos)+2*idatos-1)
            iaddi(2) = iadd(jiadd(iclaos)+2*idatos)
            if (iaddi(1) .eq. 0) then
                call utmess('A', 'JEVEUX_27', sk=noml32(1:24))
                goto 999
            end if
            call jjalty(typei, ltypi, 'L', 1, jctab)
            iadmi = iadm(jiadm(iclaos)+2*idatos-1)
        end if
        ideci = 0
        call jjimpo(unit, iadmi, ideci, 0, genri, &
                    typei, ltypi, lonoi, mess)
        if (iadmex .eq. 0) then
            call jjlide('JEIMPO', noml32, inat)
        end if
    else if (iret .eq. 2) then
!
! ----- CAS D'UNE COLLECTION
!
        lcol = .true.
        call jjallc(iclaco, idatco, 'L', ibacol)
        if (noml32(25:32) .eq. '        ') then
            inat = 2
        else
            call jjcroc(noml32(25:32), icre)
            inat = 3
        end if
    end if
    if (inat .eq. 2) then
!
! ----- CAS D'UNE COLLECTION ENTIERE
!
        ixiadd = iszon(jiszon+ibacol+idiadd)
        ixiadm = iszon(jiszon+ibacol+idiadm)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        genri = genr(jgenr(iclaco)+ixdeso)
        typei = type(jtype(iclaco)+ixdeso)
        ltypi = ltyp(jltyp(iclaco)+ixdeso)
        if (ixiadd .eq. 0) then
!
! ------- COLLECTION CONTIGUE
!
            iadmi = iadm(jiadm(iclaco)+2*ixdeso-1)
            iaddi(1) = iadd(jiadd(iclaco)+2*ixdeso-1)
            iaddi(2) = iadd(jiadd(iclaco)+2*ixdeso)
            iadmex = iadmi
            if (iadmex .eq. 0) then
                if (iaddi(1) .eq. 0) then
                    call utmess('A', 'JEVEUX_28', sk=noml32(1:24))
                    goto 999
                end if
                call jjalty(typei, ltypi, 'L', 2, jctab)
                iadmi = iadm(jiadm(iclaco)+2*ixdeso-1)
            end if
            lonoi = lono(jlono(iclaco)+ixdeso)*ltypi
            ideci = 0
            call jjimpo(unit, iadmi, ideci, -1, genri, &
                        typei, ltypi, lonoi, mess)
            if (iadmex .eq. 0) then
                call jjlide('JEIMPO', noml32, inat)
            end if
        else
!
! ------- COLLECTION DISPERSEE
!
            nbmax = iszon(jiszon+ibacol+ivnmax)
            ibiadm = iadm(jiadm(iclaco)+2*ixiadm-1)
            ibiadd = iadm(jiadm(iclaco)+2*ixiadd-1)
            ideci = 0
            do k = 1, nbmax
                iadmi = iszon(jiszon+ibiadm-1+2*k-1)
                if (iadmi .eq. 0) then
                    iaddi(1) = iszon(jiszon+ibiadd-1+2*k-1)
                    iaddi(2) = iszon(jiszon+ibiadd-1+2*k)
                    if (iaddi(1) .eq. 0) goto 10
                    call jjalty(typei, ltypi, 'L', 3, jctab)
                    iadmi = iszon(jiszon+ibiadm-1+2*k-1)
                end if
                ixlono = iszon(jiszon+ibacol+idlono)
                if (ixlono .eq. 0) then
                    lonoi = lono(jlono(iclaco)+ixdeso)*ltypi
                else
                    iblono = iadm(jiadm(iclaco)+2*ixlono-1)
                    lonoi = iszon(jiszon+iblono-1+k)*ltypi
                end if
                call jjimpo(unit, iadmi, ideci, k, genri, &
                            typei, ltypi, lonoi, mess)
                numec = k
                call jjlide('JEIMPO', noml32//'$$XNUM  ', 2)
10              continue
            end do
        end if
        call jjlide('JEIMPO', noml32, inat)
    else if (inat .eq. 3) then
!       ------ CAS D'UN OBJET DE COLLECTION  ------
        ixiadd = iszon(jiszon+ibacol+idiadd)
        ixiadm = iszon(jiszon+ibacol+idiadm)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        ixlono = iszon(jiszon+ibacol+idlono)
        genri = genr(jgenr(iclaco)+ixdeso)
        typei = type(jtype(iclaco)+ixdeso)
        ltypi = ltyp(jltyp(iclaco)+ixdeso)
        if (ixiadd .eq. 0) then
!           ----------- COLLECTION CONTIGUE
            lconst = (iszon(jiszon+ibacol+idlong) .eq. 0)
            ibdeso = iadm(jiadm(iclaco)+2*ixdeso-1)
            iaddi(1) = iadd(jiadd(iclaco)+2*ixdeso-1)
            iaddi(2) = iadd(jiadd(iclaco)+2*ixdeso)
            iadmex = ibdeso
            if (iadmex .eq. 0) then
                if (iaddi(1) .eq. 0) then
                    call utmess('A', 'JEVEUX_29', sk=noml32(1:24), si=idatoc)
                    goto 999
                end if
                call jjalty(typei, ltypi, 'L', 2, jctab)
                ibdeso = iadm(jiadm(iclaco)+2*ixdeso-1)
            end if
            if (lconst) then
                lonoi = lono(jlono(iclaco)+ixdeso)*ltypi
                lonoi = lonoi/iszon(jiszon+ibacol+ivnmax)
                iadmi = ibdeso
                ideci = (idatoc-1)*lonoi
            else
                iblono = iadm(jiadm(iclaco)+2*ixlono-1)
                lonoi = ltypi*(iszon(jiszon+iblono-1+idatoc+1)-iszon(jiszon+iblono-1+idatoc) &
                               )
                iadmi = ibdeso
                ideci = (ltypi*(iszon(jiszon+iblono-1+idatoc)-1))
            end if
            call jjimpo(unit, iadmi, ideci, idatoc, genri, &
                        typei, ltypi, lonoi, mess)
            if (iadmex .eq. 0) then
                call jjlide('JEIMPO', noml32, inat)
            end if
        else
!
! -------- COLLECTION DISPERSEE
!
            ibiadm = iadm(jiadm(iclaco)+2*ixiadm-1)
            ibiadd = iadm(jiadm(iclaco)+2*ixiadd-1)
            iadmi = iszon(jiszon+ibiadm-1+2*idatoc-1)
            iadmex = iadmi
            ideci = 0
            if (iadmex .eq. 0) then
                iaddi(1) = iszon(jiszon+ibiadd-1+2*idatoc-1)
                iaddi(2) = iszon(jiszon+ibiadd-1+2*idatoc)
                if (iaddi(1) .eq. 0) then
                    call utmess('A', 'JEVEUX_29', sk=noml32(1:24), si=idatoc)
                    goto 999
                end if
                call jjalty(typei, ltypi, 'L', inat, jctab)
                iadmi = iszon(jiszon+ibiadm-1+2*idatoc-1)
            end if
            ixlono = iszon(jiszon+ibacol+idlono)
            if (ixlono .eq. 0) then
                lonoi = lono(jlono(iclaco)+ixdeso)*ltypi
            else
                iblono = iadm(jiadm(iclaco)+2*ixlono-1)
                lonoi = iszon(jiszon+iblono+idatoc-1)*ltypi
            end if
            call jjimpo(unit, iadmi, ideci, idatoc, genri, &
                        typei, ltypi, lonoi, mess)
            if (iadmex .eq. 0) then
                call jjlide('JEIMPO', noml32, inat)
            end if
        end if
    end if
999 continue
    if (lcol) then
        call jjlide('JEIMPO', noml32(1:24), 2)
    end if
    ipgc = ipgcex
! FIN ------------------------------------------------------------------
end subroutine
