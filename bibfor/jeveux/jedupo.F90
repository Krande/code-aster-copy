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
subroutine jedupo(schin, claout, schout, dupcol)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/assert.h"
#include "asterfort/jedetr.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjallt.h"
#include "asterfort/jjcrec.h"
#include "asterfort/jjcren.h"
#include "asterfort/jjecrs.h"
#include "asterfort/jjlide.h"
#include "asterfort/jxdeps.h"
#include "asterfort/jxliro.h"
#include "asterfort/utmess.h"
    character(len=*) :: schin, claout, schout
    aster_logical :: dupcol
! ----------------------------------------------------------------------
!     RECOPIE L'OBJET DE NOM SCHIN LA CLASSE DANS
!     SCHOUT AVEC LE NOM DISTINCT SCHOUT
!
! IN  SCHIN  : SOUS-CHAINE EN ENTREE
! IN  CLAOUT : NOM DE LA CLASSE EN SORTIE
! IN  SCHOUT : SOUS-CHAINE EN SORTIE
! IN  DUPCOL : .TRUE. DUPLIQUE LES OBJETS PARTAGEABLES D'UNE COLLECTION
!              .FALSE. S'ARRETE SUR ERREUR
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmi, iadmo1, iadmo2, iadout, iadyn, iadzon, ibacol
    integer(kind=8) :: ibaout, ibiadd, ibiadm, ibiado, iblono, ibmaro, icin
    integer(kind=8) :: icout, idat, idcout, idin, idout, iret1, iret2
    integer(kind=8) :: iret3, ista1, ista2, ixdeso, ixiadd, ixiadm, ixiado
    integer(kind=8) :: ixlono, ixmaro, jcara, jdate, jdocu, jgenr, jhcod
    integer(kind=8) :: jiadd, jiadm, jlong, jlono, jltyp, jluti, jmarq
    integer(kind=8) :: jorig, jrnom, jtype, k, lonoi, n, nbl
    integer(kind=8) :: nmax
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
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: istat
    common/istaje/istat(4)
    integer(kind=8) :: numec
    common/inumje/numec
    character(len=24) :: nomec
    common/knomje/nomec
!     ------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm, idmarq, idlono, idnum
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3,&
     &               idmarq=4,&
     &               idlono=8, idnum=10)
    integer(kind=8) :: iv(idnum)
    character(len=8) :: csuffi(idnum)
! ----------------------------------------------------------------------
    integer(kind=8) :: ltypi, iaddi(2)
    character(len=32) :: nomin, nomout, nomcol
    character(len=1) :: kclas, genri, typei
    aster_logical :: idenba, libcol, x2u
    data iv/0, 0, 0, 0, 1, 0, 1, 1, 1, 1/
    data csuffi/'$$DESO  ', '$$IADD  ', '$$IADM  ',&
     &                          '$$MARQ  ', '$$NOM   ', '        ',&
     &                          '$$LONG  ', '$$LONO  ', '$$LUTI  ',&
     &                          '$$NUM   '/
! DEB ------------------------------------------------------------------
    kclas = claout
    icout = index(classe, kclas)
!
    nomin = schin
    nomout = schout
    call jjcren(nomin(1:24), 0, iret1)
    if (iret1 .eq. 1) then
        icin = iclaos
        idin = idatos
        iclas = icout
        call jjcren(nomout(1:24), 0, iret2)
        if (iret2 .ne. 0) then
            call jedetr(nomout(1:24))
        end if
        call jjcren(nomout(1:24), 1, iret2)
        idout = idatos
        genr(jgenr(icout)+idout) = genr(jgenr(icin)+idin)
        ltyp(jltyp(icout)+idout) = ltyp(jltyp(icin)+idin)
        type(jtype(icout)+idout) = type(jtype(icin)+idin)
        lono(jlono(icout)+idout) = lono(jlono(icin)+idin)
        long(jlong(icout)+idout) = long(jlong(icin)+idin)
        genri = genr(jgenr(icout)+idout)
        typei = type(jtype(icout)+idout)
        ltypi = ltyp(jltyp(icout)+idout)
        lonoi = lono(jlono(icout)+idout)
        nbl = lonoi*ltypi
        iadmi = iadm(jiadm(icin)+2*idin-1)
!
        if (iadmi .ne. 0) then
            ista1 = iszon(jiszon+iadmi-1)
            ista2 = iszon(jiszon+iszon(jiszon+iadmi-4)-4)
            if (ista1 .eq. istat(1) .and. (ista2 .eq. istat(3) .or. ista2 .eq. istat(4))) then
                x2u = .true.
                iszon(jiszon+iadmi-1) = istat(2)
            else
                x2u = .false.
            end if
        end if
        call jjallt(nbl, icout, genri, typei, ltypi, &
                    'INIT', iadout, iadyn)
        call jjecrs(iadout, icout, idout, 0, 'E', &
                    imarq(jmarq(icout)+2*idout-1))
        iadm(jiadm(icout)+2*idout-1) = iadout
        iadm(jiadm(icout)+2*idout) = iadyn
        iadmi = iadm(jiadm(icin)+2*idin-1)
        iaddi(1) = iadd(jiadd(icin)+2*idin-1)
        iaddi(2) = iadd(jiadd(icin)+2*idin)
        if (iadmi .ne. 0) then
            iadmo1 = (iadmi-1)*lois+iszon(jiszon+iadmi-3)+1
            iadmo2 = (iadout-1)*lois+iszon(jiszon+iadout-3)+1
            call jxdeps(iadmo1, iadmo2, lonoi*ltypi)
            if (x2u) iszon(jiszon+iadmi-1) = istat(1)
        else if (iaddi(1) .gt. 0) then
            call jxliro(icin, iadout, iaddi, lonoi*ltypi)
        else
            call utmess('F', 'JEVEUX1_66', sk=nomin)
        end if
        docu(jdocu(icout)+idout) = docu(jdocu(icin)+idin)
        luti(jluti(icout)+idout) = luti(jluti(icin)+idin)
        call jjlide('JELIBE', nomout(1:24), iret2)
    else if (iret1 .eq. 2) then
        nomcol = nomin
        idin = idatco
        icin = iclaco
        ibacol = iadm(jiadm(icin)+2*idin-1)
        libcol = .false.
        if (ibacol .eq. 0) then
            libcol = .true.
        else
            if (iszon(jiszon+ibacol-1) .eq. istat(1)) libcol = .true.
        end if
        call jjallc(icin, idin, 'L', ibacol)
        nmax = iszon(jiszon+ibacol+ivnmax)
        iclas = icout
        call jjcren(nomout(1:24), 0, iret2)
        if (iret2 .ne. 0) then
            call jedetr(nomout(1:24))
        end if
        call jjcren(nomout(1:24), 2, iret2)
        idout = idatco
        idcout = idatco
        call jjcrec(icout, idout, 'X', 'I', idnum+1, &
                    iadzon)
        iszon(jiszon+iadzon) = nmax
        ibaout = iadm(jiadm(icout)+2*idout-1)
        if (iszon(jiszon+ibacol+2) .eq. 0) then
            iv(1) = 1
        else
            iv(1) = 0
        end if
!
! ----- RECOPIE DES OBJETS ATTRIBUTS DE COLLECTION
!
        idenba = icin .eq. icout
        do k = 1, idnum
            idat = iszon(jiszon+ibacol+k)
            if (idat .gt. 0) then
                nomin = rnom(jrnom(icin)+idat)
                if (nomin(1:24) .ne. nomcol(1:24)) then
                    if (k .eq. 5) then
                        ASSERT(.false.)
                    end if
                    if (dupcol .or. .not. idenba) then
                        nomout = nomout(1:24)//csuffi(k)
                    else
                        nomout = nomin
                    end if
                else
                    nomout = nomout(1:24)//nomin(25:32)
                end if
                call jjcren(nomout(1:32), 0, iret3)
                if (iret3 .eq. 0) then
                    call jjcren(nomout(1:32), 1, iret3)
                    idout = idatos
                    iadmi = iadm(jiadm(icin)+2*idat-1)
                    iaddi(1) = iadd(jiadd(icin)+2*idat-1)
                    iaddi(2) = iadd(jiadd(icin)+2*idat)
                    genr(jgenr(icout)+idout) = genr(jgenr(icin)+idat)
                    ltyp(jltyp(icout)+idout) = ltyp(jltyp(icin)+idat)
                    type(jtype(icout)+idout) = type(jtype(icin)+idat)
                    lono(jlono(icout)+idout) = lono(jlono(icin)+idat)
                    long(jlong(icout)+idout) = long(jlong(icin)+idat)
                    genri = genr(jgenr(icout)+idout)
                    typei = type(jtype(icout)+idout)
                    ltypi = ltyp(jltyp(icout)+idout)
                    lonoi = lono(jlono(icout)+idout)
                    nbl = lonoi*ltypi
                    if (nbl .gt. 0) then
                        if ((k .eq. 1 .and. iv(1) .eq. 1) .or. k .gt. 1) then
                            call jjallt(nbl, icout, genri, typei, ltypi, &
                                        'INIT', iadout, iadyn)
                            call jjecrs(iadout, icout, idout, 0, 'E', &
                                        imarq(jmarq(icout)+2*idout-1))
                            iadm(jiadm(icout)+2*idout-1) = iadout
                            iadm(jiadm(icout)+2*idout) = iadyn
                        end if
                    end if
                    if (iv(k) .eq. 1) then
                        if (iadmi .ne. 0) then
                            iadmo1 = (iadmi-1)*lois+iszon(jiszon+iadmi-3)+1
                            iadmo2 = (iadout-1)*lois+iszon(jiszon+iadout-3)+1
                            call jxdeps(iadmo1, iadmo2, lonoi*ltypi)
                        else if (iaddi(1) .gt. 0) then
                            call jxliro(icin, iadout, iaddi, lonoi*ltypi)
                        end if
                    end if
                    docu(jdocu(icout)+idout) = docu(jdocu(icin)+idat)
                    luti(jluti(icout)+idout) = luti(jluti(icin)+idat)
                end if
                iszon(jiszon+ibaout+k) = idatos
            end if
        end do
!
! ----- POUR UNE COLLECTION DISPERSEE, RECOPIE DES SEGMENTS DE VALEURS
! ----- ASSOCIES AUX OBJETS DE COLLECTION
!
        if (iv(1) .eq. 0) then
            ixdeso = iszon(jiszon+ibacol+iddeso)
            ixiadm = iszon(jiszon+ibacol+idiadm)
            ibiadm = iadm(jiadm(icin)+2*ixiadm-1)
!
            ixmaro = iszon(jiszon+ibaout+idmarq)
            ibmaro = iadm(jiadm(icout)+2*ixmaro-1)
!
            ixiado = iszon(jiszon+ibaout+idiadm)
            ibiado = iadm(jiadm(icout)+2*ixiado-1)
            ixiadd = iszon(jiszon+ibacol+idiadd)
            ibiadd = iadm(jiadm(icin)+2*ixiadd-1)
            ixlono = iszon(jiszon+ibacol+idlono)
            genri = genr(jgenr(icin)+ixdeso)
            typei = type(jtype(icin)+ixdeso)
            ltypi = ltyp(jltyp(icin)+ixdeso)
            do k = 1, nmax
                iadmi = iszon(jiszon+ibiadm-1+2*k-1)
                iaddi(1) = iszon(jiszon+ibiadd-1+2*k-1)
                iaddi(2) = iszon(jiszon+ibiadd-1+2*k)
                if (iadmi .eq. 0 .and. iaddi(1) .eq. 0) goto 2
                if (ixlono .eq. 0) then
                    nbl = lono(jlono(icin)+ixdeso)*ltypi
                else
                    iblono = iadm(jiadm(icin)+2*ixlono-1)
                    nbl = iszon(jiszon+iblono-1+k)*ltypi
                end if
                if (iadmi .ne. 0) then
                    ista1 = iszon(jiszon+iadmi-1)
                    ista2 = iszon(jiszon+iszon(jiszon+iadmi-4)-4)
                  if (ista1 .eq. istat(1) .and. (ista2 .eq. istat(3) .or. ista2 .eq. istat(4))) then
                        x2u = .true.
                        iszon(jiszon+iadmi-1) = istat(2)
                    else
                        x2u = .false.
                    end if
                end if
                call jjallt(nbl, icout, genri, typei, ltypi, &
                            'INIT', iadout, iadyn)
                call jjecrs(iadout, icout, k, idcout, 'E', &
                            iszon(jiszon+ibmaro-1+2*k-1))
                iszon(jiszon+ibiado-1+2*k-1) = iadout
                iszon(jiszon+ibiado-1+2*k) = iadyn
                numec = k
                nomec = ' '
                if (iadmi .ne. 0) then
                    iadmo1 = (iadmi-1)*lois+iszon(jiszon+iadmi-3)+1
                    iadmo2 = (iadout-1)*lois+iszon(jiszon+iadout-3)+1
                    call jxdeps(iadmo1, iadmo2, nbl)
                    if (x2u) iszon(jiszon+iadmi-1) = istat(1)
                    call jjlide('JELIBE', nomin(1:24)//'$$XNUM  ', 2)
                else if (iaddi(1) .gt. 0) then
                    call jxliro(icin, iadout, iaddi, nbl)
                else
                    call utmess('F', 'JEVEUX1_65', sk=nomin, si=k)
                end if
                call jjlide('JELIBE', nomout(1:24)//'$$XNUM  ', 2)
2               continue
            end do
        end if
        if (libcol) call jjlide('JELIBE', nomin(1:24), iret1)
        call jjlide('JELIBE', nomout(1:24), iret2)
    end if
! FIN ------------------------------------------------------------------
end subroutine
