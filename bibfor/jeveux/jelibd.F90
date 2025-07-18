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
subroutine jelibd(nomlu, ltot)
!
! LIBERE DE LA MEMOIRE LE SEGMENT DE VALEURS ASSOCIES A UN OBJET
! JEVEUX ALLOUE DYNAMIQUEMENT LORSQU'IL EST DANS L'ETAT XA OU XD
!
! IN   NOMLU : NOM DE L'OBJET A LIBERER
! OUT  LTOT  : LONGUEUR EN ENTIERS LIBEREE
!
! ----------------------------------------------------------------------
! person_in_charge: j-pierre.lefebvre at edf.fr
! POUR LES VARS I4ZON, LSZON, R8ZON QUI DANS UN COMMON MAIS UNIQUEMENT
! PAR EQUIVALENCE
! aslint: disable=
    implicit none
#include "jeveux_private.h"
#include "asterfort/jjallc.h"
#include "asterfort/jjcroc.h"
#include "asterfort/jjlbsg.h"
#include "asterfort/jjvern.h"
#include "asterfort/utmess.h"
    character(len=*) :: nomlu
!
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!     -----------------------------------------------------------------
    integer(kind=8) :: iclas, iclaos, iclaco, idatos, idatco, idatoc
    common/iatcje/iclas, iclaos, iclaco, idatos, idatco, idatoc
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: iadmi, iadyn, ibacol, ibiadm
    integer(kind=8) :: inat, ixdeso, ixiadd, ixiadm, jcara
    integer(kind=8) :: jdate, jdocu, jgenr, jhcod, jiacce, jiadd, jiadm
    integer(kind=8) :: jindir, jlong, jlono, jltyp, jluti, jmarq, jorig
    integer(kind=8) :: jrnom, jtype, k, ltot, ltypi, n
    integer(kind=8) :: nbacce, nbmax
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
!
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
    common/jiacce/jiacce(n), nbacce(2*n)
    common/jindir/jindir(n)
    integer(kind=8) :: isstat
    common/iconje/isstat
    integer(kind=8) :: ldyn, lgdyn, nbdyn, nbfree
    common/idynje/ldyn, lgdyn, nbdyn, nbfree
    integer(kind=8) :: icdyn, mxltot
    common/xdynje/icdyn, mxltot
    real(kind=8) :: mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio, cuvtrav
    common/r8dyje/mxdyn, mcdyn, mldyn, vmxdyn, vmet, lgio(2), cuvtrav
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: datei
    common/iheuje/datei
! ----------------------------------------------------------------------
    integer(kind=8) :: ivnmax, iddeso, idiadd, idiadm
    parameter(ivnmax=0, iddeso=1, idiadd=2, idiadm=3)
! ----------------------------------------------------------------------
    character(len=32) :: noml32
    integer(kind=8) :: icre, iret
!
    noml32 = nomlu
    ltot = 0
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
        ltypi = ltyp(jltyp(iclaos)+idatos)
        if (iadmi .eq. 0) then
            goto 999
        end if
        iadyn = iadm(jiadm(iclaos)+2*idatos)
!
        call jjlbsg(iclaos, idatos, 0, 0, iadmi, &
                    iadyn, ltot)
!
    else if (iret .eq. 2) then
!
! ----- CAS D'UNE COLLECTION
!
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
! ----- CAS D'UNE COLLECTION ENTIERE : ON LIBERE TOUS LES OC
!
        ixiadd = iszon(jiszon+ibacol+idiadd)
        ixiadm = iszon(jiszon+ibacol+idiadm)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        ltypi = ltyp(jltyp(iclaco)+ixdeso)
        if (ixiadd .eq. 0) then
!
! ------- COLLECTION CONTIGUE
!
            iadmi = iadm(jiadm(iclaco)+2*ixdeso-1)
            if (iadmi .eq. 0) then
                goto 999
            end if
            iadyn = iadm(jiadm(iclaco)+2*ixdeso)
!
            call jjlbsg(iclaco, ixdeso, 0, 0, iadmi, &
                        iadyn, ltot)
!
        else
!
! ------- COLLECTION DISPERSEE
!
            nbmax = iszon(jiszon+ibacol+ivnmax)
            ibiadm = iadm(jiadm(iclaco)+2*ixiadm-1)
            do k = 1, nbmax
                iadmi = iszon(jiszon+ibiadm-1+2*k-1)
                if (iadmi .eq. 0) then
                    goto 10
                end if
                iadyn = iszon(jiszon+ibiadm-1+2*k)
!
                call jjlbsg(iclaco, idatco, k, ibacol, iadmi, &
                            iadyn, ltot)
!
10              continue
            end do
        end if
    else if (inat .eq. 3) then
!       ------ CAS D'UN OBJET DE COLLECTION  ------
        ixiadd = iszon(jiszon+ibacol+idiadd)
        ixiadm = iszon(jiszon+ibacol+idiadm)
        ixdeso = iszon(jiszon+ibacol+iddeso)
        ltypi = ltyp(jltyp(iclaco)+ixdeso)
        if (ixiadd .ne. 0) then
!
! -------- COLLECTION DISPERSEE
!
            ibiadm = iadm(jiadm(iclaco)+2*ixiadm-1)
            iadmi = iszon(jiszon+ibiadm-1+2*idatoc-1)
            if (iadmi .eq. 0) then
                goto 999
            end if
            iadyn = iszon(jiszon+ibiadm-1+2*idatoc)
!
            call jjlbsg(iclaco, idatco, idatoc, ibacol, iadmi, &
                        iadyn, ltot)
!
        end if
    end if
999 continue
!
end subroutine
