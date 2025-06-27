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
subroutine rvtecn(releve, absc, itcopt, itsppt, coor, &
                  nomcmp, nomnoe, nbcmp, nbpoin, docu, &
                  nomtab, iocc, xnovar, ncheff, i1, &
                  ioc, isd)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsnopa.h"
#include "asterfort/rsorac.h"
#include "asterfort/rvtec2.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
!
    integer(kind=8) :: itcopt(*), itsppt(*), nbcmp, nbpoin, iocc, i1, ioc, isd
    real(kind=8) :: releve(*), absc(*), coor(*)
    character(len=4) :: docu
    character(len=8) :: nomcmp(*), nomnoe(*)
    character(len=16) :: ncheff
    character(len=19) :: nomtab
    character(len=24) :: xnovar
!     MISE EN TABLEAU POUR UN EXTRACTION SUR UN CHAM_NO
!     ------------------------------------------------------------------
! IN  : RELEVE : TABLE DES RELEVES DE VALEURS
! IN  : ABSC   : TABLE DES ABSCISSES DES POINTS
! IN  : ITCOPT : TABLE DES NOMBRES DE COUCHES PAR POINT
! IN  : ITSPPT : TABLE DES NOMBRES DE SOUS-PT PAR POINT
! IN  : COOR   : TABLE DES COORDONNEES DES POINTS
! IN  : NOMCMP : NOM DES CMP
! IN  : NOMNOE : TABLE DES EVENTUELS NOMS DE NOEUDS
! IN  : NBCMP  : NOMBRE DE CMP
! IN  : NBPOIN : NOMBRE DE POINT D'EVALUATION
! IN  : DOCU   : 'LSTN'
!     ------------------------------------------------------------------
    integer(kind=8) :: nbvari, nbpar, ilign, ipt, nbsp, nbco, lc, ln, is, ic, i2, valei(12)
    integer(kind=8) :: n1, adrval, adracc, jacc, ik, ir, ii, lcr, lck, nbacc, nbpr, jaces, iac, iadr
    integer(kind=8) :: iord(1)
    aster_logical :: exist
    real(kind=8) :: prec
    character(len=3) :: typpar
    character(len=8) :: k8b, acces, nomres, ctype, crit
    character(len=16) :: intitu
    character(len=24) :: nomval, nomacc, nnores, nomjv
    complex(kind=8) :: c16b
    character(len=80) :: valek(11)
    character(len=24), pointer :: para(:) => null()
    real(kind=8), pointer :: vale_r(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
    if (docu .ne. 'LSTN') goto 999
!
    call jelira(jexnum(xnovar, iocc), 'LONUTI', nbvari)
    if (nbvari .ne. 0) then
        call rvtec2(releve, absc, itcopt, itsppt, coor, &
                    nomnoe, nbcmp, nbpoin, docu, nomtab, &
                    iocc, xnovar, ncheff, i1)
        goto 999
    end if
!
    call getvtx('ACTION', 'INTITULE', iocc=iocc, scal=intitu, nbret=n1)
!
    call getvr8('ACTION', 'PRECISION', iocc=iocc, scal=prec, nbret=n1)
    call getvtx('ACTION', 'CRITERE', iocc=iocc, scal=crit, nbret=n1)
!
    nomval = ncheff//'.VALACCE'
    nomacc = ncheff//'.TYPACCE'
    nnores = ncheff//'.NOMRESU'
    call jeveuo(nomacc, 'L', jacc)
!
    lcr = nbcmp+10
    AS_ALLOCATE(vr=vale_r, size=lcr)
    lck = nbcmp+23
    AS_ALLOCATE(vk24=para, size=lck)
!
    ii = 0
    ir = 0
    ik = 1
    nbpar = 1
    valek(ik) = intitu
    para(nbpar) = 'INTITULE'
!
!
    if (zk8(jacc) .eq. 'DIRECT  ') then
        call jeveuo(jexnum(ncheff//'.LSCHEFF', 1), 'L', jacc)
        nbpar = nbpar+1
        para(nbpar) = 'CHAM_GD'
        ik = ik+1
        valek(ik) = zk24(jacc) (1:8)
    else
        call jeveuo(nnores, 'L', jacc)
        nomres = zk16(jacc) (1:8)
        nbpar = nbpar+1
        para(nbpar) = 'RESU'
        ik = ik+1
        valek(ik) = nomres
        nbpar = nbpar+1
        para(nbpar) = 'NOM_CHAM'
        ik = ik+1
        valek(ik) = zk16(jacc+1)
        call jeveuo(nomacc, 'L', adracc)
        call jeveuo(nomval, 'L', adrval)
        acces = zk8(adracc)
        if (acces(1:1) .eq. 'O') then
            nbpar = nbpar+1
            para(nbpar) = 'NUME_ORDRE'
            ii = ii+1
            valei(ii) = zi(adrval+i1-1)
            nomjv = '&&RVTECN.NOMS_ACCES'
            call rsnopa(nomres, 0, nomjv, nbacc, nbpr)
            if (nbacc .ne. 0) then
                call jeveuo(nomjv, 'L', jaces)
                do iac = 1, nbacc
                    call rsadpa(nomres, 'L', 1, zk16(jaces-1+iac), zi(adrval+i1-1), &
                                1, sjv=iadr, styp=ctype, istop=0)
                    call tbexip(nomtab, zk16(jaces-1+iac), exist, typpar)
                    if (.not. exist) then
                        call tbajpa(nomtab, 1, zk16(jaces-1+iac), ctype)
                    end if
                    nbpar = nbpar+1
                    para(nbpar) = zk16(jaces-1+iac)
                    if (ctype(1:1) .eq. 'I') then
                        ii = ii+1
                        valei(ii) = zi(iadr)
                    else if (ctype(1:1) .eq. 'R') then
                        ir = ir+1
                        vale_r(ir) = zr(iadr)
                    else if (ctype(1:3) .eq. 'K80') then
                        ik = ik+1
                        valek(ik) = zk80(iadr)
                    else if (ctype(1:3) .eq. 'K32') then
                        ik = ik+1
                        valek(ik) = zk32(iadr)
                    else if (ctype(1:3) .eq. 'K24') then
                        ik = ik+1
                        valek(ik) = zk24(iadr)
                    else if (ctype(1:3) .eq. 'K16') then
                        ik = ik+1
                        valek(ik) = zk16(iadr)
                    else if (ctype(1:2) .eq. 'K8') then
                        ik = ik+1
                        valek(ik) = zk8(iadr)
                    end if
                end do
                call jedetr(nomjv)
            end if
        else if (acces(1:1) .eq. 'M') then
            nbpar = nbpar+1
            para(nbpar) = 'NUME_ORDRE'
            call rsorac(nomres, 'NUME_MODE', zi(adrval+i1-1), 0.d0, k8b, &
                        c16b, prec, crit, iord, 1, &
                        n1)
            ii = ii+1
            valei(ii) = iord(1)
            nbpar = nbpar+1
            para(nbpar) = 'NUME_MODE'
            ii = ii+1
            valei(ii) = zi(adrval+i1-1)
        else if (acces(1:1) .eq. 'I') then
            nbpar = nbpar+1
            para(nbpar) = 'NUME_ORDRE'
            call rsorac(nomres, 'INST', 0, zr(adrval+i1-1), k8b, &
                        c16b, prec, crit, iord, 1, &
                        n1)
            ii = ii+1
            valei(ii) = iord(1)
            nbpar = nbpar+1
            para(nbpar) = 'INST'
            ir = ir+1
            vale_r(ir) = zr(adrval+i1-1)
        else if (acces(1:1) .eq. 'F') then
            nbpar = nbpar+1
            para(nbpar) = 'FREQ'
            ir = ir+1
            vale_r(ir) = zr(adrval+i1-1)
        end if
    end if
    if (docu .eq. 'LSTN') then
        call tbexip(nomtab, 'NOEUD', exist, typpar)
        if (.not. exist) then
            call tbajpa(nomtab, 1, 'NOEUD', 'K8')
        end if
        nbpar = nbpar+1
        para(nbpar) = 'NOEUD'
    end if
    nbpar = nbpar+1
    para(nbpar) = 'ABSC_CURV'
    nbpar = nbpar+1
    para(nbpar) = 'COOR_X'
    nbpar = nbpar+1
    para(nbpar) = 'COOR_Y'
    nbpar = nbpar+1
    para(nbpar) = 'COOR_Z'
    ic = 0
    is = 0
    do ipt = 1, nbpoin, 1
        nbsp = itsppt(ipt)
        nbco = itcopt(ipt)
        if (nbco .gt. 1 .and. ic .eq. 0) then
            call tbexip(nomtab, 'NUME_COUCHE', exist, typpar)
            if (.not. exist) then
                call tbajpa(nomtab, 1, 'NUME_COUCHE', 'I')
            end if
            ic = 1
            nbpar = nbpar+1
            para(nbpar) = 'NUME_COUCHE'
        end if
        if (nbsp .gt. 1 .and. is .eq. 0) then
            call tbexip(nomtab, 'NUME_GAUSS', exist, typpar)
            if (.not. exist) then
                call tbajpa(nomtab, 1, 'NUME_GAUSS', 'I')
            end if
            is = 1
            nbpar = nbpar+1
            para(nbpar) = 'NUME_GAUSS'
        end if
    end do
    do i2 = 1, nbcmp, 1
        nbpar = nbpar+1
        para(nbpar) = nomcmp(i2)
    end do
!
    lc = ir+4+nbcmp
    ASSERT(nbpar .le. lck)
    ASSERT(ii+2 .le. 10)
    ASSERT(ik+1 .le. 10)
    ASSERT(lc .le. lcr)
!
    ilign = 0
!
    ik = ik+1
    do ipt = 1, nbpoin, 1
!
        nbsp = itsppt(ipt)
        nbco = itcopt(ipt)
        lc = nbcmp*nbsp
        ln = lc*nbco
!
        if (docu .eq. 'LSTN') then
            valek(ik) = nomnoe(ipt)
        end if
!
        vale_r(ir+1) = absc(ipt)
        vale_r(1+ir+1) = coor(1+(ipt-1)*3)
        vale_r(1+ir+2) = coor(2+(ipt-1)*3)
        vale_r(1+ir+3) = coor(3+(ipt-1)*3)
!
        do ic = 1, nbco, 1
!
            valei(ii+1) = ic
!
            do is = 1, nbsp, 1
!
                valei(ii+2) = is
!
                do i2 = 1, nbcmp, 1
!
                    vale_r(1+ir+3+i2) = releve((ipt-1)*ln+lc*(ic-1)+nbcmp*(is-1)+i2)
!
                end do
!
                call tbajli(nomtab, nbpar, para, valei, vale_r, &
                            [c16b], valek, ilign)
!
            end do
!
        end do
!
    end do
!
    AS_DEALLOCATE(vr=vale_r)
    AS_DEALLOCATE(vk24=para)
!
999 continue
    call jedema()
end subroutine
