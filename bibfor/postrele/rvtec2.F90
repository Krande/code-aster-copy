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
subroutine rvtec2(releve, absc, itcopt, itsppt, coor, &
                  nomnoe, nbcmp, nbpoin, docu, nomtab, &
                  iocc, xnovar, ncheff, i1)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
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
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbexip.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: itcopt(*), itsppt(*), nbcmp, nbpoin, iocc, i1
    real(kind=8) :: releve(*), absc(*), coor(*)
    character(len=4) :: docu
    character(len=8) :: nomnoe(*)
    character(len=16) :: ncheff
    character(len=19) :: nomtab
    character(len=24) :: xnovar
!     MISE EN TABLEAU POUR UN EXTRACTION SUR UN CHAM_NO
!     ------------------------------------------------------------------
! IN  : RELEVE : TABLE DES RELEVE DE VALEURS
! IN  : ABSC   : TABLE DES ABSCISES DES POINTS
! IN  : ITCOPT : TABLE DES NOMBRES DE COUCHES PAR POINT
! IN  : ITSPPT : TABLE DES NOMBRES DE SOUS-PT PAR POINT
! IN  : COOR   : TABLE DES COORDONNEES DES POINTS
! IN  : NOMNOE : TABLE DES EVENTUELS NOMS DE NOEUDS
! IN  : NBCMP  : NOMBRE DE CMP
! IN  : NBPOIN : NOMBRE DE POINT D'EVALUATION
! IN  : DOCU   : 'LSTN'
!     ------------------------------------------------------------------
    integer(kind=8) :: nbvari, nbpar, ilign, ipt, nbsp, i, nbco, lc, ln, ic, i2, valei(12)
    integer(kind=8) :: n1, adrval, ind, lck, adracc, jacc, ik, ir, ii, ivari(3000), nbcmp2, jvari
    integer(kind=8) :: nbacc, nbpr, jaces, iac, iadr, lcr, iord(1)
    real(kind=8) :: prec
    aster_logical :: exist
    character(len=3) :: typpar
    character(len=7) :: kii
    character(len=8) :: k8b, acces, nomres, ctype, crit
    character(len=16) :: intitu
    character(len=24) :: nomval, nomacc, nnores, nomjv
    complex(kind=8) :: c16b
    character(len=80) :: valek(11)
    character(len=8), pointer :: nom_para(:) => null()
    character(len=24), pointer :: para(:) => null()
    character(len=8), pointer :: typ_para(:) => null()
    real(kind=8), pointer :: vale_r(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call jelira(jexnum(xnovar, iocc), 'LONUTI', nbvari)
    call jeveuo(jexnum(xnovar, iocc), 'L', jvari)
    if (nbvari .eq. 1 .and. zi(jvari) .eq. -1) then
        nbcmp2 = itsppt(1)
    else
        nbcmp2 = nbvari
    end if
    if (nbcmp2 .gt. 3000) then
        call utmess('F', 'POSTRELE_13')
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
    lcr = nbcmp2+10
    AS_ALLOCATE(vr=vale_r, size=lcr)
    lck = nbcmp2+23
    AS_ALLOCATE(vk24=para, size=lck)
    AS_ALLOCATE(vk8=nom_para, size=nbcmp2)
    AS_ALLOCATE(vk8=typ_para, size=nbcmp2)
!
    ik = 1
    ii = 0
    ir = 0
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
            nomjv = '&&RVTEC2.NOMS_ACCES'
            call rsnopa(nomres, 0, nomjv, nbacc, nbpr)
            if (nbacc .ne. 0) then
                call jeveuo(nomjv, 'L', jaces)
                do iac = 1, nbacc
                    call rsadpa(nomres, 'L', 1, zk16(jaces-1+iac), zi(adrval+i1-1), &
                                1, sjv=iadr, styp=ctype)
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
    do ipt = 1, nbpoin, 1
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
    end do
    if (nbvari .eq. 1 .and. zi(jvari) .eq. -1) then
        do i = 1, nbcmp2, 1
            ivari(i) = i
            call codent(i, 'G', kii)
            nbpar = nbpar+1
            para(nbpar) = 'V'//kii
            nom_para(i) = 'V'//kii
            typ_para(i) = 'R'
        end do
    else
        do i = 1, nbcmp2, 1
            ivari(i) = zi(jvari+i-1)
            call codent(zi(jvari+i-1), 'G', kii)
            nbpar = nbpar+1
            para(nbpar) = 'V'//kii
            nom_para(i) = 'V'//kii
            typ_para(i) = 'R'
        end do
    end if
    call tbajpa(nomtab, nbcmp2, nom_para, typ_para)
!
    lc = ir+4+nbcmp2
    ASSERT(nbpar .le. lck)
    ASSERT(ii+2 .le. 10)
    ASSERT(lc .le. lcr)
    ASSERT(ik .le. 10)
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
            do i2 = 1, nbcmp2, 1
                ind = (ipt-1)*ln+lc*(ic-1)+ivari(i2)
                vale_r(1+ir+3+i2) = releve(ind)
            end do
!
            call tbajli(nomtab, nbpar, para, valei, vale_r, &
                        [c16b], valek, ilign)
!
        end do
!
    end do
!
    AS_DEALLOCATE(vr=vale_r)
    AS_DEALLOCATE(vk24=para)
    AS_DEALLOCATE(vk8=nom_para)
    AS_DEALLOCATE(vk8=typ_para)
!
    call jedema()
end subroutine
