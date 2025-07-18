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
subroutine rvtec0(t, co, sp, absc, x, &
                  cmp, nd, sdm, nbpoin, docu, &
                  nbcmp, padr, nomtab, iocc, xnovar, &
                  ncheff, i1)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
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
    integer(kind=8) :: co(*), sp(*), nbpoin, nbcmp, padr(*), iocc, i1
    real(kind=8) :: t(*), absc(*), x(*)
    character(len=4) :: docu
    character(len=8) :: cmp(*), nd(*)
    character(len=16) :: ncheff
    character(len=19) :: nomtab
    character(len=24) :: sdm, xnovar
!     AFFICHAGE CHAM_ELEM DE NBCMP COMPOSANTES
!     COPIE ROUTINE RVIMPK
!     ------------------------------------------------------------------
! IN  T    : R  : VAL DE TOUTES LES CMP
! IN  CO   : I  : TABLE DES NOMBRES DE COUCHES
! IN  SP   : I  : TABLE DES NOMBRES DE SOUS-PT
! IN  S    : R  : TABLE DES ABSCISSES CURVILIGNES
! IN  X    : R  : TABLES DES COORDONNEES
! IN  CMP  : K8 : NOMS DE TOUTES LES CMP
! IN  ND   : K8 : NOMS DES EVENTUELS NOEUDS
! IN  SDM  : K24: NOM OJB SD MAILLES
! IN  NBPOIN : I  : NOMBRE DE POINTS
! IN  DOCU : K4 : TYPE DE LIEU
! IN  NBCMP : I  : NOMBRE TOTAL DE CMP
!     ------------------------------------------------------------------
    integer(kind=8) :: nbpar, ilign, nbsp, i, ikk, l, jam, nbco, lc, is, ic, valei(1052), n1, adrval
    integer(kind=8) :: nbmail, j, adracc, jacc, ik, ir, ii, ivari(1000), nbcmp2, jvari, ico, lm, im
    integer(kind=8) :: nbvari, nbacc, nbpr, jaces, iac, iadr, iord(1)
    real(kind=8) :: prec, valer(1050)
    complex(kind=8) :: c16b
    aster_logical :: exist, erreur
    character(len=3) :: typpar
    character(len=7) :: kii
    character(len=8) :: k8b, acces, nomres, ctype, crit
    character(len=16) :: intitu
    character(len=24) :: nomval, nomacc, nnores, nopara(1053), nomjv
    character(len=80) :: valek(1051)
    character(len=8), pointer :: nom_para(:) => null()
    character(len=8), pointer :: typ_para(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    if (nbcmp .le. 0) goto 999
!
    if (docu .ne. 'LSTN') goto 999
!
    call jelira(jexnum(xnovar, iocc), 'LONUTI', nbvari)
    if (nbvari .ne. 0) then
        call jelira(jexnum(xnovar, iocc), 'LONUTI', nbvari)
        call jeveuo(jexnum(xnovar, iocc), 'L', jvari)
        if (nbvari .eq. 1 .and. zi(jvari) .eq. -1) then
            nbcmp2 = sp(1)
        else
            nbcmp2 = nbvari
        end if
    else
        nbcmp2 = nbcmp
        do i = 1, nbcmp2, 1
            ivari(i) = i
        end do
    end if
    if (nbcmp2 .gt. 1000) then
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
    ik = 1
    ii = 0
    ir = 0
    nbpar = 1
    valek(ik) = intitu
    nopara(nbpar) = 'INTITULE'
!
!
    if (zk8(jacc) .eq. 'DIRECT  ') then
        call jeveuo(jexnum(ncheff//'.LSCHEFF', 1), 'L', jacc)
        nbpar = nbpar+1
        nopara(nbpar) = 'CHAM_GD'
        ik = ik+1
        valek(ik) = zk24(jacc) (1:8)
    else
        call jeveuo(nnores, 'L', jacc)
        nomres = zk16(jacc) (1:8)
        nbpar = nbpar+1
        nopara(nbpar) = 'RESU'
        ik = ik+1
        valek(ik) = nomres
        nbpar = nbpar+1
        nopara(nbpar) = 'NOM_CHAM'
        ik = ik+1
        valek(ik) = zk16(jacc+1)
        call jeveuo(nomacc, 'L', adracc)
        call jeveuo(nomval, 'L', adrval)
        acces = zk8(adracc)
        if (acces(1:1) .eq. 'O') then
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_ORDRE'
            ii = ii+1
            valei(ii) = zi(adrval+i1-1)
            nomjv = '&&RVRCCM.NOMS_ACCES'
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
                    nopara(nbpar) = zk16(jaces-1+iac)
                    if (ctype(1:1) .eq. 'I') then
                        ii = ii+1
                        valei(ii) = zi(iadr)
                    else if (ctype(1:1) .eq. 'R') then
                        ir = ir+1
                        valer(ir) = zr(iadr)
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
            nopara(nbpar) = 'NUME_ORDRE'
            call rsorac(nomres, 'NUME_MODE', zi(adrval+i1-1), 0.d0, k8b, &
                        c16b, prec, crit, iord, 1, &
                        n1)
            ii = ii+1
            valei(ii) = iord(1)
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_MODE'
            ii = ii+1
            valei(ii) = zi(adrval+i1-1)
        else if (acces(1:1) .eq. 'I') then
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_ORDRE'
            call rsorac(nomres, 'INST', 0, zr(adrval+i1-1), k8b, &
                        c16b, prec, crit, iord, 1, &
                        n1)
            ii = ii+1
            valei(ii) = iord(1)
            nbpar = nbpar+1
            nopara(nbpar) = 'INST'
            ir = ir+1
            valer(ir) = zr(adrval+i1-1)
        else if (acces(1:1) .eq. 'F') then
            nbpar = nbpar+1
            nopara(nbpar) = 'FREQ'
            ir = ir+1
            valer(ir) = zr(adrval+i1-1)
        end if
    end if
    if (docu .eq. 'LSTN') then
        call tbexip(nomtab, 'NOEUD', exist, typpar)
        if (.not. exist) then
            call tbajpa(nomtab, 1, 'NOEUD', 'K8')
        end if
        nbpar = nbpar+1
        nopara(nbpar) = 'NOEUD'
    end if
    call tbexip(nomtab, 'MAILLE', exist, typpar)
    if (.not. exist) then
        call tbajpa(nomtab, 1, 'MAILLE', 'K8')
    end if
    nbpar = nbpar+1
    nopara(nbpar) = 'MAILLE'
    nbpar = nbpar+1
    nopara(nbpar) = 'ABSC_CURV'
    nbpar = nbpar+1
    nopara(nbpar) = 'COOR_X'
    nbpar = nbpar+1
    nopara(nbpar) = 'COOR_Y'
    nbpar = nbpar+1
    nopara(nbpar) = 'COOR_Z'
    ic = 0
    is = 0
    do i = 1, nbpoin, 1
        nbco = co(i)
        nbsp = sp(i)
        if (nbco .gt. 1 .and. ic .eq. 0) then
            call tbexip(nomtab, 'NUME_COUCHE', exist, typpar)
            if (.not. exist) then
                call tbajpa(nomtab, 1, 'NUME_COUCHE', 'I')
            end if
            ic = 1
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_COUCHE'
        end if
        if (nbsp .gt. 1 .and. is .eq. 0) then
            call tbexip(nomtab, 'NUME_GAUSS', exist, typpar)
            if (.not. exist) then
                call tbajpa(nomtab, 1, 'NUME_GAUSS', 'I')
            end if
            is = 1
            nbpar = nbpar+1
            nopara(nbpar) = 'NUME_GAUSS'
        end if
    end do
    if (nbvari .eq. 0) then
        do i = 1, nbcmp2, 1
            nbpar = nbpar+1
            nopara(nbpar) = cmp(i)
        end do
    else
        AS_ALLOCATE(vk8=nom_para, size=nbcmp2)
        AS_ALLOCATE(vk8=typ_para, size=nbcmp2)
        if (nbvari .eq. 1 .and. zi(jvari) .eq. -1) then
            do i = 1, nbcmp2, 1
                ivari(i) = i
                call codent(i, 'G', kii)
                nbpar = nbpar+1
                nopara(nbpar) = 'V'//kii
                nom_para(i) = 'V'//kii
                typ_para(i) = 'R'
            end do
        else
            do i = 1, nbcmp2, 1
                ivari(i) = zi(jvari+i-1)
                call codent(zi(jvari+i-1), 'G', kii)
                nbpar = nbpar+1
                nopara(nbpar) = 'V'//kii
                nom_para(i) = 'V'//kii
                typ_para(i) = 'R'
            end do
        end if
        call tbajpa(nomtab, nbcmp2, nom_para, typ_para)
    end if
!
    erreur = .false.
    if (nbpar .gt. 1050) erreur = .true.
    if (ii+2 .gt. 1050) erreur = .true.
    if (ir+4+nbcmp2 .gt. 1050) erreur = .true.
    if (ik .gt. 1050) erreur = .true.
    if (erreur) then
        call utmess('F', 'POSTRELE_12')
    end if
!
    ilign = 0
!
    do i = 1, nbpoin, 1
!
        valer(ir+1) = absc(i)
        valer(ir+2) = x(1+(i-1)*3)
        valer(ir+3) = x(2+(i-1)*3)
        valer(ir+4) = x(3+(i-1)*3)
!
        ikk = 0
        if (docu .eq. 'LSTN') then
            ikk = ikk+1
            valek(ik+ikk) = nd(i)
        end if
!
        if (docu .eq. 'LSTN') then
            call jeveuo(jexnum(sdm, i), 'L', jam)
            call jelira(jexnum(sdm, i), 'LONMAX', nbmail)
            l = 0
            nbco = co(i)
            nbsp = sp(i)
            j = padr(i)
            lm = nbcmp*nbsp
            lc = lm*nbmail
        else
            if (i .eq. 1) then
                nbsp = sp(i)
                nbco = co(i)
                j = 1
                lc = nbsp*nbcmp
                lm = lc*(2*nbco-1)
                lc = 2*lc
                nbmail = 1
                l = 1
            else if (i .eq. nbpoin) then
                j = i-1
                nbsp = sp(j)
                nbco = co(j)
                lc = nbsp*nbcmp
                lm = lc*(2*nbco-1)
                j = lc*(2*(i-2)*nbco+1)+1
                lc = 2*lc
                nbmail = 1
                l = 0
            else
                nbsp = sp(i)
                nbco = co(i)
                lc = nbsp*nbcmp
                lm = lc*(2*nbco-1)
                j = lc*(2*(i-2)*nbco+1)+1
                lc = 2*lc
                nbmail = 2
                l = 0
            end if
        end if
!        POUR UN CHAMP DE TYPE "VARI"
!        LES SOUS-POINTS SONT PRIS EN CHARGE PAR LE NBCMP
        if (nbvari .ne. 0) nbsp = 1
!
        do ico = 1, nbco, 1
!
            valei(ii+1) = ico
!
            do is = 1, nbsp, 1
!
                valei(ii+2) = is
!
                do im = 1, nbmail, 1
!
                    if (docu .eq. 'LSTN') then
                        valek(ik+ikk+1) = zk8(jam+l+im-1)
                    else
                        valek(ik+ikk+1) = '   -    '
                    end if
!
                    do ic = 1, nbcmp2, 1
                        valer(ir+4+ic) = t(j-1+(ico-1)*lc+(is-1)*nbcmp+(im-1)*lm+ivari(ic))
                    end do
!
                    call tbajli(nomtab, nbpar, nopara, valei, valer, &
                                [c16b], valek, ilign)
!
                end do
!
            end do
!
        end do
!
    end do
!
    if (nbvari .ne. 0) then
        AS_DEALLOCATE(vk8=nom_para)
        AS_DEALLOCATE(vk8=typ_para)
    end if
!
999 continue
    call jedema()
!
end subroutine
