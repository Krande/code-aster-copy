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
subroutine rcevod(csigm, cinst, cnoc, sm, lfatig, &
                  lpmpb, lsn, csno, csne, flexio, &
                  csneo, csnee, cfao, cfae, cspo, &
                  cspe, cresu, kinti, it, jt, &
                  lrocht, symax, cpres, kemixt, cspto, &
                  cspte, cspmo, cspme, lsymm, csiex)
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rcevfu.h"
#include "asterfort/rcmcrt.h"
#include "asterfort/rctres.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
    integer(kind=8) :: it, jt
    real(kind=8) :: sm, symax
    aster_logical :: lfatig, lpmpb, lsn, flexio, lrocht, kemixt, lsymm
    character(len=16) :: kinti
    character(len=24) :: csigm, cinst, cnoc, csno, csne, csneo, csnee, cfao
    character(len=24) :: cfae, cspo, cspe, cresu, cpres, cspto, cspte, cspmo
    character(len=24) :: cspme, csiex
!     OPERATEUR POST_RCCM, TYPE_RESU_MECA='EVOLUTION'
!     TYPE_RESU = 'DETAILS'
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ncmp, jsigm, jinst, nbinst, jsno, jsne, ind, i1, i2, icmp, l1, l2
    integer(kind=8) :: l3, l4, l5, npara, ik, ir, i, vaio(5), vaie(5), npar1, jresp
    integer(kind=8) :: jsnee, jspo, jspe, jfao, jfae, jnoc, jresu, jspto, jspte, jspmo
    integer(kind=8) :: jspme, jsneo, jsigtot
    parameter(ncmp=6)
    real(kind=8) :: tpm(ncmp), tpb(ncmp), tpmpbo(ncmp), tpmpbe(ncmp), dco, dce, qlin
    real(kind=8) :: tresca, valo(39), vale(39), stlin, stpar, tfpto(ncmp), tfpte(ncmp)
    complex(kind=8) :: c16b
    character(len=8) :: nomres, typara(39)
    character(len=16) :: nomcmd, concep, nopara(39), vako(5), vake(5)
!
    integer(kind=8) :: nparen, nparpm, nparsn, nparse, nparf1, nparf2, nparf3, nparrt
    integer(kind=8) :: ifm, niv
    parameter(nparen=4, nparpm=6, nparsn=5, nparse=1,&
     &             nparf1=14, nparf2=13, nparrt=6, nparf3=17)
    character(len=8) :: typaen(nparen), typapm(nparpm), typasn(nparsn)
    character(len=8) :: typase(nparse), typaf1(nparf1), typaf2(nparf2)
    character(len=8) :: typart(nparrt), typaf3(nparf3)
    character(len=16) :: nopaen(nparen), nopapm(nparpm), nopasn(nparsn)
    character(len=16) :: nopase(nparse), nopaf1(nparf1), nopaf2(nparf2)
    character(len=16) :: nopart(nparrt), nopaf3(nparf3)
!
    data nopaen/'INTITULE', 'LIEU', 'SM', '3SM'/
    data typaen/'K16', 'K8', 'R', 'R'/
    data nopart/'TABL_PRES', 'SY', 'INST', 'SIGM_M_PRES',&
     &              'VALE_MAXI_LINE', 'VALE_MAXI_PARAB'/
    data typart/'K8', 'R', 'R', 'R', 'R', 'R'/
    data nopapm/'TABL_RESU', 'INST', 'PM', 'PB', 'PMB', 'F'/
    data typapm/'K8', 'R', 'R', 'R', 'R', 'R'/
    data nopasn/'TABL_RESU_1', 'INST_1',&
     &              'TABL_RESU_2', 'INST_2', 'SN'/
    data typasn/'K8', 'R', 'K8', 'R', 'R'/
    data nopase/'SN*'/
    data typase/'R'/
    data nopaf1/'TABL_RESU_1', 'INST_1', 'NB_OCCUR_1',&
     &              'TABL_RESU_2', 'INST_2', 'NB_OCCUR_2',&
     &              'SN', 'SN*', 'SP', 'KE', 'SALT', 'NADM',&
     &              'DOMMAGE', 'DOMMAGE_CUMU'/
    data typaf1/'K8', 'R', 'I', 'K8', 'R', 'I',&
     &              'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R'/
    data nopaf2/'TABL_RESU_1', 'INST_1', 'NB_OCCUR_1',&
     &              'TABL_RESU_2', 'INST_2', 'NB_OCCUR_2',&
     &              'SN', 'SP', 'KE', 'SALT', 'NADM',&
     &              'DOMMAGE', 'DOMMAGE_CUMU'/
    data typaf2/'K8', 'R', 'I', 'K8', 'R', 'I',&
     &              'R', 'R', 'R', 'R', 'R', 'R', 'R'/
    data nopaf3/'TABL_RESU_1', 'INST_1', 'NB_OCCUR_1',&
     &            'TABL_RESU_2', 'INST_2', 'NB_OCCUR_2',&
     &            'SN', 'SN*', 'SP', 'SP_MECA', 'SP_THER', 'KE_MECA',&
     &            'KE_THER', 'SALT', 'NADM', 'DOMMAGE', 'DOMMAGE_CUMU'/
    data typaf3/'K8', 'R', 'I', 'K8', 'R', 'I', 'R', 'R',&
     &              'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R'/
! DEB ------------------------------------------------------------------
    call jemarq()
!
    c16b = (0.d0, 0.d0)
    call getres(nomres, concep, nomcmd)
    call infniv(ifm, niv)
!
    call jeveuo(cresu, 'L', jresu)
    call jeveuo(cinst, 'L', jinst)
    call jelira(cinst, 'LONMAX', nbinst)
!
! --- CREATION DE LA TABLE
!
    if (it .eq. 1 .and. jt .eq. 1) then
        npara = nparen
        do i = 1, nparen
            nopara(i) = nopaen(i)
            typara(i) = typaen(i)
        end do
        if (lrocht) then
            do i = 1, nparrt
                nopara(npara+i) = nopart(i)
                typara(npara+i) = typart(i)
            end do
            npara = npara+nparrt
        end if
        if (lpmpb) then
            do i = 1, nparpm
                nopara(npara+i) = nopapm(i)
                typara(npara+i) = typapm(i)
            end do
            npara = npara+nparpm
        end if
        if (lfatig) then
            if (kemixt) then
                do i = 1, nparf3
                    nopara(npara+i) = nopaf3(i)
                    typara(npara+i) = typaf3(i)
                end do
                npara = npara+nparf3
            else
                if (flexio) then
                    do i = 1, nparf1
                        nopara(npara+i) = nopaf1(i)
                        typara(npara+i) = typaf1(i)
                    end do
                    npara = npara+nparf1
                else
                    do i = 1, nparf2
                        nopara(npara+i) = nopaf2(i)
                        typara(npara+i) = typaf2(i)
                    end do
                    npara = npara+nparf2
                end if
            end if
        else
            if (lsn) then
                do i = 1, nparsn
                    nopara(npara+i) = nopasn(i)
                    typara(npara+i) = typasn(i)
                end do
                npara = npara+nparsn
                if (flexio) then
                    do i = 1, nparse
                        nopara(npara+i) = nopase(i)
                        typara(npara+i) = typase(i)
                    end do
                    npara = npara+nparse
                end if
            end if
        end if
!
        call tbcrsd(nomres, 'G')
        call tbajpa(nomres, npara, nopara, typara)
    end if
!
! --- LES LIGNES  LIEU ET SM
!
    ik = 0
    npara = nparen
    do i = 1, nparen
        nopara(i) = nopaen(i)
    end do
    ik = ik+1
    vako(ik) = kinti
    vake(ik) = kinti
    ik = ik+1
    vako(ik) = 'ORIG'
    vake(ik) = 'EXTR'
!
    valo(1) = sm
    vale(1) = sm
!
    valo(2) = 3*sm
    vale(2) = 3*sm
!
! --- POUR LE ROCHET THERMIQUE
!
    if (lrocht) then
        call jeveuo(csigm, 'L', jsigm)
        call jeveuo(cpres, 'L', jresp)
        do i = 1, nparrt
            nopara(npara+i) = nopart(i)
        end do
        npar1 = npara+nparrt
        vako(ik+1) = zk8(jresp-1+jt)
        vake(ik+1) = zk8(jresp-1+jt)
        ir = 2+1
        valo(ir) = symax
        vale(ir) = symax
        do i = 1, nbinst
            ir = 3+1
            valo(ir) = zr(jinst+i-1)
            vale(ir) = zr(jinst+i-1)
            do icmp = 1, ncmp
                if (lsymm) then
                    l3 = 6*ncmp*nbinst+ncmp*(i-1)+icmp
                    tpm(icmp) = zr(jsigm-1+l3)
                else
                    l3 = 4*ncmp*nbinst+ncmp*(i-1)+icmp
                    tpm(icmp) = zr(jsigm-1+l3)
                end if
            end do
            call rctres(tpm, tresca)
            call rcmcrt(symax, tresca, stlin, stpar)
!
            ir = ir+1
            valo(ir) = tresca
            vale(ir) = tresca
            ir = ir+1
            valo(ir) = stlin
            vale(ir) = stlin
            ir = ir+1
            valo(ir) = stpar
            vale(ir) = stpar
            call tbajli(nomres, npar1, nopara, vaio, valo, &
                        [c16b], vako, 0)
            call tbajli(nomres, npar1, nopara, vaie, vale, &
                        [c16b], vake, 0)
        end do
!
    end if
!
! --- POUR L'OPTION "PMPB"
!
    if (lpmpb) then
!
! --- LES CRITERES DE NIVEAU 0 VISENT A PREMUNIR LE MATERIEL CONTRE LES
!     DOMMAGES DE DEFORMATION EXCESSIVE, D'INSTABILITE PLASTIQUE ET
!     D'INSTABILITE ELASTIQUE ET ELASTOPLASTIQUE.
!     ON NE PREND QUE LA PARTIE MECANIQUE
!
        do i = 1, nparpm
            nopara(npara+i) = nopapm(i)
        end do
        npar1 = npara+nparpm
!
        call jeveuo(csigm, 'L', jsigm)
        call jeveuo(csiex, 'L', jsigtot)
        do i = 1, nbinst
            vako(ik+1) = zk8(jresu+i-1)
            ir = 2+1
            valo(ir) = zr(jinst+i-1)
            if (lsymm) then
                do icmp = 1, ncmp
                    l1 = ncmp*(i-1)+icmp
                    l2 = ncmp*nbinst+ncmp*(i-1)+icmp
                    l4 = 3*ncmp*nbinst+ncmp*(i-1)+icmp
                    l5 = 4*ncmp*nbinst+ncmp*(i-1)+icmp
                    tpm(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l4)
                    tpb(icmp) = zr(jsigm-1+l2)-zr(jsigm-1+l5)
                    tpmpbo(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l4)+(zr(jsigm-1+l2)-zr(jsigm-1&
                                   &+l5))
                    qlin = zr(jsigm-1+l4)+zr(jsigm-1+l5)
                    tfpto(icmp) = zr(jsigtot-1+l1)-tpmpbo(icmp)-qlin
                end do
            else
                do icmp = 1, ncmp
                    l1 = ncmp*(i-1)+icmp
                    l2 = ncmp*nbinst+ncmp*(i-1)+icmp
                    l3 = 2*ncmp*nbinst+ncmp*(i-1)+icmp
                    l4 = 3*ncmp*nbinst+ncmp*(i-1)+icmp
                    tpm(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l3)
                    tpb(icmp) = zr(jsigm-1+l2)-zr(jsigm-1+l4)
                    tpmpbo(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l2)-(zr(jsigm-1+l3)-zr(jsigm-1&
                                   &+l4))
                    qlin = zr(jsigm-1+l3)+zr(jsigm-1+l4)
                    tfpto(icmp) = zr(jsigtot-1+l1)-tpmpbo(icmp)-qlin
                end do
                print *, ''
            end if
            call rctres(tpm, tresca)
            ir = ir+1
            valo(ir) = tresca
            call rctres(tpb, tresca)
            ir = ir+1
            valo(ir) = tresca
            call rctres(tpmpbo, tresca)
            ir = ir+1
            valo(ir) = tresca
            call rctres(tfpto, tresca)
            ir = ir+1
            valo(ir) = tresca
            call tbajli(nomres, npar1, nopara, vaio, valo, &
                        [c16b], vako, 0)
        end do
        do i = 1, nbinst
            vake(ik+1) = zk8(jresu+i-1)
            ir = 2+1
            vale(ir) = zr(jinst+i-1)
            if (lsymm) then
                do icmp = 1, ncmp
                    l1 = ncmp*(i-1)+icmp
                    l2 = 2*ncmp*nbinst+ncmp*(i-1)+icmp
                    l4 = 3*ncmp*nbinst+ncmp*(i-1)+icmp
                    l5 = 5*ncmp*nbinst+ncmp*(i-1)+icmp
                    tpm(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l4)
                    tpb(icmp) = zr(jsigm-1+l2)-zr(jsigm-1+l5)
                    tpmpbe(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l4)+(zr(jsigm-1+l2)-zr(jsigm-1&
                                   &+l5))
                    qlin = zr(jsigm-1+l4)+zr(jsigm-1+l5)
                    tfpte(icmp) = zr(jsigtot-1+ncmp*nbinst+ncmp*(i-1)+icmp)-tpmpbe(icmp)-qlin
                end do
            else
                do icmp = 1, ncmp
                    l1 = ncmp*(i-1)+icmp
                    l2 = ncmp*nbinst+ncmp*(i-1)+icmp
                    l3 = 2*ncmp*nbinst+ncmp*(i-1)+icmp
                    l4 = 3*ncmp*nbinst+ncmp*(i-1)+icmp
                    tpm(icmp) = zr(jsigm-1+l1)-zr(jsigm-1+l3)
                    tpb(icmp) = zr(jsigm-1+l2)-zr(jsigm-1+l4)
                    tpmpbe(icmp) = zr(jsigm-1+l1)+zr(jsigm-1+l2)-(zr(jsigm-1+l3)+zr(jsigm-1&
                                   &+l4))
                    qlin = zr(jsigm-1+l3)+zr(jsigm-1+l4)
                    tfpte(icmp) = zr(jsigtot-1+l2)-tpmpbe(icmp)-qlin
                end do
            end if
            call rctres(tpm, tresca)
            ir = ir+1
            vale(ir) = tresca
            call rctres(tpb, tresca)
            ir = ir+1
            vale(ir) = tresca
            call rctres(tpmpbe, tresca)
            ir = ir+1
            vale(ir) = tresca
            call rctres(tfpte, tresca)
            ir = ir+1
            vale(ir) = tresca
            call tbajli(nomres, npar1, nopara, vaie, vale, &
                        [c16b], vake, 0)
        end do
    end if
!
! --- POUR L'OPTION "SN"
!
    if (lsn .and. .not. lfatig) then
        do i = 1, nparsn
            nopara(npara+i) = nopasn(i)
        end do
        npar1 = npara+nparsn
        if (flexio) then
            npar1 = npar1+1
            nopara(npar1) = 'SN*'
        end if
!
        call jeveuo(csno, 'L', jsno)
        if (flexio) call jeveuo(csneo, 'L', jsneo)
!
        ind = 0
        do i1 = 1, nbinst
            ind = ind+1
            do i2 = i1+1, nbinst
                ind = ind+1
                vako(ik+1) = zk8(jresu+i1-1)
                vako(ik+2) = zk8(jresu+i2-1)
                ir = 2+1
                valo(ir) = zr(jinst+i1-1)
                ir = ir+1
                valo(ir) = zr(jinst+i2-1)
                ir = ir+1
                valo(ir) = zr(jsno+ind-1)
                if (flexio) then
                    ir = ir+1
                    valo(ir) = zr(jsneo+ind-1)
                end if
                call tbajli(nomres, npar1, nopara, vaio, valo, &
                            [c16b], vako, 0)
            end do
        end do
!
        call jeveuo(csne, 'L', jsne)
        if (flexio) call jeveuo(csnee, 'L', jsnee)
        ind = 0
        do i1 = 1, nbinst
            ind = ind+1
            do i2 = i1+1, nbinst
                ind = ind+1
                vake(ik+1) = zk8(jresu+i1-1)
                vake(ik+2) = zk8(jresu+i2-1)
                ir = 2+1
                vale(ir) = zr(jinst+i1-1)
                ir = ir+1
                vale(ir) = zr(jinst+i2-1)
                ir = ir+1
                vale(ir) = zr(jsne+ind-1)
                if (flexio) then
                    ir = ir+1
                    vale(ir) = zr(jsnee+ind-1)
                end if
                call tbajli(nomres, npar1, nopara, vaie, vale, &
                            [c16b], vake, 0)
            end do
        end do
    end if
!
! --- POUR L'OPTION "FATIGUE"
!
    if (lfatig) then
        if (kemixt) then
            do i = 1, nparf3
                nopara(npara+i) = nopaf3(i)
            end do
            npar1 = npara+nparf3-1
        else if (flexio) then
            do i = 1, nparf1
                nopara(npara+i) = nopaf1(i)
            end do
            npar1 = npara+nparf1-1
        else
            do i = 1, nparf2
                nopara(npara+i) = nopaf2(i)
            end do
            npar1 = npara+nparf2-1
        end if
!
        call jeveuo(cnoc, 'L', jnoc)
        call jeveuo(csno, 'L', jsno)
        if (flexio) call jeveuo(csneo, 'L', jsneo)
        call jeveuo(cspo, 'L', jspo)
        call jeveuo(cfao, 'L', jfao)
        if (kemixt) then
            call jeveuo(cspmo, 'L', jspmo)
            call jeveuo(cspto, 'L', jspto)
        end if
        ind = 0
        do i1 = 1, nbinst
            ind = ind+1
            do i2 = i1+1, nbinst
                ind = ind+1
                vako(ik+1) = zk8(jresu+i1-1)
                vako(ik+2) = zk8(jresu+i2-1)
                vaio(1) = zi(jnoc+i1-1)
                vaio(2) = zi(jnoc+i2-1)
                ir = 2+1
                valo(ir) = zr(jinst+i1-1)
                ir = ir+1
                valo(ir) = zr(jinst+i2-1)
                ir = ir+1
                valo(ir) = zr(jsno+ind-1)
                if (flexio) then
                    ir = ir+1
                    valo(ir) = zr(jsneo+ind-1)
                end if
                ir = ir+1
                valo(ir) = zr(jspo-1+ind)
                if (kemixt) then
                    ir = ir+1
                    valo(ir) = zr(jspmo-1+ind)
                    ir = ir+1
                    valo(ir) = zr(jspto-1+ind)
                end if
                ir = ir+1
                valo(ir) = zr(jfao-1+5*(ind-1)+1)
                if (kemixt) then
                    ir = ir+1
                    valo(ir) = zr(jfao-1+5*(ind-1)+5)
                end if
                ir = ir+1
                valo(ir) = zr(jfao-1+5*(ind-1)+2)
                ir = ir+1
                valo(ir) = zr(jfao-1+5*(ind-1)+3)
                ir = ir+1
                valo(ir) = zr(jfao-1+5*(ind-1)+4)
                call tbajli(nomres, npar1, nopara, vaio, valo, &
                            [c16b], vako, 0)
            end do
        end do
!
        call jeveuo(csne, 'L', jsne)
        if (flexio) call jeveuo(csnee, 'L', jsnee)
        call jeveuo(cspe, 'L', jspe)
        call jeveuo(cfae, 'L', jfae)
        if (kemixt) then
            call jeveuo(cspme, 'L', jspme)
            call jeveuo(cspte, 'L', jspte)
        end if
        ind = 0
        do i1 = 1, nbinst
            ind = ind+1
            do i2 = i1+1, nbinst
                ind = ind+1
                vake(ik+1) = zk8(jresu+i1-1)
                vake(ik+2) = zk8(jresu+i2-1)
                vaie(1) = zi(jnoc+i1-1)
                vaie(2) = zi(jnoc+i2-1)
                ir = 2+1
                vale(ir) = zr(jinst+i1-1)
                ir = ir+1
                vale(ir) = zr(jinst+i2-1)
                ir = ir+1
                vale(ir) = zr(jsne+ind-1)
                if (flexio) then
                    ir = ir+1
                    vale(ir) = zr(jsnee+ind-1)
                end if
                ir = ir+1
                vale(ir) = zr(jspe-1+ind)
                if (kemixt) then
                    ir = ir+1
                    vale(ir) = zr(jspme-1+ind)
                    ir = ir+1
                    vale(ir) = zr(jspte-1+ind)
                end if
                ir = ir+1
                vale(ir) = zr(jfae-1+5*(ind-1)+1)
                if (kemixt) then
                    ir = ir+1
                    vale(ir) = zr(jfae-1+5*(ind-1)+5)
                end if
                ir = ir+1
                vale(ir) = zr(jfae-1+5*(ind-1)+2)
                ir = ir+1
                vale(ir) = zr(jfae-1+5*(ind-1)+3)
                ir = ir+1
                vale(ir) = zr(jfae-1+5*(ind-1)+4)
                call tbajli(nomres, npar1, nopara, vaie, vale, &
                            [c16b], vake, 0)
            end do
        end do
!
        nopara(npara+1) = 'DOMMAGE_CUMU'
        npar1 = npara+1
        if (niv .eq. 2) write (6, *) '******* ORIGINE DU SEGMENT *******'
        call rcevfu(cnoc, cfao, dco)
        if (niv .eq. 2) write (6, *) '******* EXTREMITE DU SEGMENT *******'
        call rcevfu(cnoc, cfae, dce)
        valo(3) = dco
        vale(3) = dce
        call tbajli(nomres, npar1, nopara, vaio, valo, &
                    [c16b], vako, 0)
        call tbajli(nomres, npar1, nopara, vaie, vale, &
                    [c16b], vake, 0)
!
    end if
!
    call jedema()
end subroutine
