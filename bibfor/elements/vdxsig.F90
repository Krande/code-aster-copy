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

subroutine vdxsig(nomte, option, xi, nb1, npgsr, &
                  sigmpg, effgt, nbcou)
    implicit none
#include "jeveux.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrth.h"
#include "asterfort/rcvarc.h"
#include "asterfort/trndgl.h"
#include "asterfort/utmess.h"
#include "asterfort/vdefge.h"
#include "asterfort/vdesga.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
    real(kind=8) :: sigmpg(*)
!
! CALCUL DES OPTIONS SIEF_ELGA
!                    EFGE_ELNO
!  POUR LES COQUE_3D
!
! ======================================================================
    character(len=16) :: nomte
    character(len=*) :: option
    integer(kind=8) :: nb1, nb2, npge, npgsr, npgsn
    integer(kind=8) :: nbcou, jcou, imoy, iret, iret1, iret2, iret3
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icou, indith, inte, intsn, intsr, j
    integer(kind=8) :: jcara, jdepg, k, k1, kpgs, kwgt, lzi
    integer(kind=8) :: lzr
    real(kind=8) :: tref
!-----------------------------------------------------------------------

    real(kind=8) :: xi(3, 9), sig(nbcou*162), eps(nbcou*162), tem(nbcou*27)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: depl(42), rotf(9)
    real(kind=8) :: epsiln(6, 27), sigma(6, 27), effgt(8, 9)
    real(kind=8) :: tempga(27)
    real(kind=8) :: young, nu, alpha, epais
    real(kind=8) :: xi3, p1xi3, p2xi3, p3xi3, tinf, tmoy, tsup
    real(kind=8) :: epsval(3), ksi3s2, zero, deux, hic, zmin, un
!
!     LES TABLEAUX SIG,EPS,TEM ONT ETE ALLOUES DE FACON STATIQUE POUR
!     OPTIMISER LE CPU CAR LES APPELS A WKVECT DANS LES TE SONT COUTEUX.
    zero = 0.0d0
    un = 1.0d0
    deux = 2.0d0
! NOMBRE DE POINTS DE GAUSS DANS LA TRANCHE
! (POUR RESTER COHERENT AVEC SIEF_ELGA EN PLASTICITE )
    npge = 3
    epsval(1) = -1.d0
    epsval(2) = 0.d0
    epsval(3) = 1.d0
!
!     RECUPERATION DES OBJETS
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
    call jevech('PCACOQU', 'L', jcara)
    epais = zr(jcara)
!
    call jevech('PNBSP_I', 'L', jcou)
    nbcou = zi(jcou)
    imoy = (3*nbcou+1)/2
    hic = un/nbcou
    zmin = -0.5d0
!
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tref, iret)
    call vectan(nb1, nb2, xi, zr(lzr), vecta, &
                vectn, vectpt)
!
    call jevech('PDEPLAR', 'L', jdepg)
!
    call trndgl(nb2, vectn, vectpt, zr(jdepg), depl, &
                rotf)
!
    kwgt = 0
    kpgs = 0
!
    do icou = 1, nbcou
!
        do inte = 1, npge
!
!     CALCUL DE BTDMR, BTDSR : M=MEMBRANE , S=CISAILLEMENT , R=REDUIT
!
!       -- COTE DES POINTS D'INTEGRATION
!       --------------------------------
            if (inte .eq. 1) then
                ksi3s2 = zmin+(icou-1)*hic
            else if (inte .eq. 2) then
                ksi3s2 = zmin+hic/deux+(icou-1)*hic
            else
                ksi3s2 = zmin+hic+(icou-1)*hic
            end if
!
            do intsr = 1, npgsr
                call mahsms(0, nb1, xi, ksi3s2, intsr, &
                            zr(lzr), epais, vectn, vectg, vectt, &
                            hsfm, hss)
!
                call hsj1ms(epais, vectg, vectt, hsfm, hss, &
                            hsj1m, hsj1s)
!
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            epais, vectpt, hsj1m, hsj1s, btdm, &
                            btds)
            end do
!
            do intsn = 1, npgsn
!
!     CALCUL DE BTDFN : F=FLEXION , N=NORMAL
!     ET DEFINITION DE WGT=PRODUIT DES POIDS ASSOCIES AUX PTS DE GAUSS
!                          (NORMAL) ET DU DETERMINANT DU JACOBIEN
!
                call mahsf(1, nb1, xi, ksi3s2, intsn, &
                           zr(lzr), epais, vectn, vectg, vectt, &
                           hsf)
!
                call hsj1f(intsn, zr(lzr), epais, vectg, vectt, &
                           hsf, kwgt, hsj1fx, wgt)
!
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), epais, vectpt, hsj1fx, btdf)
!
!     CALCUL DE BTDMN, BTDSN : M=MEMBRANE , S=CISAILLEMENT , N=NORMAL
!     FORMATION DE BTILD
!
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
!
!     APPEL DE MATRTH POUR RECUPERER INDITH AFIN DE SAVOIR SI
!     ALPHA EST DONNE C'EST A DIRE SI THERMIQUE
!
                call rcvarc(' ', 'TEMP', '+', 'MASS', intsn, &
                            3*icou-2, tinf, iret1)
                call rcvarc(' ', 'TEMP', '+', 'MASS', intsn, &
                            3*icou-1, tmoy, iret2)
                call rcvarc(' ', 'TEMP', '+', 'MASS', intsn, &
                            3*icou, tsup, iret3)
                if ((iret1+iret2+iret3) .eq. 0) then
                    call matrth('MASS', npgsn, young, nu, alpha, &
                                indith)
                    xi3 = epsval(inte)
                    p1xi3 = 1-xi3*xi3
                    p2xi3 = -xi3*(1-xi3)/2.d0
                    p3xi3 = xi3*(1+xi3)/2.d0
                    if (iret .eq. 1) then
                        call utmess('F', 'CALCULEL_15')
                    else
                        tem(kwgt) = tmoy*p1xi3+tinf*p2xi3+tsup*p3xi3
                        tem(kwgt) = tem(kwgt)-tref
                    end if
                else
                    indith = -1
                end if
!
                call vdesga(kwgt, nb1, nb2, depl, btild, &
                            indith, alpha, tem, eps, sig, &
                            vectt)
!
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbcou+npge*(icou-1)+inte-1)
                if (option .eq. 'EPSI_ELGA') then
                    do i = 1, 6
                        sigmpg(k1+i) = eps(i+6*(kpgs-1))
                    end do
                else if (option .eq. 'SIEF_ELGA') then
                    do i = 1, 6
                        sigmpg(k1+i) = sig(i+6*(kpgs-1))
!                  write (6,*) "sigmpg(k1+i)",sigmpg(k1+i)
                    end do
                end if
!
            end do
        end do
    end do
!
    kwgt = 0
    kpgs = 0
!
    do inte = 1, npge
!
!     CALCUL DE BTDMR, BTDSR : M=MEMBRANE , S=CISAILLEMENT , R=REDUIT
!
        ksi3s2 = epsval(inte)/2.d0
        do intsr = 1, npgsr
            call mahsms(0, nb1, xi, ksi3s2, intsr, &
                        zr(lzr), epais, vectn, vectg, vectt, &
                        hsfm, hss)
!
            call hsj1ms(epais, vectg, vectt, hsfm, hss, &
                        hsj1m, hsj1s)
!
            call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                        epais, vectpt, hsj1m, hsj1s, btdm, &
                        btds)
!
!       CALL BTDMSP(NB1,NB2,XI,INTE,INTSR,ZR(LZR),EPAIS,VECTPT,
!    &                                       HSJ1M,HSJ1S,BTDM,BTDS)
            call mahsf(0, nb1, xi, ksi3s2, intsr, &
                       zr(lzr), epais, vectn, vectg, vectt, &
                       hsf)
!
            call hsj1f(intsr, zr(lzr), epais, vectg, vectt, &
                       hsf, kwgt, hsj1fx, wgt)
!
            call btdfn(0, nb1, nb2, ksi3s2, intsr, &
                       zr(lzr), epais, vectpt, hsj1fx, btdf)
!     CALL BTDFP(0,NB1,NB2,XI,INTE,INTSR,ZR(LZR),EPAIS,VECTPT,HSJ1FX,
!    &                                                             BTDF)
!
            call btdmsn(0, nb1, intsr, npgsr, zr(lzr), &
                        btdm, btdf, btds, btild)
!
!     CALL BTILDP(0,NB1,XI,INTE,INTSR,NPGSR,ZR(LZR),BTDM,BTDF,BTDS,
!    &                                                            BTILD)
!
!     APPEL DE MATRTH POUR RECUPERER INDITH AFIN DE SAVOIR SI
!     ALPHA EST DONNE C'EST A DIRE SI THERMIQUE
!
            call rcvarc(' ', 'TEMP', '+', 'RIGI', intsr, &
                        1, tinf, iret1)
            call rcvarc(' ', 'TEMP', '+', 'RIGI', intsr, &
                        imoy, tmoy, iret2)
            call rcvarc(' ', 'TEMP', '+', 'RIGI', intsr, &
                        3*nbcou, tsup, iret3)
            if ((iret1+iret2+iret3) .eq. 0) then
                call matrth('RIGI', npgsr, young, nu, alpha, &
                            indith)
                xi3 = epsval(inte)
                p1xi3 = 1-xi3*xi3
                p2xi3 = -xi3*(1-xi3)/2.d0
                p3xi3 = xi3*(1+xi3)/2.d0
                if (iret .eq. 1) then
                    call utmess('F', 'CALCULEL_15')
                else
                    tempga(kwgt) = tmoy*p1xi3+tinf*p2xi3+tsup*p3xi3
                    tempga(kwgt) = tempga(kwgt)-tref
                end if
            else
                tempga(kwgt) = 0.d0
            end if
!
            call vdesga(kwgt, nb1, nb2, depl, btild, &
                        indith, alpha, tempga, epsiln, sigma, &
                        vectt)
!
        end do
    end do
!
    if (option(1:9) .eq. 'EFGE_ELNO') then
!
        call vdefge(nomte, nb1, npgsr, zr(lzr), epais, &
                    sigma, effgt)
!
    end if
!
! --- DETERMINATION DES REPERES  LOCAUX DE L'ELEMENT AUX POINTS
! --- D'INTEGRATION ET STOCKAGE DE CES REPERES DANS LE VECTEUR .DESR :
!     --------------------------------------------------------------
    k = 0
    do intsr = 1, npgsr
        call vectgt(0, nb1, xi, zero, intsr, &
                    zr(lzr), epais, vectn, vectg, vectt)
!
        do j = 1, 3
            do i = 1, 3
                k = k+1
                zr(lzr+2000+k-1) = vectt(i, j)
            end do
        end do
    end do
!
end subroutine
