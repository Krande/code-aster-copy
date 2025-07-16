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
! aslint: disable=W0413
! => real zero (init by calcul.F90)
subroutine fornpd(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/epseff.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "asterfort/trndgl.h"
#include "asterfort/trnflg.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "asterfort/vexpan.h"
#include "blas/daxpy.h"
    character(len=16) :: option, nomte
!
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE 3D
!     OPTION : FORC_NODA (REPRISE)
!
!
    integer(kind=8) :: icodre(26)
    character(len=10) :: phenom
!
    integer(kind=8) :: i, ib, icou, inte, intsn, intsr, j, k1, kpgs, kwgt, itab(7), iret
    integer(kind=8) :: icontm, jvDisp, imate, ivectu, jcara, jgeom, lzi, lzr
    integer(kind=8) :: nb1, nb2, npge, npgsr, npgsn, jnbspi, nbcou, nval, nbsp
!
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3), vecpt(9, 3, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: epais
    real(kind=8) :: rotfm(9)
    real(kind=8) :: deplm(42), effint(42), vecl(48), vecll(51)
    real(kind=8) :: sgmtd(5)
    real(kind=8) :: ksi3s2
    real(kind=8) :: sigtmp(5), ftemp(40), sigref
    real(kind=8) :: zero, zic, zmin, coef, hepa, hic
!
    character(len=16) :: kmess(2)
    blas_int :: b_incx, b_incy, b_n
!
    parameter(npge=3)
! DEB
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
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    if (nbcou .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PCACOQU', 'L', jcara)
    epais = zr(jcara)
    zmin = -epais/2.d0
    hic = epais/nbcou
!
    call jevech('PMATERC', 'L', imate)
!
    if (option .eq. 'FORC_NODA') then
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, &
                    itab=itab)
        icontm = itab(1)
        nbsp = itab(7)
        if (nbsp .ne. npge*nbcou) then
            call utmess('F', 'ELEMENTS_4')
        end if
    else if (option .eq. 'REFE_FORC_NODA') then
        call terefe('SIGM_REFE', 'MECA_COQUE3D', sigref)
    end if
!
    if (option .eq. "FORC_NODA") then
        call jevech('PDEPLAR', 'L', jvDisp)
    else
        call jevech('PDEPLMR', 'L', jvDisp)
    end if
!
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
!
    if (phenom .ne. 'ELAS' .and. phenom .ne. 'ELAS_ORTH') then
        call utmess('F', 'ELEMENTS_44', sk=phenom)
    end if
!
    call vectan(nb1, nb2, zr(jgeom), zr(lzr), vecta, &
                vectn, vectpt)
!
    call trndgl(nb2, vectn, vectpt, zr(jvDisp), deplm, &
                rotfm)
!
    do i = 1, 5*nb1+2
        effint(i) = 0.d0
    end do
!
    if (option .eq. 'REFE_FORC_NODA') then
        call r8inir(nb1*5, 0.d0, ftemp, 1)
    end if
!
    kwgt = 0
    kpgs = 0
    do icou = 1, nbcou
        do inte = 1, npge
            if (inte .eq. 1) then
                zic = zmin+(icou-1)*hic
                coef = 1.d0/3.d0
            else if (inte .eq. 2) then
                zic = zmin+hic/2.d0+(icou-1)*hic
                coef = 4.d0/3.d0
            else
                zic = zmin+hic+(icou-1)*hic
                coef = 1.d0/3.d0
            end if
            ksi3s2 = zic/hic
            hepa = hic
!
!   CALCUL DE BTDMR, BTDSR : M=MEMBRANE , S=CISAILLEMENT , R=REDUIT
!
            do intsr = 1, npgsr
                call mahsms(0, nb1, zr(jgeom), ksi3s2, intsr, &
                            zr(lzr), hepa, vectn, vectg, vectt, &
                            hsfm, hss)
!
                call hsj1ms(hepa, vectg, vectt, hsfm, hss, &
                            hsj1m, hsj1s)
!
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            hepa, vectpt, hsj1m, hsj1s, btdm, &
                            btds)
            end do
!
            do intsn = 1, npgsn
!
                call mahsf(1, nb1, zr(jgeom), ksi3s2, intsn, &
                           zr(lzr), hepa, vectn, vectg, vectt, &
                           hsf)
!
                call hsj1f(intsn, zr(lzr), hepa, vectg, vectt, &
                           hsf, kwgt, hsj1fx, wgt)
!
                wgt = coef*wgt
!
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), hepa, vectpt, hsj1fx, btdf)
!
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
!
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
!
                if (option .eq. 'FORC_NODA') then
                    sgmtd(1) = zr(icontm-1+k1+1)
                    sgmtd(2) = zr(icontm-1+k1+2)
                    sgmtd(3) = zr(icontm-1+k1+4)
                    sgmtd(4) = zr(icontm-1+k1+5)
                    sgmtd(5) = zr(icontm-1+k1+6)
!
                    call epseff('EFFORI', nb1, [0.d0], btild, sgmtd, &
                                [0.d0], wgt, effint)
!
                else if (option .eq. 'REFE_FORC_NODA') then
!
!      CALCUL DES FORCES NODALES DE REFERENCE
!      EN AFFECTANT LA VALEUR SIGM_REFE A CHAQUE CMP SUCCESSIVEMENT
!      POUR CHAQUE POINT D'INTEGRATION
!
                    call r8inir(5, 0.d0, sigtmp, 1)
!
                    do i = 1, 5
                        sigtmp(i) = sigref
                        call epseff('EFFORI', nb1, [0.d0], btild, sigtmp, &
                                    [0.d0], wgt, effint)
                        sigtmp(i) = 0.d0
                        do j = 1, nb1*5
                            ftemp(j) = ftemp(j)+abs(effint(j))
                        end do
                    end do
!
                end if
!
            end do
        end do
    end do
!
!      ON PREND LA VALEUR MOYENNE DES FORCES NODALES DE REFERENCE
!
    if (option .eq. 'REFE_FORC_NODA') then
        nval = nbcou*npge*npgsn*5
        b_n = to_blas_int(nb1*5)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/nval, ftemp, b_incx, effint, &
                   b_incy)
    end if
!
!
!-- EXPANSION DU CHAMP
    call vexpan(nb1, effint, vecl)
!
    do i = 1, 6*nb1
        vecll(i) = vecl(i)
    end do
    vecll(6*nb1+1) = effint(5*nb1+1)
    vecll(6*nb1+2) = effint(5*nb1+2)
!        VECLL(6*NB1+3)=0.D0
!
!     ICI PAS DE CONTRIBUTION DES DDL DE LA ROTATION FICTIVE DANS EFFINT
!
    zero = 0.d0
    do i = 1, nb1
        vecll(6*i) = zero*rotfm(i)
    end do
    i = nb2
    vecll(6*nb1+3) = zero*rotfm(nb2)
!
!     TRANFORMATION DANS REPERE GLOBAL PUIS STOCKAGE
!
    do ib = 1, nb2
        do i = 1, 2
            do j = 1, 3
                vecpt(ib, i, j) = vectpt(ib, i, j)
            end do
        end do
        vecpt(ib, 3, 1) = vectn(ib, 1)
        vecpt(ib, 3, 2) = vectn(ib, 2)
        vecpt(ib, 3, 3) = vectn(ib, 3)
    end do
!
    call jevech('PVECTUR', 'E', ivectu)
!
    call trnflg(nb2, vecpt, vecll, zr(ivectu))
!
    if (option .eq. 'REFE_FORC_NODA') then
        do j = 1, 51
            if (zr(ivectu+j-1) .eq. 0.) then
                kmess(1) = 'COQUE3D'
                kmess(2) = 'SIGM_REFE'
                call utmess('F', 'MECANONLINE5_59', nk=2, valk=kmess)
            end if
        end do
    end if
!
end subroutine
