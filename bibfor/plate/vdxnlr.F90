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
subroutine vdxnlr(option, nomte, xi, rig, nb1, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btkb.h"
#include "asterfort/epseff.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrc2.h"
#include "asterfort/matrkb.h"
#include "asterfort/moytpg.h"
#include "asterfort/nmcomp.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/trndgl.h"
#include "asterfort/trnflg.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "asterfort/vexpan.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jnbspi
    character(len=32) :: elasKeyword
    character(len=16) :: option, nomte
    integer(kind=8) :: nb1, nb2, nddle, npge, npgsr, npgsn, itab(8), codret
    integer(kind=8) :: cod, ksp
    real(kind=8) :: xi(3, 9)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3), vecpt(9, 3, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42), wmatcb(5, 42), ktildi(42, 42)
    real(kind=8) :: ktild(42, 42), rig(51, 51)
    real(kind=8) :: ctor, epais, kappa
    integer(kind=8), parameter :: nbv = 2
    character(len=16), parameter :: nomres(nbv) = (/'E ', 'NU'/)
    integer(kind=8) :: valret(nbv)
    real(kind=8) :: valres(nbv)
    real(kind=8) :: rotfcm(9), rotfcp(9)
    real(kind=8) :: deplm(42), deplp(42)
    real(kind=8) :: epsi(5), depsi(5), eps2d(4), deps2d(4)
    real(kind=8) :: dtild(5, 5), sgmtd(5), effint(42), vecl(48), vecll(51)
    real(kind=8) :: sign(4), sigma(4), dsidep(6, 6), angmas(3)
    real(kind=8) :: matc(5, 5), valpar
    integer(kind=8) :: i, ib, icarcr, icontm, icontp, icou
    integer(kind=8) :: ideplm, ideplp, iinstm, iinstp, imate, inte, intsn
    integer(kind=8) :: intsr, iret, ivarim, ivarip, ivarix, ivectu, j
    integer(kind=8) :: jcara, jcrf, k1, k2, kpgs, kwgt, lgpg
    integer(kind=8) :: lzi, lzr, nbcou, nbvari, nddlet, ndimv
    real(kind=8) :: coef, crf, gxz, gyz, hic
    real(kind=8) :: x(1), zic, zmin
    parameter(npge=3)
    real(kind=8) :: ksi3s2
    aster_logical :: lVect, lMatr, lVari, lSigm

    blas_int :: b_incx, b_incy, b_n
    real(kind=8) :: cisail
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8), parameter :: ndimLdc = 2
    character(len=8), parameter :: typmod(2) = (/"C_PLAN  ", "        "/)
    type(Behaviour_Integ) :: BEHinteg
    character(len=4), parameter :: fami = "MASS"
    character(len=16), pointer :: compor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
!
    nddle = 5*nb1+2
    ktild = 0.d0
    effint = 0.d0
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get input fields
    call jevech('PMATERC', 'L', imate)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PNBSP_I', 'L', jnbspi)
    call jevech('PCARCRI', 'L', icarcr)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMP', 'L', ivarix)
    call jevech('PCACOQU', 'L', jcara)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=itab)
    if (itab(6) .le. 1) then
        lgpg = itab(7)
    else
        lgpg = itab(6)*itab(7)
    end if
    nbcou = zi(jnbspi-1+1)
    if (nbcou .le. 0) then
        call utmess('F', 'PLATE1_10')
    end if

! - Don"t use AFFE_CARA_ELEM/MASSIF
    angmas = r8vide()

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndimLdc, typmod, option, &
                              compor, zr(icarcr), &
                              zr(iinstm), zr(iinstp), &
                              fami, zi(imate), &
                              BEHinteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Properties of behaviour
    read (compor(NVAR), '(I16)') nbvari
!
    epais = zr(jcara)
    kappa = zr(jcara+3)
    ctor = zr(jcara+4)
    zmin = -epais/2.d0
    hic = epais/nbcou

! - Get output fields
    if (option .eq. 'RAPH_MECA') then
        call jevech('PCACO3D', 'L', jcrf)
        crf = zr(jcrf)
    else
        call jevech('PCACO3D', 'E', jcrf)
    end if
    ivarip = ivarim
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
!
    ndimv = lgpg*npgsn
    b_n = to_blas_int(ndimv)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
!
    call vectan(nb1, nb2, xi, zr(lzr), vecta, &
                vectn, vectpt)
!
    call rccoma(zi(imate), 'ELAS', 1, elasKeyword, valret(1))
    if (elasKeyword .ne. 'ELAS' .and. elasKeyword .ne. 'ELAS_ORTH') then
        call utmess('F', 'PLATE1_12', sk=elasKeyword)
    end if
!
!===============================================================
!     CALCULS DES 2 DDL INTERNES
!
    call trndgl(nb2, vectn, vectpt, zr(ideplm), deplm, &
                rotfcm)
!
    call trndgl(nb2, vectn, vectpt, zr(ideplp), deplp, &
                rotfcp)
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
!
!     CALCUL DE BTDMR, BTDSR : M=MEMBRANE , S=CISAILLEMENT , R=REDUIT
!
            do intsr = 1, npgsr
                call mahsms(0, nb1, xi, ksi3s2, intsr, &
                            zr(lzr), hic, vectn, vectg, vectt, &
                            hsfm, hss)
!
                call hsj1ms(hic, vectg, vectt, hsfm, hss, &
                            hsj1m, hsj1s)
!
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            hic, vectpt, hsj1m, hsj1s, btdm, &
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
                           zr(lzr), hic, vectn, vectg, vectt, &
                           hsf)
!
                call hsj1f(intsn, zr(lzr), hic, vectg, vectt, &
                           hsf, kwgt, hsj1fx, wgt)
!
!     PRODUIT DU POIDS DES PTS DE GAUSS DANS L'EPAISSEUR ET DE WGT
!
                wgt = coef*wgt
!
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), hic, vectpt, hsj1fx, btdf)
!
!     CALCUL DE BTDMN, BTDSN
!     ET
!     FORMATION DE BTILD
!
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
!
!     CALCULS DES COMPOSANTES DE DEFORMATIONS TRIDIMENSIONNELLES :
!     EPSXX, EPSYY, EPSXY, EPSXZ, EPSYZ (CE SONT LES COMPOSANTES TILDE)
                kpgs = kpgs+1
                call epseff('DEFORM', nb1, deplm, btild, x, &
                            epsi, wgt, x)
                eps2d(1) = epsi(1)
                eps2d(2) = epsi(2)
                eps2d(3) = 0.d0
                eps2d(4) = epsi(3)/rac2
!
                call epseff('DEFORM', nb1, deplp, btild, x, &
                            depsi, wgt, x)
                deps2d(1) = depsi(1)
                deps2d(2) = depsi(2)
                deps2d(3) = 0.d0
                deps2d(4) = depsi(3)/rac2
!
                gxz = epsi(4)+depsi(4)
                gyz = epsi(5)+depsi(5)
!
                k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
                k2 = lgpg*(intsn-1)+(npge*(icou-1)+inte-1)*nbvari
                do i = 1, 3
                    sign(i) = zr(icontm-1+k1+i)
                end do
                sign(4) = zr(icontm-1+k1+4)*rac2

                cisail = 0.d0

! ------------- Index of "sub"-point
                ksp = (icou-1)*npge+inte

! ------------- Set main parameters for behaviour (on point)
                call behaviourSetParaPoin(intsn, ksp, BEHinteg)

! ------------- Integrator
                if (elasKeyword .eq. 'ELAS') then
                    sigma = 0.d0
                    call nmcomp(BEHinteg, &
                                fami, intsn, ksp, ndimLdc, typmod, &
                                zi(imate), compor, zr(icarcr), zr(iinstm), zr(iinstp), &
                                4, eps2d, deps2d, 4, sign, &
                                zr(ivarim+k2), option, angmas, &
                                sigma, zr(ivarip+k2), 36, dsidep, cod)
!           COD=1 : ECHEC INTEGRATION LOI DE COMPORTEMENT
!           COD=3 : C_PLAN DEBORST SIGZZ NON NUL
                    if (cod .ne. 0) then
                        if (codret .ne. 1) then
                            codret = cod
                        end if
                        if (cod .eq. 1) goto 999
                    end if
!
                    call rcvalb(fami, intsn, ksp, '+', zi(imate), &
                                ' ', elasKeyword, 0, ' ', [0.d0], &
                                nbv, nomres, valres, valret, 1)
!
                    cisail = valres(1)/(1.d0+valres(2))
!
                else if (elasKeyword .eq. 'ELAS_ORTH') then
                    call moytpg('RIGI', intsn, 3, '+', valpar, &
                                iret)
                    call matrc2(1, 'TEMP    ', [valpar], kappa, matc, &
                                vectt)
                end if
!
!    CALCULS DE LA MATRICE TANGENTE : BOUCLE SUR L'EPAISSEUR
                if (lMatr) then
                    if (elasKeyword .eq. 'ELAS') then
                        dtild(1, 1) = dsidep(1, 1)
                        dtild(1, 2) = dsidep(1, 2)
                        dtild(1, 3) = dsidep(1, 4)/rac2
                        dtild(1, 4) = 0.d0
                        dtild(1, 5) = 0.d0
                        dtild(2, 1) = dsidep(2, 1)
                        dtild(2, 2) = dsidep(2, 2)
                        dtild(2, 3) = dsidep(2, 4)/rac2
                        dtild(2, 4) = 0.d0
                        dtild(2, 5) = 0.d0
                        dtild(3, 1) = dsidep(4, 1)/rac2
                        dtild(3, 2) = dsidep(4, 2)/rac2
                        dtild(3, 3) = dsidep(4, 4)/2.d0
                        dtild(3, 4) = 0.d0
                        dtild(3, 5) = 0.d0
                        dtild(4, 1) = 0.d0
                        dtild(4, 2) = 0.d0
                        dtild(4, 3) = 0.d0
                        dtild(4, 4) = cisail*kappa/2.d0
                        dtild(4, 5) = 0.d0
                        dtild(5, 1) = 0.d0
                        dtild(5, 2) = 0.d0
                        dtild(5, 3) = 0.d0
                        dtild(5, 4) = 0.d0
                        dtild(5, 5) = cisail*kappa/2.d0
                    else if (elasKeyword .eq. 'ELAS_ORTH') then
                        dtild(1, 1) = matc(1, 1)
                        dtild(1, 2) = matc(1, 2)
                        dtild(1, 3) = matc(1, 3)
                        dtild(1, 4) = 0.d0
                        dtild(1, 5) = 0.d0
                        dtild(2, 1) = matc(2, 1)
                        dtild(2, 2) = matc(2, 2)
                        dtild(2, 3) = matc(2, 3)
                        dtild(2, 4) = 0.d0
                        dtild(2, 5) = 0.d0
                        dtild(3, 1) = matc(3, 1)
                        dtild(3, 2) = matc(3, 2)
                        dtild(3, 3) = matc(3, 3)
                        dtild(3, 4) = 0.d0
                        dtild(3, 5) = 0.d0
                        dtild(4, 1) = 0.d0
                        dtild(4, 2) = 0.d0
                        dtild(4, 3) = 0.d0
                        dtild(4, 4) = matc(4, 4)
                        dtild(4, 5) = matc(4, 5)
                        dtild(5, 1) = 0.d0
                        dtild(5, 2) = 0.d0
                        dtild(5, 3) = 0.d0
                        dtild(5, 4) = matc(5, 4)
                        dtild(5, 5) = matc(5, 5)
                    else
                        ASSERT(ASTER_FALSE)
                    end if
!
                    b_n = to_blas_int(25)
                    b_incx = to_blas_int(1)
                    call dscal(b_n, wgt, dtild, b_incx)
!
                    call btkb(5, 42, nddle, dtild, btild, &
                              wmatcb, ktildi)
!
                    do i = 1, nddle
                        do j = 1, nddle
                            ktild(i, j) = ktild(i, j)+ktildi(i, j)
                        end do
                    end do
                end if
!
                if (lSigm) then
                    ASSERT(lVect)
                    if (elasKeyword .eq. 'ELAS') then
                        do i = 1, 3
                            zr(icontp-1+k1+i) = sigma(i)
                        end do
                        zr(icontp-1+k1+4) = sigma(4)/rac2
                        zr(icontp-1+k1+5) = cisail*kappa*gxz/2.d0
                        zr(icontp-1+k1+6) = cisail*kappa*gyz/2.d0
!
!    CALCULS DES EFFORTS INTERIEURS
                        sgmtd(1) = zr(icontp-1+k1+1)
                        sgmtd(2) = zr(icontp-1+k1+2)
                        sgmtd(3) = zr(icontp-1+k1+4)
                        sgmtd(4) = cisail*kappa*gxz/2.d0
                        sgmtd(5) = cisail*kappa*gyz/2.d0
!
                    else if (elasKeyword .eq. 'ELAS_ORTH') then
                        zr(icontp-1+k1+1) = (epsi(1)+depsi(1))*matc(1, 1)+ &
                                            (epsi(2)+depsi(2))*matc(1, 2)+ &
                                            (epsi(3)+depsi(3))*matc(1, 3)
                        zr(icontp-1+k1+2) = (epsi(1)+depsi(1))*matc(2, 1)+ &
                                            (epsi(2)+depsi(2))*matc(2, 2)+ &
                                            (epsi(3)+depsi(3))*matc(2, 3)
                        zr(icontp-1+k1+3) = 0.d0
                        zr(icontp-1+k1+4) = (epsi(1)+depsi(1))*matc(3, 1)+ &
                                            (epsi(2)+depsi(2))*matc(3, 2)+ &
                                            (epsi(3)+depsi(3))*matc(3, 3)
                        zr(icontp-1+k1+5) = matc(4, 4)*gxz+matc(4, 5)*gyz
                        zr(icontp-1+k1+6) = matc(5, 4)*gxz+matc(5, 5)*gyz
!
!    CALCULS DES EFFORTS INTERIEURS
                        sgmtd(1) = zr(icontp-1+k1+1)
                        sgmtd(2) = zr(icontp-1+k1+2)
                        sgmtd(3) = zr(icontp-1+k1+4)
                        sgmtd(4) = dtild(4, 4)*gxz
                        sgmtd(5) = dtild(5, 5)*gyz
                    end if
!
                    call epseff('EFFORI', nb1, x, btild, sgmtd, &
                                x, wgt, effint)
!
                end if
!
            end do
        end do
    end do
!
    if (lMatr) then
!
!     EXPANSION DE LA MATRICE : AJOUTER DE LA ROTATION FICTIVE
!
        nddlet = 6*nb1+3
        call matrkb(nb1, 42, 51, nddlet, ktild, &
                    ctor, rig, crf)
        zr(jcrf) = crf
!
!     AJOUTER DES 3 TRANSLATIONS FICTIVES ASSOCIEES AU NOEUD INTERNE
!     LES 3 TERMES DE RAIDEUR (FICTIVE) ASSOCIEES ONT POUR VALEUR CELLE
!     DES ROTATION FICTIVE
!
    end if
!
    if (lVect) then
!
        call vexpan(nb1, effint, vecl)
!
        do i = 1, 6*nb1
            vecll(i) = vecl(i)
        end do
        vecll(6*nb1+1) = effint(5*nb1+1)
        vecll(6*nb1+2) = effint(5*nb1+2)
        vecll(6*nb1+3) = 0.d0
!
!     CONTRIBUTION DES DDL DE LA ROTATION FICTIVE DANS EFFINT
!
        do i = 1, nb1
            vecll(6*i) = crf*(rotfcm(i)+rotfcp(i))
        end do
        i = nb2
        vecll(6*nb1+3) = crf*(rotfcm(nb2)+rotfcp(nb2))
!     TRANFORMATION DANS REPERE GLOBAL PUIS STOCKAGE
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
        call trnflg(nb2, vecpt, vecll, zr(ivectu))
    end if
!
999 continue
end subroutine
