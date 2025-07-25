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
subroutine enelpg(fami, iadmat, instan, igau, angl_naut, &
                  compor, f, sigma, nbvari, vari, &
                  enelas)
!.......................................................................
    implicit none
! -----------------------------------------------------------------
!
!  OPTION ENEL_ELGA : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  ================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!  EN HPP
!   ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES DE CAUCHY
!            .D         EST LE TENSEUR DE HOOKE
!
!  EN GRANDES DEFORMATIONS SIMO MIEHE POUR ELAS OU VMIS_ISOT
!   ENERLAS = ENERGIE ELASTIQUE SPECIFIQUE
!           = K/2(0.5(J^2-1)-lnJ)+0.5mu(tr(J^(-2/3)be)-3)
!           SI PRESENCE DE THERMIQUE, ON AJOUTE UNE CORRECTION
!           SPECIFIQUE PRESENTEE DANS LA DOC R
!  EN GRANDES DEFORMATIONS GDEF_LOG
!   ENERELAS = SOMME_VOLUME((T_T*(1/D)*T).DV)
!        OU  .T       EST LE TENSEUR DES CONTRAINTES DU FORMALISME
!            .D         EST LE TENSEUR DE HOOKE
! -----------------------------------------------------------------
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
!-----------------------------------------------------------------------
#include "asterfort/d1mamc.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nrsmt1.h"
#include "asterfort/nrsmtb.h"
#include "asterfort/nrsmtt.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/zerop3.h"
#include "blas/dcopy.h"
    integer(kind=8) :: nbsig, nbvari, nsol, i, iadmat, igau, icodre(2), isig, jsig
    real(kind=8) :: c1, c2, deux, vari(*), enelas, trt, un, undemi, zero
    real(kind=8) :: sol(3), sigma(6), jzero, uzero, mzero, instan, epsi(6)
    real(kind=8) :: mjac, ujac, wbe, be(6), e, nu, f(3, 3), angl_naut(3)
    real(kind=8) :: mu, troisk, jac, tau(6), trtau, eqtau, dvtau(6), tlog(6)
    real(kind=8) :: trbe, epsthe, kr(6), pdtsca(6), d1(36), valres(2)
    character(len=4) :: fami
    character(len=16) :: nomres(2), compor(*)
    blas_int :: b_incx, b_incy, b_n
    data kr/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
    data pdtsca/1.d0, 1.d0, 1.d0, 2.d0, 2.d0, 2.d0/
!-----------------------------------------------------------------------
!
    zero = 0.0d0
    deux = 2.0d0
    undemi = 0.5d0
    un = 1.0d0
    nbsig = nbsigm()
    enelas = 0.d0
    nsol = 0
    jzero = zero
    uzero = zero
    ujac = zero
    mzero = zero
    mjac = zero
    wbe = zero
    do i = 1, 3
        sol(i) = zero
    end do
!
!
! --- CAS EN GRANDES DEFORMATIONS SIMO_MIEHE
    if ((compor(3) .eq. 'SIMO_MIEHE') .and. &
        ((compor(1) (1:9) .eq. 'VMIS_ISOT') .or. (compor(1) .eq. 'ELAS'))) then
!
! ---    RECUPERATION DES CARACTERISTIQUES DU MATERIAU :
        nomres(1) = 'E'
        nomres(2) = 'NU'
        call rcvalb(fami, igau, 1, '+', iadmat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 2)
        e = valres(1)
        nu = valres(2)
        mu = e/(2.d0*(1.d0+nu))
        troisk = e/(1.d0-2.d0*nu)
!
        jac = f(1, 1)*(f(2, 2)*f(3, 3)-f(2, 3)*f(3, 2))-f(2, 1)*(f(1, 2)*f(3, 3)-f(1, 3)*f(3, 2)&
              &)+f(3, 1)*(f(1, 2)*f(2, 3)-f(1, 3)*f(2, 2))
!
! ---    CALCUL DE TAU TEL QUE TAU=JAC*SIGMA
!
        tau(5) = 0.d0
        tau(6) = 0.d0
!
        do i = 1, nbsig
            tau(i) = jac*sigma(i)
        end do
!
! ---    CALCUL DE LA TRACE DE TAU- TAU EQUIVALENT ET TAU DEVIATORIQUE
!
        trtau = tau(1)+tau(2)+tau(3)
        eqtau = 0.d0
        do i = 1, 6
            dvtau(i) = tau(i)-kr(i)*trtau/3.d0
            eqtau = eqtau+pdtsca(i)*(dvtau(i)**2.d0)
        end do
        eqtau = sqrt(1.5d0*eqtau)
!
! ---    CALCUL DE LA TRACE DES DEFORMATIONS ELASTIQUES BE
!
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vari(3), b_incx, be, b_incy)
        trbe = be(1)+be(2)+be(3)
        trbe = jac**(-2.d0/3.d0)*(3.d0-2.d0*trbe)
!
!  ---   DEFORMATION THERMIQUE AU POINT D'INTEGRATION COURANT :
!
        call verift(fami, igau, 1, '+', iadmat, &
                    epsth_=epsthe)
!
!
! ---    ATTENTION, EN PRESENCE DE THERMIQUE, CA MET LE BAZARD...
        if (epsthe .ne. 0) then
            call zerop3(-3.d0*epsthe, -1.d0, -3.d0*epsthe, sol, nsol)
            jzero = sol(1)
            call nrsmt1(troisk/3.d0, jzero, uzero)
            call nrsmtt(troisk, jzero, epsthe, mzero)
            call nrsmtt(troisk, jac, epsthe, mjac)
        end if
!
! ---    CALCUL DES TERMES DE L'ENERGIE
        call nrsmt1(troisk/3.d0, jac, ujac)
        call nrsmtb(mu, trbe, wbe)
!
        enelas = ujac+wbe+mjac-uzero-mzero
!
!
! --- CAS EN GRANDES DEFORMATIONS GDEF_LOG
!
    else if ((compor(3) (1:8) .eq. 'GDEF_LOG')) then
!
! ---    RECUPERATION DES CARACTERISTIQUES DU MATERIAU :
        nomres(1) = 'E'
        nomres(2) = 'NU'
        call rcvalb(fami, igau, 1, '+', iadmat, &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 2)
        e = valres(1)
        nu = valres(2)
!
        c1 = (un+nu)/e
        c2 = nu/e
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vari(nbvari-5), b_incx, tlog, b_incy)
!
!        CAS 3D
        if (lteatt('DIM_TOPO_MAILLE', '3')) then
!
            trt = tlog(1)+tlog(2)+tlog(3)
            enelas = undemi*( &
                     tlog(1)*(c1*tlog(1)-c2*trt)+tlog(2)*(c1*tlog(2)-c2*trt)+tlog(3)*(c1*tlog(3)-&
                     &c2*trt)+(tlog(4)*c1*tlog(4)+tlog(5)*c1*tlog(5)+tlog(6)*c1*tlog(6)) &
                     )
!
!
! ---    CAS DES CONTRAINTES PLANES :
        else if (lteatt('C_PLAN', 'OUI')) then
            trt = tlog(1)+tlog(2)
!
            enelas = undemi*( &
                     tlog(1)*(c1*tlog(1)-c2*trt)+tlog(2)*(c1*tlog(2)-c2*trt)+deux*tlog(4)*c1*tlog&
                     &(4) &
                     )
! ---    CAS AXI ET DEFORMATIONS PLANES :
        else
            trt = tlog(1)+tlog(2)+tlog(3)
!
            enelas = undemi*( &
                     tlog(1)*(c1*tlog(1)-c2*trt)+tlog(2)*(c1*tlog(2)-c2*trt)+tlog(3)*(c1*tlog(3)-&
                     &c2*trt)+deux*tlog(4)*c1*tlog(4) &
                     )
        end if
! --- EN HPP SI ON CONSIDERE LE MATERIAU ISOTROPE
! --- E_ELAS = 1/2*SIGMA*1/D*SIGMA :
!
! --- CAS EN GRANDES DEFORMATIONS SIMO_MIEHE
    else if (compor(3) (1:5) .eq. 'PETIT') then
!
!  --    CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE (LE MATERIAU
!  --    POUVANT ETRE ISOTROPE, ISOTROPE-TRANSVERSE OU ORTHOTROPE)
!        ---------------------------------------------------------
        call d1mamc(fami, iadmat, instan, '+', igau, &
                    1, angl_naut, nbsig, d1)
!
!  --    DENSITE D'ENERGIE POTENTIELLE ELASTIQUE AU POINT
!  --    D'INTEGRATION COURANT
!        ---------------------
        call r8inir(6, 0.d0, epsi, 1)
        do isig = 1, nbsig
            do jsig = 1, nbsig
                epsi(isig) = epsi(isig)+d1(nbsig*(isig-1)+jsig)*sigma(jsig)
            end do
            enelas = enelas+undemi*sigma(isig)*epsi(isig)
        end do
!
    else
        call utmess('F', 'COMPOR1_77', sk=compor(3))
    end if
!
end subroutine
