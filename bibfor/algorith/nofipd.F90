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
! aslint: disable=W1306,W1504
!
subroutine nofipd(ndim, nnod, nnop, nnog, npg, &
                  iw, vffd, vffp, vffg, idffd, &
                  vu, vp, vpi, geomi, typmod, &
                  option, nomte, mate, compor, lgpg, &
                  carcri, instm, instp, ddlm, ddld, &
                  angmas, sigm, vim, sigp, vip, &
                  vect, matr, codret, lSigm, lVect, &
                  lMatr)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codere.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmepsi.h"
#include "asterfort/tanbul.h"
#include "asterfort/uthk.h"
#include "blas/ddot.h"
#include "asterfort/Behaviour_type.h"
!
    integer(kind=8) :: ndim, nnod, nnop, nnog, npg, iw, idffd, lgpg
    integer(kind=8) :: mate
    integer(kind=8) :: vu(3, 27), vp(27), vpi(3, 27)
    integer(kind=8) :: codret
    real(kind=8) :: vffd(nnod, npg), vffp(nnop, npg), vffg(nnog, npg)
    real(kind=8) :: instm, instp
    real(kind=8) :: geomi(ndim, nnod), ddlm(*), ddld(*), angmas(*)
    real(kind=8) :: sigm(2*ndim+1, npg), sigp(2*ndim+1, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    real(kind=8) :: vect(*), matr(*)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=8), intent(in)  :: typmod(2)
    character(len=16), intent(in)  :: compor(COMPOR_SIZE), option
    character(len=16) :: nomte
    aster_logical, intent(in) :: lSigm, lVect, lMatr
!
! --------------------------------------------------------------------------------------------------
!
!          CALCUL DES FORCES INTERNES POUR LES ELEMENTS
!          INCOMPRESSIBLES POUR LES PETITES DEFORMATIONS
!          3D/D_PLAN/AXIS
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AUX DEPLACEMENTS
! IN  NNO2    : NOMBRE DE NOEUDS DE L'ELEMENT LIES AU GONFLEMENT
! IN  NNO3    : NOMBRE DE NOEUDS DE L'ELEMENT LIES A LA PRESSION
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS
! IN  VFF1    : VALEUR  DES FONCTIONS DE FORME LIES AUX DEPLACEMENTS
! IN  VFF2    : VALEUR  DES FONCTIONS DE FORME LIES AU GONFLEMENT
! IN  VFF3    : VALEUR  DES FONCTIONS DE FORME LIES A LA PRESSION
! IN  IDFF1   : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  VU      : TABLEAU DES INDICES DES DDL DE DEPLACEMENTS
! IN  VP      : TABLEAU DES INDICES DES DDL DE PRESSION
! IN  VIP     : TABLEAU DES INDICES DES DDL DE GRADIENT DE PRESSION PROJETE
! IN  GEOMI   : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  MATE    : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTM   : INSTANT PRECEDENT
! IN  INSTP   : INSTANT DE CALCUL
! IN  DDLM    : DEGRES DE LIBERTE A L'INSTANT PRECEDENT
! IN  DDLD    : INCREMENT DES DEGRES DE LIBERTE
! IN  ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT VECT    : FORCES INTERNES
! OUT MATR    : MATRICE DE RIGIDITE (RIGI_MECA_TANG ET FULL_MECA)
! OUT CODRET  : CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    character(len=4), parameter :: fami = "RIGI"
    aster_logical, parameter :: mini = ASTER_FALSE, grand = ASTER_FALSE
    aster_logical :: axi
    integer(kind=8) :: kpg, nddl
    integer(kind=8) :: ia, na, ra, sa, ib, nb, rb, sb, ja, jb
    integer(kind=8) :: os, kk
    integer(kind=8) :: vuiana, vpiana, vpsa
    integer(kind=8) :: cod(npg)
    character(len=16) :: rela_comp
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: deplm(3*27), depld(3*27)
    real(kind=8) :: r, w, dff1(nnod, ndim)
    real(kind=8) :: presm(27), presd(27)
    real(kind=8) :: gpresm(3*27), gpresd(3*27)
    real(kind=8) :: pm, pd, gpm(ndim), gpd(ndim), pim(ndim), pid(ndim)
    real(kind=8) :: fm(3, 3), epsm(6), deps(6)
    real(kind=8) :: sigma(6), sigmPrep(6), sigtr
    real(kind=8) :: dsidep(6, 6)
    real(kind=8) :: def(2*ndim, nnod, ndim), deftr(nnod, ndim), ddivu, divum
    real(kind=8) :: ddev(6, 6), devd(6, 6), dddev(6, 6)
    real(kind=8) :: t1, t2
    real(kind=8) :: alpha, trepst
    real(kind=8) :: dsbdep(2*ndim, 2*ndim)
    real(kind=8) :: stab, hk
    type(Behaviour_Integ) :: BEHinteg
    real(kind=8), parameter :: idev(6, 6) = reshape((/2.d0, -1.d0, -1.d0, 0.d0, 0.d0, 0.d0, &
                                                      -1.d0, 2.d0, -1.d0, 0.d0, 0.d0, 0.d0, &
                                                      -1.d0, -1.d0, 2.d0, 0.d0, 0.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 3.d0, 0.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 3.d0, 0.d0, &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 3.d0/), &
                                                    (/6, 6/))
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    cod = 0
    nddl = nnod*ndim+nnop+nnog*ndim
    dsidep = 0.d0
    codret = 0
    if (lVect) then
        vect(1:nddl) = 0.d0
    end if
    if (lMatr) then
        matr(1:nddl*(nddl+1)/2) = 0.d0
    end if

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, carcri, &
                              instm, instp, &
                              fami, mate, &
                              BEHinteg)

! - Compute stabilization
    call uthk(nomte, geomi, hk, ndim, 1)
    stab = 1.d-4*hk*hk
!
! - Extract for fields
!
    do na = 1, nnod
        do ia = 1, ndim
            deplm(ia+ndim*(na-1)) = ddlm(vu(ia, na))
            depld(ia+ndim*(na-1)) = ddld(vu(ia, na))
        end do
    end do
    do sa = 1, nnop
        presm(sa) = ddlm(vp(sa))
        presd(sa) = ddld(vp(sa))
    end do
    do ra = 1, nnog
        do ia = 1, ndim
            gpresm(ia+ndim*(ra-1)) = ddlm(vpi(ia, ra))
            gpresd(ia+ndim*(ra-1)) = ddld(vpi(ia, ra))
        end do
    end do
!
! - Properties of behaviour
!
    rela_comp = compor(RELA_NAME)

! - Loop on Gauss points
    do kpg = 1, npg
        epsm = 0.d0
        deps = 0.d0
! ----- Kinematic - Previous strains
        call dfdmip(ndim, nnod, axi, geomi, kpg, &
                    iw, vffd(1, kpg), idffd, r, w, &
                    dff1)
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dff1, deplm, fm, epsm)
        divum = epsm(1)+epsm(2)+epsm(3)
! ----- Kinematic - Increment of strains
        call nmepsi(ndim, nnod, axi, grand, vffd(1, kpg), &
                    r, dff1, depld, fm, deps)
        ddivu = deps(1)+deps(2)+deps(3)
! ----- Pressure and "gonflement"
        b_n = to_blas_int(nnop)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pm = ddot(b_n, vffp(1, kpg), b_incx, presm, b_incy)
        b_n = to_blas_int(nnop)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        pd = ddot(b_n, vffp(1, kpg), b_incx, presd, b_incy)
        do ia = 1, ndim
            b_n = to_blas_int(nnog)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(ndim)
            pim(ia) = ddot(b_n, vffg(1, kpg), b_incx, gpresm(ia), b_incy)
            b_n = to_blas_int(nnog)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(ndim)
            pid(ia) = ddot(b_n, vffg(1, kpg), b_incx, gpresd(ia), b_incy)
            b_n = to_blas_int(nnop)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            gpm(ia) = ddot(b_n, dff1(1, ia), b_incx, presm, b_incy)
            b_n = to_blas_int(nnop)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            gpd(ia) = ddot(b_n, dff1(1, ia), b_incx, presd, b_incy)
        end do
! ----- Kinematic - Product [F].[B]
        if (ndim .eq. 2) then
            do na = 1, nnod
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = 0.d0
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                end do
            end do
            if (axi) then
                do na = 1, nnod
                    def(3, na, 1) = fm(3, 3)*vffd(na, kpg)/r
                end do
            end if
        else if (ndim .eq. 3) then
            do na = 1, nnod
                do ia = 1, ndim
                    def(1, na, ia) = fm(ia, 1)*dff1(na, 1)
                    def(2, na, ia) = fm(ia, 2)*dff1(na, 2)
                    def(3, na, ia) = fm(ia, 3)*dff1(na, 3)
                    def(4, na, ia) = (fm(ia, 1)*dff1(na, 2)+fm(ia, 2)*dff1(na, 1))/rac2
                    def(5, na, ia) = (fm(ia, 1)*dff1(na, 3)+fm(ia, 3)*dff1(na, 1))/rac2
                    def(6, na, ia) = (fm(ia, 2)*dff1(na, 3)+fm(ia, 3)*dff1(na, 2))/rac2
                end do
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
! ----- CALCUL DE TRACE(B)
        do na = 1, nnod
            do ia = 1, ndim
                deftr(na, ia) = def(1, na, ia)+def(2, na, ia)+def(3, na, ia)
            end do
        end do
! ----- Prepare stresses
        do ia = 1, 3
            sigmPrep(ia) = sigm(ia, kpg)+sigm(2*ndim+1, kpg)
        end do
        do ia = 4, 2*ndim
            sigmPrep(ia) = sigm(ia, kpg)*rac2
        end do

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHinteg)

! ----- Integrator
        sigma = 0.d0
        call nmcomp(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    mate, compor, carcri, instm, instp, &
                    6, epsm, deps, 6, sigmPrep, &
                    vim(1, kpg), option, angmas, &
                    sigma, vip(1, kpg), 36, dsidep, cod(kpg))
        if (cod(kpg) .eq. 1) then
            goto 999
        end if
! ----- Compute "bubble" matrix
        call tanbul(option, ndim, kpg, mate, rela_comp, &
                    lVect, mini, alpha, dsbdep, trepst)
! ----- Internal forces
        if (lVect) then
            sigtr = sigma(1)+sigma(2)+sigma(3)
            do ia = 1, 3
                sigma(ia) = sigma(ia)-sigtr/3.d0+(pm+pd)
            end do
            do na = 1, nnod
                do ia = 1, ndim
                    kk = vu(ia, na)
                    b_n = to_blas_int(2*ndim)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    t1 = ddot(b_n, sigma, b_incx, def(1, na, ia), b_incy)
                    vect(kk) = vect(kk)+w*t1
                end do
            end do
            t2 = (divum+ddivu-(pm+pd)*alpha-trepst)
            do sa = 1, nnop
                kk = vp(sa)
                t1 = 0.d0
                do ia = 1, ndim
                    t1 = t1+dff1(sa, ia)*(gpm(ia)+gpd(ia)-pim(ia)-pid(ia))
                end do
                t1 = vffp(sa, kpg)*t2-stab*t1
                vect(kk) = vect(kk)+w*t1
            end do
            do ra = 1, nnog
                do ia = 1, ndim
                    kk = vpi(ia, ra)
                    t1 = stab*vffg(ra, kpg)*(gpm(ia)+gpd(ia)-pim(ia)-pid(ia))
                    vect(kk) = vect(kk)+w*t1
                end do
            end do
        end if
! ----- Cauchy stresses
        if (lSigm) then
            do ia = 1, 3
                sigp(ia, kpg) = sigma(ia)
            end do
            do ia = 4, 2*ndim
                sigp(ia, kpg) = sigma(ia)/rac2
            end do
            sigp(2*ndim+1, kpg) = sigtr/3.d0-(pm+pd)
        end if
! ----- Rigidity matrix
        if (lMatr) then
            devd = matmul(idev/3.d0, dsidep)
            ddev = matmul(dsidep, idev/3.d0)
            dddev = matmul(devd, idev/3.d0)
! - TERME K:UX
            do na = 1, nnod
                do ia = 1, ndim
                    vuiana = vu(ia, na)
                    os = (vuiana-1)*vuiana/2
! - TERME K:UU      KUU(NDIM,NNO1,NDIM,NNO1)
                    do nb = 1, nnod
                        do ib = 1, ndim
                            if (vu(ib, nb) .le. vuiana) then
                                kk = os+vu(ib, nb)
                                t1 = 0.d0
                                do ja = 1, 2*ndim
                                    do jb = 1, 2*ndim
                                        t1 = t1+def(ja, na, ia)*dddev(ja, jb)*def(jb, nb, ib)
                                    end do
                                end do
                                matr(kk) = matr(kk)+w*t1
                            end if
                        end do
                    end do
! - TERME K:UP      KUP(NDIM,NNO1,NNO2)
                    do sb = 1, nnop
                        if (vp(sb) .lt. vuiana) then
                            kk = os+vp(sb)
                            t1 = deftr(na, ia)*vffp(sb, kpg)
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
! - TERME K:UPI = 0.D0     KUP(NDIM,NNO1,NDIM,NNO3)
                end do
            end do
! - TERME K:PX
            do sa = 1, nnop
                vpsa = vp(sa)
                os = (vpsa-1)*vpsa/2
! - TERME K:PU      KPU(NDIM,NNO2,NNO1)
                do nb = 1, nnod
                    do ib = 1, ndim
                        if (vu(ib, nb) .lt. vpsa) then
                            kk = os+vu(ib, nb)
                            t1 = vffp(sa, kpg)*deftr(nb, ib)
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
                end do
! - TERME K:PP      KPP(NNO2,NNO2)
                do sb = 1, nnop
                    if (vp(sb) .le. vpsa) then
                        kk = os+vp(sb)
                        t1 = -vffp(sa, kpg)*vffp(sb, kpg)*alpha
                        t2 = 0.d0
! - PRODUIT SCALAIRE DES GRAD DE FONCTIONS DE FORME
                        do ia = 1, ndim
                            t2 = t2-dff1(sa, ia)*dff1(sb, ia)
                        end do
                        matr(kk) = matr(kk)+w*(t2*stab+t1)
                    end if
                end do
! - TERME K:PPI     KPPI(NNO2,NDIM,NNO3)
                do rb = 1, nnog
                    do ib = 1, ndim
                        if (vpi(ib, rb) .lt. vpsa) then
                            kk = os+vpi(ib, rb)
                            t1 = vffg(rb, kpg)*deftr(sa, ib)*stab
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
                end do
            end do
! - TERME K:PIX
            do ra = 1, nnog
                do ia = 1, ndim
                    vpiana = vpi(ia, ra)
                    os = (vpiana-1)*vpiana/2
!
! - TERME K:PIU = 0.D0     KUU(NDIM,NNO3,NDIM,NNO1)
!
! - TERME K:PIP     KPIP(NDIM,NNO3,NNO2)
                    do sb = 1, nnop
                        if (vp(sb) .lt. vpiana) then
                            kk = os+vp(sb)
                            t1 = vffg(ra, kpg)*deftr(sb, ia)*stab
                            matr(kk) = matr(kk)+w*t1
                        end if
                    end do
!
! - TERME K:PIPI    KPIPI(NDIM,NNO3,NDIM,NNO3)
! - REMARQUE : MATRICE NON NULLE QUE SI I=J DONC K:PIPI(I,N,J,M)=0 SINON
                    do rb = 1, nnog
                        do ib = 1, ndim
                            if (vpi(ib, rb) .le. vpiana .and. (ia .eq. ib)) then
                                kk = os+vpi(ib, rb)
                                t1 = -vffg(ra, kpg)*vffg(rb, kpg)*stab
                                matr(kk) = matr(kk)+w*t1
                            end if
                        end do
                    end do
                end do
            end do
        end if
    end do
!
999 continue

! - Return code summary
    call codere(cod, npg, codret)
!
end subroutine
