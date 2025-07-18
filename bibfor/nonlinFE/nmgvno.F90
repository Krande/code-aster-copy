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
! aslint: disable=W1504,W0413
!
subroutine nmgvno(fami, ndim, nno1, nno2, npg, &
                  iw, vff1, vff2, idfde1, idfde2, &
                  geom, typmod, option, mat, compor, &
                  lgpg, carcri, instam, instap, ddlm, &
                  ddld, angmas, sigm, vim, sigp, &
                  vip, matr, vect, codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/coefdg.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmepsi.h"
#include "asterfort/nmgvdn.h"
#include "asterfort/nmmabu.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "jeveux.h"
!
    character(len=8) :: typmod(2)
    character(len=*) :: fami
    character(len=16) :: option, compor(COMPOR_SIZE)
    integer(kind=8) :: ndim, nno1, nno2, npg, idfde1, idfde2, iw, mat, lgpg, codret
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg)
    real(kind=8) :: geom(ndim, nno1)
    real(kind=8) :: carcri(CARCRI_SIZE), instam, instap
    real(kind=8) :: ddlm(*), ddld(*), sigm(2*ndim+1, npg), sigp(2*ndim+1, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), matr(*), vect(*)
    real(kind=8) :: angmas(3)
!
! --------------------------------------------------------------------------------------------------
!
!     RAPH_MECA, RIGI_MECA_* ET FULL_MECA_* POUR GRAD_VARI (2D ET 3D)
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DES ELEMENTS
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  IDFDE1  : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE E)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE E)
! IN  IDFDE2  : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE E)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  IW      : POIDS DES POINTS DE GAUSS DE REFERENCE (INDICE)
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODEELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  MAT     : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT PRECEDENT
! IN  INSTAP  : INSTANT DE CALCUL
! IN  TEMPM   : TEMPERATURE AUX NOEUDS A L'INSTANT PRECEDENT
! IN  TEMPP   : TEMPERATURE AUX NOEUDS A L'INSTANT DE CALCUL
! IN  TREF    : TEMPERATURE DE REFERENCE
! IN  DDLM    : DDL A L'INSTANT PRECEDENT
! IN  DDLD    : INCREMENT DES DDL
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  LGPG    : LONGUEUR DU TABLEAU DES VARIABLES INTERNES
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA   ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA   ET FULL_MECA_*)
! OUT MATR    : MATRICE DE RIGIDITE   (RIGI_MECA_* ET FULL_MECA_*)
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR
! MEM DFDI2   :
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: k2(1)
    character(len=16), parameter :: nom(1) = (/'C_GRAD_VARI'/)
    character(len=8), parameter ::  famiNonLocal = "FPG1"
    character(len=8) :: poum
    aster_logical :: grand, axi, lElas, lMatrPred
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: nddl, ndimsi, kpg, cod(27), n, i, m, j, kl, pq, os, osa, kk
    integer(kind=8) :: iu(3*27), ia(8), spt
    real(kind=8) :: rac2, c, val(1)
    real(kind=8) :: deplm(3*27), depld(3*27), dfdi1(27, 3)
    real(kind=8) :: avm, avd, avp, agm(3), agd(3), agp(3), bp
    real(kind=8) :: r, wg, epsm(2*ndim+1), epsd(2*ndim+1), f(3, 3), b(6, 3, 27)
    real(kind=8) :: sigmam(2*ndim+1), sigma(2*ndim+1), dsidep(2*ndim+1, 2*ndim+1), t1, t2
    real(kind=8) :: di, char
    real(kind=8) :: dfdi2(8*3)
    real(kind=8) :: critd(20)
    type(Behaviour_Integ) :: BEHinteg
!
! --------------------------------------------------------------------------------------------------
!

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, carcri, &
                              instam, instap, &
                              fami, mat, &
                              BEHinteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm)
    lElas = option .eq. 'FULL_MECA_ELAS' .or. option .eq. 'RIGI_MECA_ELAS'
    lMatrPred = option .eq. 'RIGI_MECA_TANG'
!
    rac2 = sqrt(2.d0)
    grand = .false.
    axi = typmod(1) .eq. 'AXIS'
    nddl = nno1*ndim+nno2
    ndimsi = 2*ndim

    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(famiNonLocal, kpg, spt, poum, mat, &
                ' ', 'NON_LOCAL', 0, ' ', [0.d0], &
                1, nom, val, k2, 2)
    call coefdg(compor(1), mat, di)
!
    c = val(1)
!
    cod = 0
!
    if (lMatr) then
        call r8inir((nddl*(nddl+1))/2, 0.d0, matr, 1)
    end if
    if (lVect) then
        call r8inir(nddl, 0.d0, vect, 1)
    end if
!
    call nmgvdn(ndim, nno1, nno2, iu, ia)
!
!    EXTRACTION DES DEPLACEMENTS
!
    do n = 1, nno1
        do i = 1, ndim
            deplm(i+(n-1)*ndim) = ddlm(iu(nno1*(i-1)+n))
            if (lMatrPred) then
                depld(i+(n-1)*ndim) = 0.d0
            else
                depld(i+(n-1)*ndim) = ddld(iu(nno1*(i-1)+n))
            end if
        end do
    end do
!
! - CREATION D'UN VECTEUR VALANT 0 POUR ABSENCE DE DEPLACEMENT
!
    do n = 1, nno2
        critd(n) = 0.d0
        do i = 1, ndim
            critd(n) = critd(n)+abs(ddld(iu(nno1*(i-1)+n)))
        end do
    end do
!
! - CALCUL POUR CHAQUE POINT DE GAUSS
!
    do kpg = 1, npg
!
!      CALCUL DES ELEMENTS GEOMETRIQUES DE L'EF POUR A
!
        call dfdmip(ndim, nno2, axi, geom, kpg, &
                    iw, vff2(1, kpg), idfde2, r, wg, &
                    dfdi2)
        avm = 0
        avd = 0
        do n = 1, nno2
            avm = avm+vff2(n, kpg)*ddlm(ia(n))
            avd = avd+vff2(n, kpg)*ddld(ia(n))
            if (lMatrPred) then
                avd = 0.d0
            end if
        end do
        avp = avm+avd
!

!
        do i = 1, ndim
            agm(i) = 0
            agd(i) = 0
            do n = 1, nno2
                agm(i) = agm(i)+dfdi2(nno2*(i-1)+n)*ddlm(ia(n))
                agd(i) = agd(i)+dfdi2(nno2*(i-1)+n)*ddld(ia(n))
                if (lMatrPred) then
                    agd(i) = 0.d0
                end if
            end do
            agp(i) = agm(i)+agd(i)
        end do
!
!      CALCUL DES ELEMENTS GEOMETRIQUES DE L'EF POUR U
!
        call dfdmip(ndim, nno1, axi, geom, kpg, &
                    iw, vff1(1, kpg), idfde1, r, wg, &
                    dfdi1)

!
        call nmepsi(ndim, nno1, axi, grand, vff1(1, kpg), &
                    r, dfdi1, deplm, f, epsm(1:ndimsi))
        call nmepsi(ndim, nno1, axi, grand, vff1(1, kpg), &
                    r, dfdi1, depld, f, epsd(1:ndimsi))
        call nmmabu(ndim, nno1, .false._1, grand, dfdi1, &
                    b)
        if (axi) then
            do n = 1, nno1
                b(3, 1, n) = vff1(n, kpg)/r
            end do
        end if

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHinteg)

! ----- Integrator
        do kl = 1, 3
            sigmam(kl) = sigm(kl, kpg)
        end do
        do kl = 4, ndimsi
            sigmam(kl) = sigm(kl, kpg)*rac2
        end do
        sigmam(ndimsi+1) = sigm(ndimsi+1, kpg)
        epsm(2*ndim+1) = avm
        epsd(2*ndim+1) = avd

        call nmcomp(BEHinteg, &
                    fami, kpg, ksp, ndim, typmod, &
                    mat, compor, carcri, instam, instap, &
                    ndimsi+1, epsm, epsd, ndimsi+1, sigmam, &
                    vim(1, kpg), option, angmas, &
                    sigma, vip(1, kpg), (ndimsi+1)*(ndimsi+1), dsidep, cod(kpg))
!
        if (cod(kpg) .eq. 1) goto 900
!
!      FORCE INTERIEURE ET CONTRAINTES DE CAUCHY
!
        if (lSigm) then
!
!        CONTRAINTES
!
            do kl = 1, 3
                sigp(kl, kpg) = sigma(kl)
            end do
            do kl = 4, ndimsi
                sigp(kl, kpg) = sigma(kl)/rac2
            end do
            sigp(ndimsi+1, kpg) = sigma(ndimsi+1)
        end if

        if (lVect) then
            bp = sigp(ndimsi+1, kpg)
!
!        VECTEUR FINT:U
!
            do n = 1, nno1
                do i = 1, ndim
                    kk = iu(nno1*(i-1)+n)
                    t1 = 0
                    do kl = 1, ndimsi
                        t1 = t1+sigma(kl)*b(kl, i, n)
                    end do
                    vect(kk) = vect(kk)+wg*t1
                end do
            end do
!
!        VECTEUR FINT:A
!
            do n = 1, nno2
                t1 = vff2(n, kpg)*bp
                t2 = 0
                do i = 1, ndim
                    t2 = t2+c*dfdi2(nno2*(i-1)+n)*agp(i)
                end do
                kk = ia(n)
                vect(kk) = vect(kk)+wg*(t2+t1)
            end do
        end if
!
!   CALCUL DE LA MATRICE DE RIGIDITE
!   STOCKAGE TRIANGLE INFERIEUR LIGNE DE DFI/DUJ
!
        if (lMatr) then
!
!        MATRICE K:U(I,N),U(J,M)
!
            do n = 1, nno1
                do i = 1, ndim
                    os = ((iu(nno1*(i-1)+n)-1)*iu(nno1*(i-1)+n))/2
                    do m = 1, nno1
                        do j = 1, ndim
                            if (iu(nno1*(j-1)+m) .gt. iu(nno1*(i-1)+n)) goto 821
                            kk = os+iu(nno1*(j-1)+m)
                            t1 = 0
                            do kl = 1, ndimsi
                                do pq = 1, ndimsi
                                    t1 = t1+dsidep(kl, pq)*b(pq, j, m)*b(kl, i, n)
                                end do
                            end do
                            matr(kk) = matr(kk)+wg*t1
                        end do
                    end do
821                 continue
                end do
            end do
!
!        MATRICES K:A(N),A(M) SI ENDO NON-NUL
!
            do n = 1, nno2
                osa = ((ia(n)-1)*ia(n))/2
                do m = 1, nno2
                    t1 = vff2(n, kpg)*vff2(m, kpg)*dsidep(ndimsi+1, ndimsi+1)
                    t2 = 0
                    do i = 1, ndim
                        t2 = t2+dfdi2(nno2*(i-1)+n)*dfdi2(nno2*(i-1)+m)
                    end do
                    t2 = c*t2
                    if (ia(m) .le. ia(n)) then
                        kk = osa+ia(m)
                        matr(kk) = matr(kk)+wg*(t2+t1)
                    end if
                end do
            end do
        end if
!
        if (lMatr .and. .not. lElas) then
!
!        MATRICES K:A(N),U(J,M)
!
            do n = 1, nno2
                do m = 1, nno1
                    do j = 1, ndim
                        t1 = 0
                        do kl = 1, ndimsi
                            t1 = t1+dsidep(kl, ndimsi+1)*b(kl, j, m)
                        end do
                        t1 = vff2(n, kpg)*t1
                        if (ia(n) .ge. iu(nno1*(j-1)+m)) then
                            kk = ((ia(n)-1)*ia(n))/2+iu(nno1*(j-1)+m)
                        else
                            kk = ((iu(nno1*(j-1)+m)-1)*iu(nno1*(j-1)+m))/2+ia(n)
                        end if
                        matr(kk) = matr(kk)+wg*t1
                    end do
                end do
            end do
!
        end if
!
        if (lMatr) then
!
            do n = 1, nno2
                osa = ((ia(n)-1)*ia(n))/2
                do m = 1, nno2
                    if (ia(m) .le. ia(n)) then
                        kk = osa+ia(m)
!
                        char = ddld(ia(m))
!
                        if (char .eq. 0.d0 .and. critd(m) .ne. 0.d0) then
                            if (ia(m) .eq. ia(n)) then
                                matr(kk) = di
                            else
                                matr(kk) = 0.d0
                            end if
                        end if
!
                        char = ddld(ia(n))
!
                        if (char .eq. 0.d0 .and. critd(n) .ne. 0.d0) then
                            if (ia(m) .eq. ia(n)) then
                                matr(kk) = di
                            else
                                matr(kk) = 0.d0
                            end if
                        end if
!
                        if (ia(m) .eq. ia(n)) then
                            if (critd(n) .eq. 0.d0 .and. avp .eq. 0.d0) then
                                matr(kk) = di
                            end if
                        end if
                    end if
                end do
            end do
!
        end if
!
        if (lMatr .and. .not. lElas) then
!
            do n = 1, nno2
!
                char = ddld(ia(n))
!
                if (char .eq. 0.d0 .and. critd(n) .ne. 0.d0) then
                    do m = 1, nno1
                        do j = 1, ndim
                            if (ia(n) .ge. iu(nno1*(j-1)+m)) then
                                kk = ((ia(n)-1)*ia(n))/2+iu(nno1*(j-1) &
                                                            +m)
                            else
                                kk = ((iu(nno1*(j-1)+m)-1)*iu(nno1*(j-1) &
                                                              +m))/2+ia(n)
                            end if
                            matr(kk) = 0.d0
                        end do
                    end do
                end if
            end do
        end if
    end do
!
900 continue
!
    call codere(cod, npg, codret)
!
end subroutine
