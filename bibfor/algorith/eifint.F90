! --------------------------------------------------------------------
! Copyright (C) 1991 - 2024 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1504,W1306
!
subroutine eifint(ndim, axi, nno1, nno2, npg, &
                  wref, vff1, vff2, dffr2, geom, &
                  ang, typmod, option, mat, compor, &
                  lgpg, carcri, instam, instap, ddlm, &
                  ddld, iu, im, vim, sigp, &
                  vip, matr, vect, &
                  lMatr, lVect, lSigm, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/codere.h"
#include "asterfort/eicine.h"
#include "asterfort/nmcomp.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=8) :: typmod(2)
    character(len=16) :: option, compor(COMPOR_SIZE)
    aster_logical :: axi
    integer :: ndim, nno1, nno2, npg, mat, lgpg, iu(3, 18), im(3, 9)
    integer :: codret
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno2)
    real(kind=8) :: wref(npg)
    real(kind=8) :: carcri(CARCRI_SIZE), instam, instap
    real(kind=8) :: ddlm(2*nno1*ndim+nno2*ndim), ddld(2*nno1*ndim+nno2*ndim)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), vect(2*nno1*ndim+nno2*ndim)
    real(kind=8) :: dffr2(ndim-1, nno2, npg), ang(*), sigp(2*ndim, npg), matr(*)
    aster_logical, intent(in) :: lMatr, lVect, lSigm
!
! --------------------------------------------------------------------------------------------------
!
!   RAPH_MECA, RIGI_MECA_* ET FULL_MECA_* POUR L'ELEMENT D'INTERFACE
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  AXI     : .TRUE. SI AXISYMETRIE
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  IDFF1   : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE X)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE X)
! IN  DFFR2   : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE L)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  WREF    : POIDS DES POINTS DE GAUSS DE REFERENCE
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  MAT     : MATERIAU CODE
! IN  COMPOR  : COMPORTEMENT
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT PRECEDENT
! IN  INSTAP  : INSTANT DE CALCUL
! IN  DDLM    : DDL A L'INSTANT PRECEDENT
! IN  DDLD    : INCREMENT DES DDL
! IN  IU      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! IN  IM      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  LGPG    : LONGUEUR DU TABLEAU DES VARIABLES INTERNES
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA   ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA   ET FULL_MECA_*)
! OUT MATR    : MATRICE DE RIGIDITE   (RIGI_MECA_* ET FULL_MECA_*)
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    type(Behaviour_Integ) :: BEHinteg
    character(len=4), parameter :: fami = 'RIGI'
    integer, parameter :: ksp = 1
    integer :: nddl, g, cod(npg), n, i, m, j, k, l, os, kk
    real(kind=8) :: angmas(3), r, wg, b(3, 3, 18), t1
    real(kind=8) :: mu_m(ndim), su_m(ndim), mu_d(ndim), su_d(ndim), mu_p(ndim), su_p(ndim)
    real(kind=8) :: eps_m(2*ndim), eps_d(2*ndim), sig_m(2*ndim), de(2*ndim), ddedt(2*ndim, 2*ndim)
    real(kind=8) :: cl(ndim)
!
! TODO: passer de(3) et ddedt(3,3) avec impact sur les lois de comportement
! cf. exemple lc7077
! --------------------------------------------------------------------------------------------------
!
    nddl = nno1*2*ndim+nno2*ndim
    cod = 0
    codret = 0
    sig_m = 0
    angmas = 0

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, carcri, &
                              instam, instap, &
                              fami, mat, &
                              BEHinteg)

    if (lMatr) matr(1:(nddl*(nddl+1))/2) = 0.d0
    if (lVect) vect(1:nddl) = 0.d0

!
! - Loop on Gauss points
!
    do g = 1, npg
!
!      CALCUL DES ELEMENTS GEOMETRIQUES DE L'EF
!
        call eicine(ndim, axi, nno1, nno2, vff1(1, g), &
                    vff2(1, g), wref(g), dffr2(1, 1, g), geom, ang, &
                    wg, b)

        su_m = 0
        su_d = 0
        do i = 1, ndim
            do j = 1, ndim
                do n = 1, 2*nno1
                    su_m(i) = su_m(i)+b(i, j, n)*ddlm(iu(j, n))
                    su_d(i) = su_d(i)+b(i, j, n)*ddld(iu(j, n))
                end do
            end do
        end do
        su_p = su_m+su_d
!
        mu_m = 0
        mu_d = 0
        do i = 1, ndim
            do n = 1, nno2
                mu_m(i) = mu_m(i)+vff2(n, g)*ddlm(im(i, n))
                mu_d(i) = mu_d(i)+vff2(n, g)*ddld(im(i, n))
            end do
        end do
        mu_p = mu_m+mu_d

! ----- Integration of behaviour
!      conventions :
!       1. mu est range dans eps(1:ndim)
!       2. su est range dans eps(ndim+1:2*ndim)
!       3. delta est renvoye dans sig(1:ndim)           : de
!       4. d(delta)/dt est renvoye dans dsidep(1:3,1:3) : ddedt
!       5. r (penalisation) est renvoye dans tampon(1)  : r

        eps_m(1:ndim) = mu_m
        eps_m(ndim+1:2*ndim) = su_m
        eps_d(1:ndim) = mu_d
        eps_d(ndim+1:2*ndim) = su_d

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(g, ksp, BEHinteg)

! ----- Integrator
        call nmcomp(BEHinteg, &
                    fami, g, ksp, ndim, typmod, &
                    mat, compor, carcri, instam, instap, &
                    2*ndim, eps_m, eps_d, 2*ndim, sig_m, &
                    vim(1, g), option, angmas, &
                    de, vip(1, g), 2*ndim*2*ndim, ddedt, cod(g))
        if (cod(g) .eq. 1) goto 999

        r = BEHinteg%behavESVA%behavESVAOther%r
        cl = su_p(1:ndim)-de(1:ndim)

! ----- Stresses
        if (lSigm) then
            sigp(1:ndim, g) = mu_p+r*cl
            sigp(ndim+1:2*ndim, g) = cl
        end if

! ----- Internal forces
        if (lVect) then
! --------- Vector (DOF: U)
            do n = 1, 2*nno1
                do i = 1, ndim
                    kk = iu(i, n)
                    t1 = dot_product(b(1:ndim, i, n), mu_p+r*cl)
                    vect(kk) = vect(kk)+wg*t1
                end do
            end do

! --------- Vector (DOF: MU)
            do n = 1, nno2
                do i = 1, ndim
                    kk = im(i, n)
                    t1 = vff2(n, g)*cl(i)
                    vect(kk) = vect(kk)+wg*t1
                end do
            end do
        end if

! ----- Matrix
        if (lMatr) then
! --------- Matrix [U(I,N),U(J,M)]
            do n = 1, 2*nno1
                do i = 1, ndim
                    os = ((iu(i, n)-1)*iu(i, n))/2
                    do m = 1, 2*nno1
                        do j = 1, ndim
                            if (iu(j, m) .gt. iu(i, n)) cycle
                            kk = os+iu(j, m)
                            t1 = 0
                            do k = 1, ndim
                                t1 = t1+b(k, i, n)*b(k, j, m)
                                do l = 1, ndim
                                    t1 = t1-b(k, i, n)*r*ddedt(k, l)*b(l, j, m)
                                end do
                            end do
                            t1 = t1*r
                            matr(kk) = matr(kk)+wg*t1
                        end do
                    end do
                end do
            end do

! --------- Matrix [MU(I,N),U(J,M)]
            do n = 1, nno2
                do i = 1, ndim
                    do m = 1, 2*nno1
                        do j = 1, ndim
                            if (im(i, n) .ge. iu(j, m)) then
                                kk = ((im(i, n)-1)*im(i, n))/2+iu(j, m)
                            else
                                kk = ((iu(j, m)-1)*iu(j, m))/2+im(i, n)
                            end if
                            t1 = vff2(n, g)*b(i, j, m)
                            do l = 1, ndim
                                t1 = t1-vff2(n, g)*r*ddedt(i, l)*b(l, j, m)
                            end do
                            matr(kk) = matr(kk)+wg*t1
                        end do
                    end do
                end do
            end do

! --------- Matrix [MU(I,N),MU(J,M)]
            do n = 1, nno2
                do i = 1, ndim
                    os = ((im(i, n)-1)*im(i, n))/2
                    do m = 1, nno2
                        do j = 1, ndim
                            if (im(j, m) .gt. im(i, n)) cycle
                            kk = os+im(j, m)
                            t1 = -vff2(n, g)*ddedt(i, j)*vff2(m, g)
                            matr(kk) = matr(kk)+wg*t1
                        end do
                    end do
                end do
            end do
        end if
    end do
!
999 continue
!
! - Save return code
!
    if (lSigm) then
        call codere(cod, npg, codret)
    end if
!
end subroutine
