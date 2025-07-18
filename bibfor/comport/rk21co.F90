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

subroutine rk21co(fami, kpg, ksp, rela_comp, mod, &
                  imat, matcst, nbcomm, cpmono, nfs, &
                  nsg, toutms, nvi, nmat, y, &
                  kp, ee, a, h, pgl, &
                  nbphas, cothe, coeff, dcothe, dcoeff, &
                  coel, x, pas, neps, epsd, &
                  detot, nhsr, numhsr, hsr, itmax, &
                  toler, iret)
! aslint: disable=W1306,W1504
    implicit none
!     INTEGRATION DE LOIS DE COMPORTEMENT ELASTO-VISCOPLASTIQUE
!     PAR UNE METHODE DE RUNGE KUTTA D'ORDRE 2 A 2 PAS EMBOITES
!
!     CALCUL DE LA SOLUTION A L ORDRE 1 ET A L ORDRE 2
!     ----------------------------------------------------------------
!     IN  FAMI    :  FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!         KPG,KSP :  NUMERO DU (SOUS)POINT DE GAUSS
!         RELA_COMP    :  NOM DU MODELE DE COMPORTEMENT
!         MOD     :  TYPE DE MODELISATION
!         IMAT    :  CODE DU MATERIAU CODE
!         MATCST  : 'OUI'  'NAP'  'NON'
!         NVI     :  NOMBRE DE VARIABLES INTERNES
!         NMAT    :  NOMBRE DE PARAMETRES MATERIAU INELASTIQUE
!     VAR Y       :  VARIABLES INTERNES
!         KP      :  INDICE POUR L'INTEGRATION
!                  KP=1 AUGMENTATION DU PAS DE TEMPS
!                  KP=0 DIMINUTION DU PAS DE TEMPS
!     VAR EE      :  ERREUR=DIFF ENTRE Y2 ET Y1
!     VAR A       :  F(Y,T)
!         H       :  PAS DE TEMPS TESTE
!         COTHE   :  COEFFICIENTS MATERIAU ELAS A T
!         COEFF   :  COEFFICIENTS MATERIAU INELAS A T
!         DCOTHE  :  COEFFICIENTS MATERIAU ELAS A T+DT
!         DCOEFF  :  COEFFICIENTS MATERIAU INELAS A T+DT
!     VAR X       :  INTERVALE DE TEMPS ADAPTATIF
!         PAS     :  PAS DE TEMPS
!         EPSD    :  DEFORMATION TOTALE A T
!         DETOT   :  INCREMENT DE DEFORMATION TOTALE
!     ----------------------------------------------------------------
!
#include "asterfort/rdif01.h"
    integer(kind=8) :: kpg, ksp, nmat, imat, nbcomm(nmat, 3), kp, nvi, i, nfs, nsg
    integer(kind=8) :: nbphas, itmax, iret, nhsr, numhsr(*), neps
    character(len=16) :: rela_comp
    character(len=24) :: cpmono(5*nmat+1)
    character(len=8) :: mod
    character(len=*) :: fami
    character(len=3) :: matcst
    real(kind=8) :: pgl(3, 3), coel(nmat), cothe(nmat), dcothe(nmat)
    real(kind=8) :: x, pas, h, hs2, epsd(6), detot(6), y(nvi)
    real(kind=8) :: f(nvi), hsr(nsg, nsg, nhsr), toler
    real(kind=8) :: coeff(nmat), dcoeff(nmat), ee(nvi), a(nvi)
!      POUR GAGNER EN TEMPS CPU
    real(kind=8) :: toutms(*)
!
    do i = 1, nvi
        ee(i) = 0.d0
        f(i) = 0.d0
    end do
!
    if (kp .eq. 1) then
!        INTEGRATION Y1=F(Y,T)
        call rdif01(fami, kpg, ksp, rela_comp, mod, &
                    imat, matcst, nbcomm, cpmono, nfs, &
                    nsg, toutms, nvi, nmat, y, &
                    cothe, coeff, dcothe, dcoeff, pgl, &
                    nbphas, coel, x, pas, neps, &
                    epsd, detot, f, nhsr, numhsr, &
                    hsr, itmax, toler, iret)
        do i = 1, nvi
            a(i) = f(i)
            y(i) = y(i)+a(i)*h
        end do
    else
        do i = 1, nvi
            y(i) = y(i)+a(i)*h
        end do
    end if
!
    x = x+h
!     INTEGRATION Y2=F(Y1,T+H)
!
    call rdif01(fami, kpg, ksp, rela_comp, mod, &
                imat, matcst, nbcomm, cpmono, nfs, &
                nsg, toutms, nvi, nmat, y, &
                cothe, coeff, dcothe, dcoeff, pgl, &
                nbphas, coel, x, pas, neps, &
                epsd, detot, f, nhsr, numhsr, &
                hsr, itmax, toler, iret)
!
    hs2 = 0.5d0*h
!
    do i = 1, nvi
        ee(i) = (f(i)-a(i))*hs2
        y(i) = y(i)+ee(i)
    end do
!
end subroutine
