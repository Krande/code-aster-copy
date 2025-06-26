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
subroutine lcmmin(typess, essai, mod, nmat, materf, &
                  nr, nvi, yd, deps, dy, &
                  comp, nbcomm, cpmono, pgl, nfs, &
                  nsg, toutms, timed, timef, vind, &
                  sigd, epstr)
! aslint: disable=W1504
    implicit none
! person_in_charge: jean-michel.proix at edf.fr
!       MONOCRISTAL : CALCUL SOLUTION INITIALE
!
!       IN  ESSAI  :  VALEUR DE LA SOLUTION D ESSAI
!           MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           MATERF :  COEFFICIENTS MATERIAU A T+DT
!           NR     :  DIMENSION DECLAREE DRDY
!           NVI    :  NOMBRE DE VARIABLES INTERNES
!           YD     :  VARIABLES A T
!           DY     :  SOLUTION  A L'ITERATION PRECEDENTE
!           COMP   :  NOM COMPORTEMENT
!           NBCOMM :  INCIDES DES COEF MATERIAU
!           CPMONO :  NOM DES COMPORTEMENTS
!           PGL    :  MATRICE DE PASSAGE
!           TOUTMS :  TENSEURS D'ORIENTATION
!           VIND   :  VARIABLES INTERNES A L'INSTANT PRECEDENT
!           SIGD   :  CONTRAINTE A T
!       VAR DEPS   :  INCREMENT DE DEFORMATION (ACTUALISE EN C_PLAN)
!           TYPESS :  TYPE DE SOLUTION D ESSAI
!                               0 = NUL(0)
!                               1 = ELASTIQUE
!                               2 = EXPLICITE (=-1 INITIALEMENT)
!                               3 = ESSAI
!       OUT DY     :  SOLUTION ESSAI  = ( DSIG DVIN (DEPS3) )
!       ----------------------------------------------------------------
!
#include "asterfort/lcmmsg.h"
#include "asterfort/lcopil.h"
#include "asterfort/lcopli.h"
#include "asterfort/r8inir.h"
#include "asterfort/tnsvec.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
    integer(kind=8) :: ndt, ndi, typess, nmat, nr, nvi, types0, nfs, nsg
!
    real(kind=8) :: yd(nr), dy(nr), essai
    real(kind=8) :: hook(6, 6)
    real(kind=8) :: deps(6)
    real(kind=8) :: dsig(6)
    real(kind=8) :: epstr(6), dkooh(6, 6), epsed(6)
!
    real(kind=8) :: materf(nmat, 2)
    real(kind=8) :: toutms(nfs, nsg, 6)
!
    character(len=8) :: mod
!     ----------------------------------------------------------------
    common/tdim/ndt, ndi
!     ----------------------------------------------------------------
    integer(kind=8) :: i, nbfsys, nbsys, is, nbcomm(nmat, 3), ifa, nums
    real(kind=8) :: evp(6), fe(3, 3), df(3, 3), fe1(3, 3), fe1t(3, 3)
    real(kind=8) :: pgl(3, 3), ms(6), ng(3), q(3, 3), lg(3), fetfe(3, 3)
    character(len=16) :: comp(*)
    character(len=24) :: cpmono(5*nmat+1)
    character(len=24) :: nomfam
    real(kind=8) :: timed, timef, vind(*), sigd(6), sigdn(6)
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
!
! - SOLUTION INITIALE = NUL
!
    types0 = typess
!
    typess = 0
!         TYPESS=7
!
!     POUR LE CRITERE DE CONVERGENCE CF LCMMCV
    if (materf(nmat, 1) .eq. 0) then
        call lcopil('ISOTROPE', mod, materf(1, 1), dkooh)
    else if (materf(nmat, 1) .eq. 1) then
        call lcopil('ORTHOTRO', mod, materf(1, 1), dkooh)
    end if
    epsed(1:ndt) = matmul(dkooh(1:ndt, 1:ndt), sigd(1:ndt))
    epstr(1:ndt) = epsed(1:ndt)+deps(1:ndt)
!
    if (typess .eq. 0) then
        dy(:) = 0.d0
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = 0.d0
        end if
! Les autres intitialisations ci-dessous ne sont pas utilisées
! actuellement pour la loi MONOCRISTAL
!
! - SOLUTION INITIALE = ELASTIQUE
!
    else if (typess .eq. 1 .or. typess .eq. -1) then
        if (materf(nmat, 1) .eq. 0) then
            call lcopli('ISOTROPE', mod, materf(1, 1), hook)
        else if (materf(nmat, 1) .eq. 1) then
            call lcopli('ORTHOTRO', mod, materf(1, 1), hook)
        end if
!        GDEF : INITIALISATION PAR FET.DFT.DF.FE
        if (gdef .eq. 1) then
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vind(nvi-3-18+10), b_incx, fe, b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, deps, b_incx, df, b_incy)
            fe1 = matmul(df, fe)
            fe1t = transpose(fe1)
            fetfe = matmul(fe1t, fe1)
            call tnsvec(3, 3, fetfe, dy, 1.d0)
        else
            hook(1:ndt, 1:ndt) = transpose(hook(1:ndt, 1:ndt))
            dsig(1:ndt) = matmul(hook(1:ndt, 1:ndt), deps(1:ndt))
            dy(1:ndt) = dsig(1:ndt)
        end if
!
! - SOLUTION INITIALE = EXPLICITE
!
!      ELSEIF ( TYPESS .EQ. 2 ) THEN
!
! - SOLUTION INITIALE = VALEUR ESSAI POUR TOUTES LES COMPOSANTES
!
    else if (typess .eq. 3) then
        dy(:) = essai
        if (mod(1:6) .eq. 'C_PLAN') then
            deps(3) = essai
            dy(3) = 0.d0
        end if
!
    else if (typess .eq. 7) then
!
        nbfsys = nbcomm(nmat, 2)
        nums = 0
!
        do ifa = 1, nbfsys
!
            nomfam = cpmono(5*(ifa-1)+1)
!       RECUPERATION DU NOMBRE DE SYSTEME DE GLISSEMENT NBSYS
            call lcmmsg(nomfam, nbsys, 0, pgl, ms, &
                        ng, lg, 0, q)
            if (nbsys .eq. 0) then
                call utmess('F', 'ALGORITH_70')
            end if
!
            call r8inir(6, 0.d0, evp, 1)
!
            do is = 1, nbsys
                nums = nums+1
                dy(ndt+6+3*ifa*(is-1)+3) = vind(6+3*ifa*(is-1)+3)*(timef-timed)/timef
                dy(ndt+6+3*ifa*(is-1)+2) = abs(vind(6+3*ifa*(is-1)+2))*(timef-timed)/timef
                dy(ndt+6+3*ifa*(is-1)+1) = vind(6+3*ifa*(is-1)+1)*(timef-timed)/timef
!           RECUPERATION DE MS ET CALCUL DE EVP
                call lcmmsg(nomfam, nbsys, is, pgl, ms, &
                            ng, lg, 0, q)
                do i = 1, 6
                    evp(i) = evp(i)+ms(i)*dy(ndt+6+3*ifa*(is-1)+2)
                end do
            end do
        end do
!      ATTRIBUTIION A DY LA VALEUR DE EVP CALCULEE
        dy(ndt+1:ndt+6) = evp(1:6)
!
        do i = 1, 6
            sigdn(i) = sigd(i)*(timef-timed)/timef
        end do
        dy(1:ndt) = sigdn(1:ndt)
!
!
        if (mod(1:6) .eq. 'C_PLAN') then
            dy(3) = 0.d0
        end if
!
    end if
!
    typess = types0
end subroutine
