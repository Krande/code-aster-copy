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
subroutine lkdgde(val, vintr, dt, seuive, ucrim, &
                  im, sm, vinm, nbmat, mater, &
                  depsv, dgamv, retcom)
!
    implicit none
#include "asterfort/lcdevi.h"
#include "asterfort/lkbpri.h"
#include "asterfort/lkcalg.h"
#include "asterfort/lkcaln.h"
#include "asterfort/lkdfds.h"
#include "asterfort/lkdhds.h"
#include "asterfort/lkds2h.h"
#include "asterfort/lkvacv.h"
#include "asterfort/lkvarv.h"
    integer(kind=8) :: nbmat, retcom, val
    real(kind=8) :: seuive, ucrim, im, sm(6), vintr
    real(kind=8) :: mater(nbmat, 2), vinm(7), depsv(6), dgamv
    real(kind=8) :: dt
! --- MODELE LETK : LAIGLE VISCOPLASTIQUE--------------------------
! =================================================================
! --- BUT : DEFINITION DE LA DEFORMATION VISQUEUSE ET DU PARAMETRE
! ---- D ECROUISSAGE VISQUEUX
! =================================================================
! IN  : VAL   :  INDICATEUR POUR LES LOIS DE DILATANCE ------------
! --- : VINTR  :  INDICATEUR CONTRACTANCE OU  DILATANCE ------------
! --- : DT    :  PAS DE TEMPS -------------------------------------
! --- : SEUIVE:  SEUIL VISQUEUX EN FONCTION DE LA PREDICITION------
!---- : UCRIM :  EN FONCTION DES CONTRAINTES A LINSTANT MOINS------
! --- : IM    :  INVARIANT DES CONTRAINTES A L INSTANT MOINS-------
! --- : SM    :  DEVIATEUR DES CONTRAINTES A L INSTANT MOINS-------
! --- : VINM  :  VARIABLES INTERNES -------------------------------
! --- : NBMAT :  NOMBRE DE PARAMETRES MATERIAU --------------------
! --- : MATER :  COEFFICIENTS MATERIAU A T+DT ---------------------
! ----------- :  MATER(*,1) = CARACTERISTIQUES ELASTIQUES ---------
! ----------- :  MATER(*,2) = CARACTERISTIQUES PLASTIQUES ---------
! OUT : DEPSV : DEFORMATIONS VISQUEUSES ---------------------------
!     : DGAMV : PARAMETRE D ECROUISSAGE VISQUEUX ------------------
! --- : RETCOM: CODE RETOUR POUR REDECOUPAGE DU PAS DE TEMPS-------
! =================================================================
    common/tdim/ndt, ndi
    integer(kind=8) :: i, ndi, ndt
    real(kind=8) :: a, n, pa
    real(kind=8) :: bidon, deux, trois, zero
    real(kind=8) :: paravi(3), varvi(4)
    real(kind=8) :: dhds(6), ds2hds(6), dfdsv(6)
    real(kind=8) :: bprime, vecnv(6), gv(6)
    real(kind=8) :: ddepsv(6)
! =================================================================
! --- INITIALISATION DE PARAMETRES --------------------------------
! =================================================================
    parameter(zero=0.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
    dfdsv = 0.d0
    vecnv = 0.d0
! =================================================================
! --- RECUPERATION DES DONNEES MATERIAUX --------------------------
! =================================================================
    pa = mater(1, 2)
    a = mater(21, 2)
    n = mater(22, 2)
! =================================================================
! --- CALCUL DE DF/DSIG ------------------------------------
! =================================================================
!
    call lkdhds(nbmat, mater, im, sm, dhds, &
                retcom)
    call lkds2h(nbmat, mater, im, sm, dhds, &
                ds2hds, retcom)
    call lkvarv(vintr, nbmat, mater, paravi)
!
    call lkvacv(nbmat, mater, paravi, varvi)
    call lkdfds(nbmat, mater, sm, paravi, varvi, &
                ds2hds, ucrim, dfdsv)
!
    bprime = lkbpri(val, vinm, nbmat, mater, paravi, im, sm)
!
    call lkcaln(sm, bprime, vecnv, retcom)
! =================================================================
! --- CALCUL DE GVISC ------------------------------------
! =================================================================
    call lkcalg(dfdsv, vecnv, gv, bidon)
! =================================================================
! --- CALCUL DE DEPSV ------------------------------------
! =================================================================
    depsv = 0.d0
    do i = 1, ndt
        if (seuive .le. zero) then
            depsv(i) = zero
        else
            depsv(i) = a*(seuive/pa)**n*gv(i)*dt
        end if
    end do
!
! =================================================================
! --- CALCUL DU DEVIATEUR DU TENSEUR DES DEFORMATIONS VISQUEUSES -
! =================================================================
    call lcdevi(depsv, ddepsv)
!
! =================================================================
! --- CALCUL DE DGAMV ------------------------------------
! =================================================================
!
    dgamv = 0.d0
!
    do i = 1, ndt
        dgamv = dgamv+ddepsv(i)**2
    end do
    dgamv = sqrt(deux/trois*dgamv)
! =================================================================
end subroutine
