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

subroutine cadldp(vp, sigeqe, nbmat, materf, parame, &
                  derive, sig3, eta, dg, detadg, &
                  dgdl, ddldsp)
    implicit none
    integer(kind=8) :: nbmat
    real(kind=8) :: ddldsp, materf(nbmat, 2), parame(4), derive(5)
    real(kind=8) :: vp(3), sigeqe, eta, dg, detadg, sig3, dgdl
! ======================================================================
! --- LOI DE HOEK BROWN : CALCUL DE DDLAMBDA/DSIP (CONT. TOTALES) -----
! ======================================================================
! IN  VP      VALEURS PROPRES DU DEVIATEUR ELASTIQUE SE ----------------
! IN  NBMAT   NOMBRE DE DONNEES MATERIAU -------------------------------
! IN  MATERF  DONNEES MATERIAU -----------------------------------------
! IN  PARAME  VALEUR DES PARAMETRES DE LA LOI S*SIG, M*SIG, B ----------
! IN  DERIVE  VALEUR DES DERIVEES DES PARAMETRES PAR RAPPORT A GAMMA ---
! IN  SIG3    CONTRAINTE PRINCIPALE SIG3 -------------------------------
! IN  DG      INCREMENT DU PARAMETRE D ECROUISSAGE GAMMA ---------------
! IN  DETADG  DERIVEE DE ETA PAR RAPPORT A GAMMA -----------------------
! IN  DGDL    DERIVEE  DE GAMMA PAR RAPPORT A LAMBDA -------------------
! OUT DDLDSP  DDLAMDA/DSIP ---------------------------------------------
! ======================================================================
    real(kind=8) :: un, deux, trois, k, dl, mu
    real(kind=8) :: a2, a3, a4, c5, a6, aux1, aux2, aux3, denom
    integer(kind=8) :: ndt, ndi
! ======================================================================
    parameter(un=1.0d0)
    parameter(deux=2.0d0)
    parameter(trois=3.0d0)
! ======================================================================
    common/tdim/ndt, ndi
! ======================================================================
! --- INITIALISATIONS --------------------------------------------------
! ======================================================================
    mu = materf(4, 1)
    k = materf(5, 1)
! ======================================================================
    a2 = vp(3)-vp(1)
    a3 = trois*mu/sigeqe
    a4 = trois*k*eta
    c5 = un/materf(14, 2)
    a6 = a3*vp(3)
    dl = dg/(eta+un)
! ======================================================================
! --- CALCUL DU DENOMINATEUR -------------------------------------------
! ======================================================================
    aux1 = parame(1)-parame(2)*sig3
    aux2 = dl*trois*k*detadg*dgdl+a6+a4
    aux3 = dgdl*(derive(1)-sig3*derive(2))+parame(2)*aux2
    denom = -a2*a3-derive(3)*dgdl*(un+c5*sig3)+parame(3)*c5*aux2-aux3/(sqrt(aux1)*deux)
! ======================================================================
! --- CALCUL DE DDL/DSIGP ----------------------------------------------
! ======================================================================
    ddldsp = (parame(3)*c5-parame(2)*0.5d0/sqrt(aux1))/denom
! ======================================================================
end subroutine
