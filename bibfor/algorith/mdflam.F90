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

subroutine mdflam(dnorm, vitloc, cost, sint, &
                  critfl, rigifl, amorfl, defpla, &
                  fnorma, flocal, vnorm, critamor)
    implicit none
!

!***********************************************************************
! 01/01/91    G.JACQUART AMV/P61 47 65 49 41
!***********************************************************************
!     FONCTION  : CALCULE L'EFFORT DE FLAMBAGE ASSOCIE A UN OBSTACLE
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
!    DNORM          <--   DISTANCE NORMALE A L'OBSTACLE
!    VNORM          <--   VITESSE NORMALE A L'OBSTACLE
!    VITLOC         <--   VITESSE DANS LE REPERE LOCAL
!    COST,SINT      <--   DIRECTION NORMALE A L'OBSTACLE
!    CRITFL         <--   EFFORT SEUIL DE FLAMBAGE
!    RIGIFL         <--   RAIDEUR NORMALE DE CHOC/FLAMBAGE
!    AMORFL         <--   AMORTISSEMENT NORMAL DE CHOC/FLAMBAGE
!    DEFPLA         <--   DEFORMATION PLASTIQUE
!    CRITAMOR       <--   Amortissement inclus (1) ou exclus (0) au critere
!    FNORMA          -->  FORCE NORMALE DE CHOC  (MODULE)
!    FLOCAL          -->  FORCE NORMALE DE CHOC REP. LOCAL
!-----------------------------------------------------------------------
    real(kind=8) :: vitloc(3), flocal(3), fnorma
!-----------------------------------------------------------------------
    real(kind=8) :: cost, sint, dnorm, vnorm
    real(kind=8) :: defpla, critfl, rigifl, amorfl
    real(kind=8) :: amorfl2
    integer(kind=8) :: critamor

!-----------------------------------------------------------------------
! --- Calcul de la vitesse normale
    vnorm = vitloc(2)*cost+vitloc(3)*sint

! --- Calcul de l'effort normal et de la déformation plastique
    if (critamor .eq. 0) then
        ! --- Amortissement exclus au critere ---
        amorfl2 = 0.0d0
    else if (critamor .eq. 1) then
        ! --- Amortissement inclus au critere ---
        amorfl2 = amorfl
    end if

    fnorma = -rigifl*(dnorm+defpla)-amorfl2*vnorm
    if (fnorma .ge. critfl) then
        fnorma = critfl
        defpla = -dnorm-critfl/rigifl
    end if

! --- Mise à jour de la déformation plastique ---
    if (defpla .lt. 0.d0) defpla = 0.d0

! --- Mise à jour de l'effort normal ---
    if (critamor .eq. 0) then
        fnorma = fnorma-amorfl*vnorm
    end if

    if (fnorma .lt. 0.d0) fnorma = 0.d0
    if (-dnorm .lt. defpla) fnorma = 0.0d0

! --- Calcul de l'effort de flambage ---
    flocal(1) = 0.d0
    flocal(2) = fnorma*cost
    flocal(3) = fnorma*sint
end subroutine
