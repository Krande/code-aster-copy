! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine mdflam(dnorm, vitloc, knorm, cnorm, cost, sint, &
                  flim, fseuil, rigifl, defpla, fnorma, &
                  flocal, vnorm, defmax, enfo_fl, def, &
                  deft0, deft, amor, cfl, critamor)
    implicit none
!

#include "asterfort/utmess.h"

!***********************************************************************
! 01/01/91    G.JACQUART AMV/P61 47 65 49 41
!***********************************************************************
!     FONCTION  : CALCULE LA DISTANCE NORMALE A L'OBSTACLE (<0 SI CHOC)
!
!-----------------------------------------------------------------------
!                             ARGUMENTS
! .________________.____.______________________________________________.
!    DNORM          <--   DISTANCE NORMALE A L'OBSTACLE
!    VITLOC         <--   VITESSE DANS LE REPERE LOCAL
!    COST,SINT      <--   DIRECTION NORMALE A L'OBSTACLE
!    KNORM          <--   RAIDEUR NORMALE DE CHOC
!    CNORM          <--   AMORTISSEMENT NORMAL DE CHOC
!    FLIM           <--   EFFORT MAXIMAL DE CHOC
!    FSEUIL         <--   EFFORT MAXIMAL DE CHOC POST FLAMBAGE
!    RIGIFL         <--   RAIDEUR NORMALE DE CHOC POST FLAMBAGE
!    CFL            <--   AMORTISSEMENT NORMAL POST FLAMBAGE
!    DEFPLA         <--   DEFORMATION PLASTIQUE
!    DEFMAX         <--   DEFORMATION TOTALE MAXIMALE
!    ENFO_FL        <--   ENFONCEMENT AU FLAMBAGE
!    DEF            <--   LISTE DE DEFORMATIONS PLASTIQUES POST FLAMBAGE
!    DEFT0          <--   DEFORMATION TOTALE A LA FIN DU PLATEAU
!    DEFT           <--   LISTE DE DEFORMATIONS TOTALES POST FLAMBAGE
!    AMOR           <--   LISTE DES AMORTISSEMENTS POST FLAMABGE
!    CRITAMOR       <--   0  ou 1 Amortissement inclus ou exclus au critere
!    FNORMA          -->  FORCE NORMALE DE CHOC  (MODULE)
!    FLOCAL          -->  FORCE NORMALE DE CHOC REP. LOCAL
!-----------------------------------------------------------------------
    real(kind=8) :: vitloc(3), flocal(3), knorm, fnorma
!-----------------------------------------------------------------------
    real(kind=8) :: cost, defpla, dnorm, flim, fseuil, rigifl, sint
    real(kind=8) :: vnorm, enfo_fl, defmax, cnorm, deft0, cfl
    real(kind=8) :: alpha, cfl2, cnorm2
    integer :: j, critamor

    real(kind=8), pointer  :: def(:)
    real(kind=8), pointer  :: deft(:)
    real(kind=8), pointer  :: amor(:)

!-----------------------------------------------------------------------
    vnorm = vitloc(2)*cost+vitloc(3)*sint

    if (critamor .eq. 0) then
! --- Amortissement exclus au critere ---
        cnorm2 = 0.0d0
        cfl2 = 0.0d0
    else if (critamor .eq. 1) then
! --- Amortissement inclus au critere ---
        cnorm2 = cnorm
        cfl2 = cfl
    end if

    if (defpla .le. 0.d0) then
!     --- FLAMBAGE NON ENCORE RENCONTRE ---
        if (-dnorm .lt. 0.d0) then
            fnorma = 0.0d0
            rigifl = knorm
            cfl = cnorm
        else
            if (-dnorm .lt. (flim+cnorm2*vnorm)/knorm) then
                fnorma = -knorm*dnorm-cnorm2*vnorm
                rigifl = knorm
                cfl = cnorm
                if (fnorma .lt. 0.d0) fnorma = 0.d0
            else
!           --- DEBUT DU FLAMBAGE ---
                fnorma = flim
                defpla = 1.d-20
                rigifl = knorm
                cfl = cnorm
            end if
        end if
    else
!     --- LE FLAMBAGE A DEJA EU LIEU --- vnorm - a la charge puis + a la decharge
!     --- Si charge inferieure au jeu mis a jour
        if (-dnorm .lt. defpla) then
            fnorma = 0.0d0
        else
!     --- Si decharge ou charge inferieure a la limite
            if (vnorm .gt. 0.d0 .or. -dnorm .le. defmax) then
                fnorma = -rigifl*(dnorm+defpla)-cfl2*vnorm
                if (critamor .eq. 1) then
                    if ((-dnorm .lt. deft0) .and. (fnorma .ge. flim)) then
                        fnorma = flim
                    else if ((fnorma .ge. flim+((fseuil-flim)/enfo_fl)*(-dnorm-deft0)) &
                             .and. (-dnorm .lt. deft(1))) then
                        fnorma = flim+((fseuil-flim)/enfo_fl)*(-dnorm-deft0)
                    else if ((fnorma .ge. fseuil) .and. (-dnorm .ge. deft(1))) then
                        fnorma = fseuil
                    end if
                    if (fnorma .lt. 0.d0) fnorma = 0.d0
                end if
            else
!     --- Deformation pendant le plateau
                if (-dnorm .lt. deft0) then
                    fnorma = flim
                    defpla = -dnorm-flim/rigifl
                end if
!     --- Deformation pendant le flambage
                if ((-dnorm .ge. deft0) .and. (-dnorm .lt. deft(1))) then
                    fnorma = flim-((fseuil-flim)/enfo_fl)*(dnorm+deft0)
                    rigifl = fnorma/(-dnorm-defpla)
                    defpla = def(1)
                    cfl = cnorm-((amor(1)-cnorm)/enfo_fl)*(dnorm+deft0)
                end if
!     --- Deformation post flambage
                if (-dnorm .ge. deft(1) .and. (size(def) .lt. 2)) then
                    fnorma = fseuil
                    rigifl = fseuil/(deft(1)-def(1))
                    defpla = -dnorm-fseuil/rigifl
                    cfl = amor(1)
                else
                    do j = 1, (size(def)-1)
                        if (-dnorm .ge. deft(j) .and. -dnorm .lt. deft(j+1)) then
                            fnorma = fseuil
                            alpha = (-dnorm-deft(j))/(deft(j+1)-deft(j))
                            defpla = def(j)+alpha*(def(j+1)-def(j))
                            cfl = amor(j)+alpha*(amor(j+1)-amor(j))
                            rigifl = fseuil/(-dnorm-defpla)
                        end if
                    end do
                    if (-dnorm .gt. deft(size(deft))) then
                        fnorma = fseuil
                        rigifl = fseuil/(deft(size(deft))-def(size(def)))
                        defpla = -dnorm-fseuil/rigifl
                        cfl = amor(size(amor))
                        call utmess('A', 'ALGORITH5_86')
                    end if
                end if
            end if
        end if
    end if

    if (defpla .lt. 0.d0) defpla = 0.d0

    if (-dnorm .gt. defmax) defmax = -dnorm

    if (critamor .eq. 0) then
        fnorma = fnorma-cfl*vnorm
        if (fnorma .lt. 0.d0) fnorma = 0.d0
        if (-dnorm .lt. defpla) fnorma = 0.0d0
    end if

    flocal(1) = 0.d0
    flocal(2) = fnorma*cost
    flocal(3) = fnorma*sint
end subroutine
