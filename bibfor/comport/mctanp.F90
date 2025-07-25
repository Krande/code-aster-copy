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

subroutine mctanp(dpstrs, rprops, ii, jj, mm, edge, right, apex)
!***********************************************************************
!
! OBJECT:
!
! COMPUTATION OF CONSISTENT TANGENT MODULUS FOR MOHR-COULOMB TYPE
! ELASTO-PLASTIC MATERIAL WITH ASSOCIATIVE/NON-ASSOCIATIVE FLOW RULE
!
! ----------------------------------------------------------------------
!
!     LOI DE COMPORTEMENT DE MOHR-COULOMB
!
! IN  PSTRS   : CONTRAINTES PRINCIPALES
! IN  II      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MINEURE
! IN  JJ      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE MAJEURE
! IN  MM      : COMPOSANTE DE LA CONTRAINTE PRINCIPALE INTERMEDIAIRE
! IN  EDGE    : Y-A-T-IL DEUX  MECANISMES ACTIFS
! IN  RIGHT   : SI EDGE, LA PROJECTION S'EFFECTUE-T-IL A DROITE?
! IN  APEX    : Y-A-T-IL TROIS MECANISMES ACTIFS
!
! OUT DPSTRS  : MATRICE TANGENTE DANS L'ESPACE DES DIRECTIONS PRINCIPALE
!
!***********************************************************************
    implicit none
! ======================================================================
!
    real(kind=8) :: dpstrs(3, 3)
    real(kind=8) :: rprops(*)
    real(kind=8) :: edge
    real(kind=8) :: right
    real(kind=8) :: apex
!
#include "asterf_types.h"
!
! Declaration of integer type variables
    integer(kind=8) :: ii, jj, mm
    real(kind=8) :: degr
!
!     aster_logical :: epflag, epflag0
!
    parameter(degr=0.017453292519943295d0)
! Real arrays and variables
    real(kind=8) :: young, poiss, sinphi, sinpsi, constb, gmodu, bulk
    real(kind=8) :: r2g, r4g, r2bulk, r1d3, r2d3, r2gd3, r4gd3
    real(kind=8) :: sphsps, consta
    real(kind=8) :: denom, b1, b2, b3, drvaa, drvab, drvba, drvbb, aux1
    real(kind=8) :: aux2, aux3, r1ddet, r1, r2, r3, r4, r0
!
    data r0, r1, r2, r3, r4/&
     &    0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0/
!
!
! Set material properties
    young = rprops(2)
    poiss = rprops(3)
    sinphi = sin(degr*rprops(4))
    sinpsi = sin(degr*rprops(5))
!
! Set some constants
    gmodu = young/(r2*(r1+poiss))
    bulk = young/(r3*(r1-r2*poiss))
    r2g = r2*gmodu
    r4g = r4*gmodu
    r2bulk = r2*bulk
    r1d3 = r1/r3
    r2d3 = r2*r1d3
    r2gd3 = r2g*r1d3
    r4gd3 = r4g*r1d3
!
! Compute elastoplastic consistent tangent
! ----------------------------------------
    if (edge .eq. r1) then
! -------------------------------------------------------------------------
!
! Tangent consistent with 2-vector return to edge
!
! -------------------------------------------------------------------------
        sphsps = sinphi*sinpsi
        consta = r4g*(r1+r1d3*sphsps)+r4*bulk*sphsps
        if (right .eq. r1) then
            constb = r2g*(r1+sinphi+sinpsi-r1d3*sphsps)+r4*bulk*sphsps
        else
            constb = r2g*(r1-sinphi-sinpsi-r1d3*sphsps)+r4*bulk*sphsps
        end if
        drvaa = -consta
        drvab = -constb
        drvba = -constb
        drvbb = -consta
        aux1 = r2g*(r1+r1d3*sinpsi)+r2bulk*sinpsi
        aux2 = (r4gd3-r2bulk)*sinpsi
        aux3 = r2g*(r1-r1d3*sinpsi)-r2bulk*sinpsi
        r1ddet = r1/(drvaa*drvbb-drvab*drvba)
!
        if (right .eq. r1) then
!
! ...returned to right edge
! -------------------------
            dpstrs(ii, ii) = bulk+r4gd3+aux1*(-drvab+drvbb+drvaa-drvba)* &
                             (r2g+(r2bulk+r2gd3)*sinphi)*r1ddet
!
            dpstrs(ii, mm) = bulk-r2gd3+aux1*(r2g*(drvab-drvaa)+((- &
                                            drvab+drvbb+drvaa-drvba)*(r2bulk+r2gd3)+(drvba-drvbb)* &
                                                                 r2g)*sinphi)*r1ddet
!
            dpstrs(ii, jj) = bulk-r2gd3+aux1*(r2g*(drvba-drvbb)+((- &
                                            drvab+drvbb+drvaa-drvba)*(r2bulk+r2gd3)+(drvab-drvaa)* &
                                                                 r2g)*sinphi)*r1ddet
!
            dpstrs(mm, ii) = bulk-r2gd3+(aux2*(drvab-drvbb)+aux3*(drvba- &
                                                          drvaa))*(r2g+(r2bulk+r2gd3)*sinphi)*r1ddet
!
            dpstrs(mm, mm) = bulk+r4gd3+(aux2*((r2bulk*(drvab-drvbb)+ &
                                         (drvab*r2gd3+drvbb*r4gd3))*sinphi-drvab*r2g)+aux3*(drvaa* &
                                             r2g+(r2bulk*(drvba-drvaa)-(drvaa*r2gd3+drvba*r4gd3))* &
                                                                                     sinphi))*r1ddet
!
            dpstrs(mm, jj) = bulk-r2gd3+(aux2*((r2bulk*(drvab-drvbb)- &
                                              (drvbb*r2gd3+drvab*r4gd3))*sinphi+drvbb*r2g)+aux3*(( &
                                           r2bulk*(drvba-drvaa)+(drvaa*r4gd3+drvba*r2gd3))*sinphi- &
                                                                                  drvba*r2g))*r1ddet
!
            dpstrs(jj, ii) = bulk-r2gd3+((aux2*(drvba-drvaa)+aux3*(drvab- &
                                                         drvbb))*((r2bulk+r2gd3)*sinphi+r2g))*r1ddet
!
            dpstrs(jj, mm) = bulk-r2gd3+(aux2*(((r2bulk*(drvba-drvaa)- &
                                            (drvba*r4gd3+drvaa*r2gd3))*sinphi)+drvaa*r2g)+aux3*((( &
                                          r2bulk*(drvab-drvbb)+(drvab*r2gd3+drvbb*r4gd3))*sinphi)- &
                                                                                  drvab*r2g))*r1ddet
!
            dpstrs(jj, jj) = bulk+r4gd3+(aux2*(((r2bulk*(drvba-drvaa)+ &
                                            (drvaa*r4gd3+drvba*r2gd3))*sinphi)-drvba*r2g)+aux3*((( &
                                          r2bulk*(drvab-drvbb)-(drvab*r4gd3+drvbb*r2gd3))*sinphi)+ &
                                                                                  drvbb*r2g))*r1ddet
        else
!
! ...returned to left edge
! -------------------------
            dpstrs(ii, ii) = bulk+r4gd3+(aux1*(((r2bulk*(drvbb-drvab)+ &
                                            (drvab*r4gd3+drvbb*r2gd3))*sinphi)+drvbb*r2g)+aux2*((( &
                                          r2bulk*(drvba-drvaa)+(drvaa*r4gd3+drvba*r2gd3))*sinphi)+ &
                                                                                  drvba*r2g))*r1ddet
!
            dpstrs(ii, mm) = bulk-r2gd3+(aux1*(((r2bulk*(drvbb-drvab)- &
                                            (drvab*r2gd3+drvbb*r4gd3))*sinphi)-drvab*r2g)+aux2*((( &
                                          r2bulk*(drvba-drvaa)-(drvaa*r2gd3+drvba*r4gd3))*sinphi)- &
                                                                                  drvaa*r2g))*r1ddet
!
            dpstrs(ii, jj) = bulk-r2gd3+((aux1*(drvbb-drvab)+aux2*(drvba- &
                                                       drvaa))*(((r2bulk+r2gd3)*sinphi)-r2g))*r1ddet
!
            dpstrs(mm, ii) = bulk-r2gd3+(aux1*(((r2bulk*(drvaa-drvba)- &
                                            (drvaa*r4gd3+drvba*r2gd3))*sinphi)-drvba*r2g)+aux2*((( &
                                          r2bulk*(drvab-drvbb)-(drvab*r4gd3+drvbb*r2gd3))*sinphi)- &
                                                                                  drvbb*r2g))*r1ddet
!
            dpstrs(mm, mm) = bulk+r4gd3+(aux1*(((r2bulk*(drvaa-drvba)+ &
                                            (drvaa*r2gd3+drvba*r4gd3))*sinphi)+drvaa*r2g)+aux2*((( &
                                          r2bulk*(drvab-drvbb)+(drvab*r2gd3+drvbb*r4gd3))*sinphi)+ &
                                                                                  drvab*r2g))*r1ddet
!
            dpstrs(mm, jj) = bulk-r2gd3+((aux1*(drvaa-drvba)+aux2*(drvab- &
                                                       drvbb))*(((r2bulk+r2gd3)*sinphi)-r2g))*r1ddet
!
            dpstrs(jj, ii) = bulk-r2gd3+(aux3*(((r2bulk*(drvab-drvbb- &
                                            drvaa+drvba)+(drvaa-drvab)*r4gd3+(drvba-drvbb)*r2gd3)* &
                                                sinphi)+(drvba-drvbb)*r2g))*r1ddet
!
            dpstrs(jj, mm) = bulk-r2gd3+(aux3*(((r2bulk*(drvab-drvbb- &
                                            drvaa+drvba)+(drvab-drvaa)*r2gd3+(drvbb-drvba)*r4gd3)* &
                                                sinphi)+(drvab-drvaa)*r2g))*r1ddet
!
            dpstrs(jj, jj) = bulk+r4gd3+(aux3*(drvab-drvbb-drvaa+drvba)* &
                                         (((r2bulk+r2gd3)*sinphi)-r2g))*r1ddet
        end if
    else if (apex .eq. r1) then
! -------------------------------------------------------------------------
!
! Tangent consistent with multi-vector return to apex
!
! -------------------------------------------------------------------------
        dpstrs(:, :) = r0
    else
! -------------------------------------------------------------------------
!
! Tangent consistent with 1-vector return to main active plane
!
! -------------------------------------------------------------------------
        sphsps = sinphi*sinpsi
        consta = r4g*(r1+r1d3*sphsps)+r4*bulk*sphsps
        denom = -consta
!
        b1 = (r2g*(r1+r1d3*sinpsi)+r2bulk*sinpsi)/denom
        b2 = (r4g*r1d3-r2bulk)*sinpsi/denom
        b3 = (r2g*(r1-r1d3*sinpsi)-r2bulk*sinpsi)/denom
!
        dpstrs(ii, ii) = r2g*(r2d3+b1*(r1+r1d3*sinphi))+bulk*(r1+r2*b1* &
                                                              sinphi)
!
        dpstrs(ii, mm) = r1d3*(r3*bulk-r2g)*(r1+r2*b1*sinphi)
!
        dpstrs(ii, jj) = r2g*(-r1d3-b1*(r1-r1d3*sinphi))+bulk*(r1+r2*b1* &
                                                               sinphi)
!
        dpstrs(mm, ii) = r2g*(-r1d3-b2*(r1+r1d3*sinphi))+bulk*(r1-r2*b2* &
                                                               sinphi)
!
        dpstrs(mm, mm) = r4g*r1d3*(r1+b2*sinphi)+bulk*(r1-r2*b2*sinphi)
!
        dpstrs(mm, jj) = r2g*(-r1d3+b2*(r1-r1d3*sinphi))+bulk*(r1-r2*b2* &
                                                               sinphi)
!
        dpstrs(jj, ii) = r2g*(-r1d3-b3*(r1+r1d3*sinphi))+bulk*(r1-r2*b3* &
                                                               sinphi)
!
        dpstrs(jj, mm) = r1d3*(r3*bulk-r2g)*(r1-r2*b3*sinphi)
!
        dpstrs(jj, jj) = r2g*(r2d3+b3*(r1-r1d3*sinphi))+bulk*(r1-r2*b3* &
                                                              sinphi)
    end if
!
end subroutine
