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

subroutine pmfitg(typfib, nf, ncarf, vf, vs)
!
!
! --------------------------------------------------------------------------------------------------
!
!                           Intégrations sur la section des PMF
!
! person_in_charge: jean-luc.flejou at edf.fr
! --------------------------------------------------------------------------------------------------
!
! IN
!   typfib  : type des fibres : 1 ou 2
!   nf      : nombre de fibres
!   ncarf   : nombre de caractéristiques par fibre
!           !!! Quand "pmfitg" est appelé par pmfd00 avec vf qui pointe sur la SD fibres
!               le nombre de composantes doit être en relation avec le type de groupe de fibres
!           !!! Quand "pmfitg" est appelé sous un TE avec PFIBRES des catalogues de poutres
!               le nombre de composantes est le maximum de tous les types (info dans PNBSP_I)
!   vf(*)   : positions des fibres
!       Types 1 et 2
!          vf(1,*) : Y fibres
!          vf(2,*) : Z fibres
!          vf(3,*) : Aire fibres
!       Types 2
!          vf(4,*) : Yp groupes de fibres
!          vf(5,*) : Zp groupes de fibres
!          vf(6,*) : GX groupes de fibres
!
! OUT
!   vs(1) : surface totale              : somme(ds)
!   vs(2) : moment statique / oz        : somme(y.ds)
!   vs(3) : moment statique / oy        : somme(z.ds)
!   vs(4) : moment quadratique / oz     : somme(y.y.ds)
!   vs(5) : moment quadratique / oy     : somme(z.z.ds)
!   vs(6) : moment produit              : somme(y.z.ds)
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterfort/utmess.h"
!
    integer(kind=8) :: typfib, nf, ncarf
    real(kind=8) :: vf(ncarf, nf), vs(6)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ii
    real(kind=8) :: yy, zz, aire
!
! --------------------------------------------------------------------------------------------------
!
    vs(:) = 0.0d0
!
    if (typfib .eq. 1) then
!       caractéristiques utiles par fibre : Y  Z  AIRE
        do ii = 1, nf
            yy = vf(1, ii)
            zz = vf(2, ii)
            aire = vf(3, ii)
!
            vs(1) = vs(1)+aire
            vs(2) = vs(2)+yy*aire
            vs(3) = vs(3)+zz*aire
            vs(4) = vs(4)+yy*yy*aire
            vs(5) = vs(5)+zz*zz*aire
            vs(6) = vs(6)+yy*zz*aire
        end do
    else if (typfib .eq. 2) then
!       caractéristiques utiles par fibre : Y  Z  AIRE  YP  ZP  GX
        do ii = 1, nf
            yy = vf(1, ii)
            zz = vf(2, ii)
            aire = vf(3, ii)
!
            vs(1) = vs(1)+aire
            vs(2) = vs(2)+yy*aire
            vs(3) = vs(3)+zz*aire
            vs(4) = vs(4)+yy*yy*aire
            vs(5) = vs(5)+zz*zz*aire
            vs(6) = vs(6)+yy*zz*aire
        end do
    else
        call utmess('F', 'ELEMENTS2_40', si=typfib)
    end if
!
end subroutine
