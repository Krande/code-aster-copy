! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
! aslint: disable=W1403
!
! ==================================================================================================
!
! Module for management of coordinates system
!
! ==================================================================================================
module coorSyst_module
! ==================================================================================================

! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: getCoorSyst
! ==================================================================================================
    private
#include "asterf_types.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getCoorSyst
!
! Get parameters of local coordinate system
!
! In  mesh             : name of mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine getCoorSyst()
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        ! character(len=*), intent(in) :: meshz
        ! type(MESH_OPER_MODI_PARA), intent(out) :: meshOperModiPara
! ----- Local
        ! integer(kind=8) :: iret, jvAngMas, coorSystType
!   ------------------------------------------------------------------------------------------------
!
!         call tecach('NNO', 'PCAMASS', 'L', iret, iad=jvAngMas)
!         anglNaut = 0.d0

!         !         C : indice de definition du repere d'orthotropie (=1 definition par 3
!         !   angles nautiques, = -1 definition par un axe et un point sur cet axe,
!         !   = 2 définition par 3 angles d'Euler ou par un champ d'orientation)
!         !   ALPHA : 1er angle nautique
!         !   BETA :  2eme angle nautique
!         !   KAPPA : 3eme angle nautique
!         !   X : nul si C=1, sinon 1ere coordonnee du point de l'axe
!         !   Y : nul si C=1, sinon 2eme coordonnee du point de l'axe
!         !   Z : nul si C=1, sinon 3eme coordonnee du point de l'axe
! !
!         if (iret .eq. 0) then
!             coorSystType = nint(zr(jvAngMas))

!             if (zr(jvAngMas) .gt. 0.d0) then
!                 angl_naut(1) = zr(jvAngMas+1)*r8dgrd()
!                 if (ndim .eq. 3) then
!                     angl_naut(2) = zr(jvAngMas+2)*r8dgrd()
!                     angl_naut(3) = zr(jvAngMas+3)*r8dgrd()
!                 end if
! !
!             else if (abs(zr(jvAngMas)+1.d0) .lt. 1.d-3) then
! !
! ! ON TRANSFORME LA DONNEE DU REPERE CYLINDRIQUE EN ANGLE NAUTIQUE
! !
!                 orig(1:ndim) = zr(jvAngMas+3+1:jvAngMas+3+ndim)
!                 if (ndim .eq. 3) then
!                     alpha = zr(jvAngMas+1)*r8dgrd()
!                     beta = zr(jvAngMas+2)*r8dgrd()
!                     dire(1) = cos(alpha)*cos(beta)
!                     dire(2) = sin(alpha)*cos(beta)
!                     dire(3) = -sin(beta)
!                     call utrcyl(coor, dire, orig, p)
!                     do i = 1, 3
!                         xg(i) = p(1, i)
!                         yg(i) = p(2, i)
!                     end do
!                     call angvxy(xg, yg, angl_naut)
!                 else
!                     xu = coor(1)-orig(1)
!                     yu = coor(2)-orig(2)
!                     xnorm = sqrt(xu**2+yu**2)
!                     xu = xu/xnorm
!                     yu = yu/xnorm
!                     p(1, 1) = xu
!                     p(2, 1) = yu
!                     p(1, 2) = -yu
!                     p(2, 2) = xu
!                     xg(1) = xu
!                     xg(2) = yu
!                     xg(3) = 0.d0
!                     call angvx(xg, alpha, beta)
!                     angl_naut(1) = alpha
!                 end if
!             end if
!         end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!===================================================================================================
!===================================================================================================
    ! public :: MESH_OPER_ORIE_SHELL, MESH_OPER_MODI_PARA
!===================================================================================================
!===================================================================================================
!
end module coorSyst_module
