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
#include "asterf_types.h"
interface 
      subroutine fluendo3d(xmat,sig0,sigf,deps,&
               nstrs,var0,varf,nvari,nbelas3d,&
               teta1,teta2,dt,epstf,ierr1,&
               iso,mfr,end3d,fl3d,local,&
               ndim,nmatbe2,iteflumax,sech,&
               matrEndo, matdech)
        real(kind=8), intent(in) :: xmat(:)
        real(kind=8) :: sig0(6)
        real(kind=8) :: sigf(:)
        real(kind=8) :: deps(:)
        integer(kind=8), intent(in) :: nstrs
        real(kind=8) :: var0(:)
        real(kind=8) :: varf(:)
        integer(kind=8), intent(in) :: nvari
        integer(kind=8), intent(in) :: nbelas3d 
        real(kind=8) :: teta1
        real(kind=8) :: teta2
        real(kind=8) :: dt
        real(kind=8) :: epstf(6)
        integer(kind=8), intent(out) :: ierr1
        aster_logical :: iso
        integer(kind=8), intent(in) :: mfr
        aster_logical, intent(in) :: end3d
        aster_logical, intent(in) :: fl3d
        aster_logical, intent(in) :: local
        integer(kind=8), intent(in) :: ndim
        integer(kind=8), intent(in) :: nmatbe2
        integer(kind=8), intent(in) :: iteflumax
        real(kind=8) :: sech
        aster_logical :: matrEndo
        real(kind=8) :: matdech(6,6)
    end subroutine fluendo3d
end interface
