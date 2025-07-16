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
      subroutine hydramat3d(hyd0,hydr,hyds,young00,young,&
                nu00,nu,rt00,rt,ref00,&
                ref,rc00,rc,delta00,delta,&
                beta00,beta,gft00,gft,ept00,&
                ept,pglim,epsm00,epsm,xnsat00,&
                xnsat,biotw00,biotw,&
                krgi00,krgi,iso,lambda,mu,&
                rt33,rtg33,ref33,raideur66,souplesse66,&
                xmt,dtiso,err1)

        real(kind=8), intent(in) :: hyd0
        real(kind=8), intent(in) :: hydr
        real(kind=8), intent(in) :: hyds
        real(kind=8), intent(in) :: young00
        real(kind=8), intent(out) :: young
        real(kind=8), intent(in) :: nu00
        real(kind=8), intent(out) :: nu
        real(kind=8), intent(in) :: rt00
        real(kind=8), intent(out) :: rt
        real(kind=8), intent(in) :: ref00
        real(kind=8), intent(out) :: ref
        real(kind=8), intent(in) :: rc00
        real(kind=8), intent(out) :: rc
        real(kind=8), intent(in) :: delta00
        real(kind=8), intent(out) :: delta
        real(kind=8), intent(in) :: beta00
        real(kind=8), intent(out) :: beta
        real(kind=8), intent(in) :: gft00
        real(kind=8), intent(out) :: gft
        real(kind=8), intent(in) :: ept00
        real(kind=8), intent(out) :: ept
        real(kind=8), intent(out) :: pglim
        real(kind=8), intent(in) :: epsm00
        real(kind=8), intent(out) :: epsm
        real(kind=8), intent(in) :: xnsat00
        real(kind=8), intent(out) :: xnsat
        real(kind=8), intent(in) :: biotw00
        real(kind=8), intent(out) :: biotw
        real(kind=8), intent(in) :: krgi00
        real(kind=8), intent(out) :: krgi
        aster_logical, intent(in) :: iso
        real(kind=8) :: lambda
        real(kind=8) :: mu
        real(kind=8) :: rt33(3,3)
        real(kind=8) :: rtg33(3,3)
        real(kind=8) :: ref33(3,3)
        real(kind=8) :: raideur66(6,6)
        real(kind=8) :: souplesse66(6,6)
        real(kind=8), intent(out) :: xmt
        aster_logical, intent(out) :: dtiso
        integer(kind=8) :: err1
    end subroutine hydramat3d
end interface
