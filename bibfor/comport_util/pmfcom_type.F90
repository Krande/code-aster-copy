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
! person_in_charge: jean-luc.flejou at edf.fr
!
module pmfcom_type
!
    implicit none
!
    type :: pmfcom_user
        ! kpg     : numéro de point de gauss
        ! debsp   : numéro de sous-point de la première fibre du groupe
        ! option  : option de calcul
        ! instam  : instant du calcul précédent
        ! instap  : instant du calcul
        ! icdmat  : code matériau
        ! epsm    : déformation a l'instant précédent sur l'élément de structure
        !
        integer(kind=8)             :: kpg
        integer(kind=8)             :: debsp
        integer(kind=8)             :: icdmat
        real(kind=8)        :: instam
        real(kind=8)        :: instap
        real(kind=8)        :: epsm
        character(len=16)   :: option
        !
    end type pmfcom_user
!
end module pmfcom_type
