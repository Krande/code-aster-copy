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
interface 
    subroutine thmCpl005(ds_thm ,&
                         lMatr, lSigm, lVari, angl_naut,&
                         j_mater,&
                         ndim   , nbvari,&
                         dimdef , dimcon,&
                         adcome , adcote, adcp11, adcp21,&
                         addeme , addete, addep1, addep2,&
                         temp   , p1    , p2    ,&
                         dtemp  , dp1   , dp2   ,&
                         deps   , epsv  , depsv ,&
                         tbiot  ,&
                         phi    , rho11 , satur ,&
                         congem , congep,&
                         vintm  , vintp , dsde  ,&
                         retcom)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        aster_logical, intent(in) :: lMatr, lSigm, lVari
        real(kind=8), intent(in) :: angl_naut(3)
        integer(kind=8), intent(in) :: j_mater, ndim, nbvari
        integer(kind=8), intent(in) :: dimdef, dimcon
        integer(kind=8), intent(in) :: adcome, adcote, adcp11, adcp21
        integer(kind=8), intent(in) :: addeme, addete, addep1, addep2
        real(kind=8), intent(in) :: temp, p1, p2
        real(kind=8), intent(in) :: dtemp, dp1, dp2
        real(kind=8), intent(in) :: epsv, depsv, deps(6), tbiot(6)
        real(kind=8), intent(out) :: phi, rho11, satur
        real(kind=8), intent(in) :: congem(dimcon)
        real(kind=8), intent(inout) :: congep(dimcon)
        real(kind=8), intent(in) :: vintm(nbvari)
        real(kind=8), intent(inout) :: vintp(nbvari)
        real(kind=8), intent(inout) :: dsde(dimcon, dimdef)
        integer(kind=8), intent(out)  :: retcom
    end subroutine thmCpl005
end interface 
