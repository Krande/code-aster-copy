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
!
interface
    subroutine cftels(typco, typstru, effrts, effm, effn, efft, effmt,&
                      dnsinf, dnssup,&
                      sigmsi, sigmss, sigmci, sigmcs, alpha,&
                      ht, bw, enrobi, enrobs, facier, fbeton,&
                      scmaxi, scmaxs, ssmax, uc, um,&
                      compress, dnstra, thetab, ak, uk, ierr)
        integer(kind=8) ::typco
        integer(kind=8) :: typstru
        real(kind=8) :: effrts(8)
        real(kind=8) :: effm
        real(kind=8) :: effn
        real(kind=8) :: efft
        real(kind=8) :: effmt
        real(kind=8) :: dnsinf
        real(kind=8) :: dnssup
        real(kind=8) :: sigmsi
        real(kind=8) :: sigmss
        real(kind=8) :: sigmci
        real(kind=8) :: sigmcs
        real(kind=8) :: alpha
        real(kind=8) :: ht
        real(kind=8) :: bw
        real(kind=8) :: enrobi
        real(kind=8) :: enrobs
        real(kind=8) :: facier
        real(kind=8) :: fbeton
        real(kind=8) :: scmaxi
        real(kind=8) :: scmaxs
        real(kind=8) :: ssmax
        integer(kind=8) :: uc
        integer(kind=8) :: um
        integer(kind=8) :: compress
        real(kind=8) :: dnstra
        real(kind=8) :: thetab
        real(kind=8) :: ak
        real(kind=8) :: uk
        integer(kind=8) :: ierr
    end subroutine cftels
end interface
