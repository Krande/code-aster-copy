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
   subroutine cafelsiter(cequi, effm, effn, ht, bw,&
                         enrobi, enrobs, scmaxi, scmaxs, ssmax, ferrcomp,&
                         precs, ferrsyme, slsyme, uc, um, condns,&
                         astend, ascomp, sstend, sscomp,&
                         sctend, sccomp,&
                         alpha, pivot, etat, ierr)  
       real(kind=8) :: cequi
       real(kind=8) :: effm
       real(kind=8) :: effn
       real(kind=8) :: ht
       real(kind=8) :: bw
       real(kind=8) :: enrobi
       real(kind=8) :: enrobs
       real(kind=8) :: scmaxi
       real(kind=8) :: scmaxs
       real(kind=8) :: ssmax
       integer(kind=8) :: ferrcomp
       integer(kind=8) :: precs
       integer(kind=8) :: ferrsyme
       real(kind=8) :: slsyme
       integer(kind=8) :: uc
       integer(kind=8) :: um
       logical :: condns
       real(kind=8) :: astend
       real(kind=8) :: ascomp
       real(kind=8) :: sstend
       real(kind=8) :: sscomp
       real(kind=8) :: sctend
       real(kind=8) :: sccomp
       real(kind=8) :: alpha
       integer(kind=8) :: pivot
       integer(kind=8) :: etat
       integer(kind=8) :: ierr
    end subroutine cafelsiter
end interface
