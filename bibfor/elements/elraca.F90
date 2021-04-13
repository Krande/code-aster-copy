! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
subroutine elraca(elrefz   ,&
                  nbfpg_   , fapg_    , nbpg_,&
                  ndim_    , nno_     , nnos_,&
                  nodeCoor_, cellVolu_)
!
implicit none
!
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "asterfort/elrfno.h"
!
character(len=*), intent(in)            :: elrefz
integer, optional, intent(out)          :: ndim_, nno_, nnos_, nbfpg_, nbpg_(MT_NBFAMX)
real(kind=8), optional, intent(out)     :: nodeCoor_(3*MT_NNOMAX), cellVolu_
character(len=8), optional, intent(out) :: fapg_(MT_NBFAMX)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Get list of integration schemes of geometric support 
!
! --------------------------------------------------------------------------------------------------
!
! In  elrefa           : name of geometric support
! Out ndim             : topological dimension (0/1/2/3)
! Out nno              : number of nodes
! Out nnos             : number of middle nodes
! Out nbfpg            : number of families of integration schemes
! Out fapg             : name of families for all integration schemes
! Out nbpg             : number of points for all integration schemes
! Out nodeCoor         : coordinates of node of geometric support in parametric space
! Out cellVolu         : volume of geometric support in parametric space
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: elrefa
    integer :: i, deca
    integer :: ndim, nno, nnos, nbfpg, nbpg(MT_NBFAMX)
    character(len=8) :: fapg(MT_NBFAMX)
    real(kind=8) :: coorno(3, MT_NNOMAX), nodeCoor(3*MT_NNOMAX), cellVolu
!
! --------------------------------------------------------------------------------------------------
!
    elrefa = elrefz
    ndim     = 0
    nno      = 0
    nnos     = 0
    nbfpg    = 0
    nbpg     = 0
    nodeCoor = 0.d0
    fapg     = ' '
    cellVolu = 0.d0

! - Get main parameters of geometric support in parametric space
    call elrfno(elrefa, nno, nnos, ndim, coorno, cellVolu)

    select case (elrefa)
        case('HE8')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 8, 27, 16, 64]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG8'
            fapg(5) = 'FPG27'
            fapg(6) = 'FPG8NOS'
            fapg(7) = 'FPG64'

        case('HE9')
            nbfpg = 6
            nbpg(1:nbfpg) = [nno, nnos, 1, 5, 7, 8]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'LOB5'
            fapg(5) = 'LOB7'
            fapg(6) = 'FPG8'

        case('H20')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 8, 27, 16, 64]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG8'
            fapg(5) = 'FPG27'
            fapg(6) = 'FPG8NOS'
            fapg(7) = 'FPG64'

        case('H27')
            nbfpg = 6
            nbpg(1:nbfpg) = [nno, nnos, 1, 8, 27, 64]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG8'
            fapg(5) = 'FPG27'
            fapg(6) = 'FPG64'

        case('TE4')
            nbfpg = 9
            nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 11, 15, 23, 8]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG4'
            fapg(5) = 'FPG5'
            fapg(6) = 'FPG11'
            fapg(7) = 'FPG15'
            fapg(8) = 'FPG23'
            fapg(9) = 'FPG4NOS'

        case('T10')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 15, 8]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG4'
            fapg(5) = 'FPG5'
            fapg(6) = 'FPG15'
            fapg(7) = 'FPG4NOS'

        case('T15')
            nbfpg = 9
            nbpg(1:nbfpg) = [nno, nnos, 1, 4, 5, 11, 15, 23, 8]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG4'
            fapg(5) = 'FPG5'
            fapg(6) = 'FPG11'
            fapg(7) = 'FPG15'
            fapg(8) = 'FPG23'
            fapg(9) = 'FPG4NOS'

        case('PE6')
            nbfpg = 8
            nbpg(1:nbfpg) = [nno, nnos, 1, 6, 6, 8, 21, 12]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG6'
            fapg(5) = 'FPG6B'
            fapg(6) = 'FPG8'
            fapg(7) = 'FPG21'
            fapg(8) = 'FPG6NOS'

        case('PE7')
            nbfpg = 5
            nbpg(1:nbfpg) = [nno, nnos, 1, 5, 7]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'LOB5'
            fapg(5) = 'LOB7'

        case('P15')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 6, 8, 21, 12]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG6'
            fapg(5) = 'FPG8'
            fapg(6) = 'FPG21'
            fapg(7) = 'FPG6NOS'

        case('P18')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 6, 8, 21, 12]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG6'
            fapg(5) = 'FPG8'
            fapg(6) = 'FPG21'
            fapg(7) = 'FPG6NOS'

        case('P21')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 6, 8, 21, 12]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG6'
            fapg(5) = 'FPG8'
            fapg(6) = 'FPG21'
            fapg(7) = 'FPG6NOS'

        case('PY5')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 5, 6, 27, 10]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG5'
            fapg(5) = 'FPG6'
            fapg(6) = 'FPG27'
            fapg(7) = 'FPG5NOS'

        case('P13')
            nbfpg = 8
            nbpg(1:nbfpg) = [nno, nnos, 1, 5, 6, 10, 27, 10]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG5'
            fapg(5) = 'FPG6'
            fapg(6) = 'FPG10'
            fapg(7) = 'FPG27'
            fapg(8) = 'FPG5NOS'

        case('P19')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 5, 6, 27, 10]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG5'
            fapg(5) = 'FPG6'
            fapg(6) = 'FPG27'
            fapg(7) = 'FPG5NOS'

        case('TR3')
            nbfpg = 13
            nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 3, 6, 13, 16, 6]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG3'
            fapg(5) = 'FPG4'
            fapg(6) = 'FPG6'
            fapg(7) = 'FPG7'
            fapg(8) = 'FPG12'
            fapg(9) = 'COT3'
            fapg(10) = 'FPG3NOS'
            fapg(11) = 'FPG13'
            fapg(12) = 'FPG16'
            fapg(13) = 'SIMP'

        case('TR6')
            nbfpg = 11
            nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 6, 13, 16]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG3'
            fapg(5) = 'FPG4'
            fapg(6) = 'FPG6'
            fapg(7) = 'FPG7'
            fapg(8) = 'FPG12'
            fapg(9) = 'FPG3NOS'
            fapg(10) = 'FPG13'
            fapg(11) = 'FPG16'

        case('TR7')
            nbfpg = 10
            nbpg(1:nbfpg) = [nno, nnos, 1, 3, 4, 6, 7, 12, 13, 16]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG3'
            fapg(5) = 'FPG4'
            fapg(6) = 'FPG6'
            fapg(7) = 'FPG7'
            fapg(8) = 'FPG12'
            fapg(9) = 'FPG13'
            fapg(10) = 'FPG16'

        case('QU4')
            nbfpg = 8
            nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 16, 2, 8]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG4'
            fapg(5) = 'FPG9'
            fapg(6) = 'FPG16'
            fapg(7) = 'FIS2'
            fapg(8) = 'FPG4NOS'

        case('QU8')
            nbfpg = 7
            nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 9, 8]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG4'
            fapg(5) = 'FPG9'
            fapg(6) = 'FPG9COQ'
            fapg(7) = 'FPG4NOS'

        case('QU9')
            nbfpg = 6
            nbpg(1:nbfpg) = [nno, nnos, 1, 4, 9, 9]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG4'
            fapg(5) = 'FPG9'
            fapg(6) = 'FPG9COQ'

        case('SE2')
            nbfpg = 13
            nbpg(1:nbfpg) = [nno, nnos, 1, 2, 3, 4, 3, 5, 4, 5, 10, nnos+2, nnos+3]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG2'
            fapg(5) = 'FPG3'
            fapg(6) = 'FPG4'
            fapg(7) = 'SIMP'
            fapg(8) = 'SIMP1'
            fapg(9) = 'COTES'
            fapg(10) = 'COTES1'
            fapg(11) = 'COTES2'
            fapg(12) = 'FPG2NOS'
            fapg(13) = 'FPG3NOS'

        case('SE3')
            nbfpg = 10
            nbpg(1:nbfpg) = [nno, nnos, 1, 2, 3, 4, 3, 4, nnos+2, nnos+3]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG2'
            fapg(5) = 'FPG3'
            fapg(6) = 'FPG4'
            fapg(7) = 'SIMP'
            fapg(8) = 'COTES'
            fapg(9) = 'FPG2NOS'
            fapg(10) = 'FPG3NOS'

        case('SE4')
            nbfpg = 8
            nbpg(1:nbfpg) = [nno, nnos, 1, 2, 3, 4, 3, 4]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'
            fapg(4) = 'FPG2'
            fapg(5) = 'FPG3'
            fapg(6) = 'FPG4'
            fapg(7) = 'SIMP'
            fapg(8) = 'COTES'

        case('PO1')
            nbfpg = 3
            nbpg(1:nbfpg) = [1, 1, 1]
            fapg(1) = 'NOEU'
            fapg(2) = 'NOEU_S'
            fapg(3) = 'FPG1'

        case default
            ASSERT(ASTER_FALSE)

    end select
!
    if (present(nodeCoor_)) then
        do i = 1, nno
            deca = ndim * (i -1)
            if (ndim .ge. 1) nodeCoor_(deca+1) = coorno(1,i)
            if (ndim .ge. 2) nodeCoor_(deca+2) = coorno(2,i)
            if (ndim .eq. 3) nodeCoor_(deca+3) = coorno(3,i)
        end do
    endif
    if (present(ndim_)) then
        ndim_ = ndim
    endif
    if (present(nno_)) then
        nno_ = nno
    endif
    if (present(nnos_)) then
        nnos_ = nnos
    endif
    if (present(nbfpg_)) then
        nbfpg_ = nbfpg
    endif
    if (present(nbpg_)) then
        nbpg_ = nbpg
    endif
    if (present(cellVolu_)) then
        cellVolu_ = cellVolu
    endif
    if (present(fapg_)) then
        fapg_ = fapg
    endif
!
end subroutine
