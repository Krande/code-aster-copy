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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine evalFaceSpeedVale(lFunc    , lTime   , time  ,&
                             nbNode   , cellDime, ipg   ,&
                             jvShape  , jvGeom  , jvLoad,&
                             speedVale,&
                             x        , y       , z_)
!
implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
!
aster_logical, intent(in) :: lFunc, lTime
integer, intent(in) :: cellDime, nbNode, ipg
integer, intent(in) :: jvGeom, jvShape, jvLoad
real(kind=8), intent(in) :: time
real(kind=8), intent(out) :: speedVale
real(kind=8), intent(out) :: x, y
real(kind=8), optional, intent(out) :: z_
!
! --------------------------------------------------------------------------------------------------
!
! Utilities for fluid
!
! Evaluation of value of speed for VITE_FACE
!
! --------------------------------------------------------------------------------------------------
!
! In  lFunc            : flag if VITE_FACE is function
! In  lTime            : flag if have time for function
! In  time             : value of current time
! In  nbNode           : total number of nodes
! In  cellDime         : dimension of cell (2 or 3)
! In  ipg              : current index of Gauss point
! In  jvShape          : JEVEUX adress for shape functions
! In  jvGeom           : JEVEUX adress for geometry (coordinates of nodes)
! In  jvLoad           : JEVEUX adress for field with parameters for load
! Out speedVale        : value of speed
! Out x, y, z          : coordinates of current Gauss point
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: z
    character(len=8) :: funcName
    integer :: i_node, ldec, iret
    integer :: nbPara
    character(len=8), parameter :: paraName(4) = (/'X   ', 'Y   ',  'Z   ', 'INST'/)
    real(kind=8) :: paraVale(4)
!
! --------------------------------------------------------------------------------------------------
!
    speedVale = 0.d0
    x = 0.d0
    y = 0.d0
    z = 0.d0

!
    if (lFunc) then
        ldec = (ipg-1)*nbNode

! ----- Coordinates of current Gauss point
        do i_node = 1, nbNode
            x = x + zr(jvGeom+(cellDime+1)*(i_node-1)-1+1) * zr(jvShape+ldec-1+i_node)
            y = y + zr(jvGeom+(cellDime+1)*(i_node-1)-1+2) * zr(jvShape+ldec-1+i_node)
            if (cellDime .eq. 2) then
                z = z + zr(jvGeom+(cellDime+1)*(i_node-1)-1+3) * zr(jvShape+ldec-1+i_node)
            endif
        end do

! ----- Evaluation of function
        nbPara = 2
        paraVale(1) = x
        paraVale(2) = y
        if (cellDime .eq. 2) then
            nbPara = 3
            paraVale(3) = z
        endif
        if (lTime) then
            nbPara = nbPara + 1
            paraVale(nbPara) = time
        endif
        funcName = zk8(jvLoad)
        call fointe('FM', funcName, nbPara, paraName, paraVale, speedVale, iret)

    else
        speedVale = zr(jvLoad)

    endif

    if (present(z_)) then
        z_ = z
    endif
!
end subroutine
