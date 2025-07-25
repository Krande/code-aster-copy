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
subroutine cjsncv(roucjs, nitimp, iter, ndt, nvi, &
                  umess, erimp, epsd, deps, sigd, &
                  vind)
!
!  DUMP EN CAS NON CONVERGENCE ITE INTERNES CJS
!
    implicit none
#include "asterfort/utmess.h"
    character(len=*) :: roucjs
    integer(kind=8) :: nitimp, iter, ndt, nvi, umess
    real(kind=8) :: erimp(nitimp, 3)
    real(kind=8) :: epsd(ndt), deps(ndt), sigd(ndt), vind(nvi)
!
    integer(kind=8) :: i
    write (umess, 2001)
2001 format(&
   &       t3, ' ITER', t10, ' ERR1=DDY',&
   &       t30, 'ERR2=DY', t50, 'ERR=DDY/DY')
    do i = 1, min(nitimp, iter)
        write (umess, 1000) i, erimp(i, 1), erimp(i, 2), erimp(i, 3)
    end do
1000 format(&
   &       t3, i4, t10, e12.5,&
   &       t30, e12.5, t50, e12.5)
    call utmess('F', 'ALGORITH2_18')
    write (6, 1002) (i, epsd(i), i=1, ndt)
    write (umess, *) ' DEPS '
    write (6, 1002) (i, deps(i), i=1, ndt)
    write (umess, *) ' SIGD '
    write (6, 1002) (i, sigd(i), i=1, ndt)
    write (umess, *) ' VIND '
    write (6, 1002) (i, vind(i), i=1, nvi)
1002 format(2x, i5, 2x, e12.5)
end subroutine
