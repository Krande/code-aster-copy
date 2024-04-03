# ZK16

The following is an error message raised when compiling with `ifix`. The error message is 
not present when compiling using gfortran on linux.

`error #7938: Character length argument mismatch.   [ZK16]`

Here is how it is used in the offending fortran file

`bibfor/plate/vdxnlr.F90`

```
call behaviourOption(option, zk16(icompo), lMatr, lVect, &
lVari, lSigm, codret)
```

`icompo` is an integer variable.

here is the `behaviourOption` function (and just the first to variable declarations, whereas
the second one `compor` is related to the error) 

`bibfor/comport_util/Behaviour_module.F90`.

```
! 
subroutine behaviourOption(option, compor, &
                               lMatr, lVect, &
                               lVari, lSigm, &
                               codret_)
        !   -----------------------------------------------------------------------
        ! - Parameters
        character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
```        

Here is the definition of `COMPOR_SIZE` in the fortran file
`#define COMPOR_SIZE 25`

Here is the full file containing the definition of zk16

`bibfor/include/jeveux.h`

```fortran
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

#ifndef JEVEUX_H_
#define JEVEUX_H_
!
!---------- DEBUT COMMUNS NORMALISES  JEVEUX  --------------------------
#include "asterf_types.h"
!
volatile           zi4, zi, zr, zc, zl
volatile           zk8, zk16, zk24, zk32, zk80
!
integer(kind=4)           :: zi4
common  / i4vaje / zi4(1)
integer :: zi
common  / ivarje / zi(1)
real(kind=8)              :: zr
common  / rvarje / zr(1)
complex(kind=8)          :: zc
common  / cvarje / zc(1)
aster_logical          :: zl
common  / lvarje / zl(1)
character(len=8)         :: zk8
character(len=16)                :: zk16
character(len=24)                         :: zk24
character(len=32)                                  :: zk32
character(len=80)                                           :: zk80
common  / kvarje / zk8(1), zk16(1), zk24(1), zk32(1), zk80(1)
!---------- FIN  COMMUNS NORMALISES  JEVEUX ----------------------------
#endif
```