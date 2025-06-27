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
#include "asterf_types.h"
!
interface
    subroutine vppara(modes, typcon, knega, lraide, lmasse,&
                      lamor, mxresf, neq, nfreq, omecor,&
                      dlagr, dbloq, vectr, vectc, nbpari,&
                      nbparr, nbpark, nopara, mod45, resui,&
                      resur, resuk, ktyp, lcomod, icom1,&
                      icom2, typres, nfreqg)
        character(len=8) :: modes
        character(len=16) :: typcon
        character(len=8) :: knega
        integer(kind=8) :: lraide
        integer(kind=8) :: lmasse
        integer(kind=8) :: lamor
        integer(kind=8) :: mxresf
        integer(kind=8) :: neq
        integer(kind=8) :: nfreq
        real(kind=8) :: omecor
        integer(kind=8) :: dlagr(*)
        integer(kind=8) :: dbloq(*)
        real(kind=8) :: vectr(*)
        complex(kind=8) :: vectc(*)
        integer(kind=8) :: nbpari
        integer(kind=8) :: nbparr
        integer(kind=8) :: nbpark
        character(len=*) :: nopara(*)
        character(len=4) :: mod45
        integer(kind=8) :: resui(*)
        real(kind=8) :: resur(*)
        character(len=*) :: resuk(*)
        character(len=1) :: ktyp
        aster_logical :: lcomod
        integer(kind=8) :: icom1
        integer(kind=8) :: icom2
        character(len=16) :: typres
        integer(kind=8) :: nfreqg
    end subroutine vppara
end interface
