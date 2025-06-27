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
#include "asterf_types.h"
!
interface
    subroutine assma3(lmasym, lmesym, tt, igr, iel,&
                      c1, rang, jnueq, numsd, jresl,&
                      nbvel, nnoe, ldist, ldgrel,&
                      ilima, jadli, jadne, jprn1, jprn2,&
                      jnulo1, jposd1, admodl,&
                      lcmodl, mode, nec, nmxcmp, ncmp,&
                      jsmhc, jsmdi, iconx1, iconx2, jtmp2,&
                      lgtmp2, jvalm, ilinu, ellagr, nbeltb, ti1, ti2, &
                      v_crco, lcontact_par)
        aster_logical :: lmasym
        aster_logical :: lmesym
        character(len=2) :: tt
        integer(kind=8) :: igr
        integer(kind=8) :: iel
        real(kind=8) :: c1
        integer(kind=8) :: rang
        integer(kind=8) :: jnueq
        integer(kind=8), pointer :: numsd(:)
        integer(kind=8) :: jresl
        integer(kind=8) :: nbvel
        integer(kind=8) :: nnoe
        aster_logical :: ldist
        aster_logical :: ldgrel
        integer(kind=8) :: ilima
        integer(kind=8) :: jadli
        integer(kind=8) :: jadne
        integer(kind=8) :: jprn1
        integer(kind=8) :: jprn2
        integer(kind=8) :: jnulo1
        integer(kind=8) :: jposd1
        integer(kind=8) :: admodl
        integer(kind=8) :: lcmodl
        integer(kind=8) :: mode
        integer(kind=8) :: nec
        integer(kind=8) :: nmxcmp
        integer(kind=8) :: ncmp
        integer(kind=8) :: jsmhc
        integer(kind=8) :: jsmdi
        integer(kind=8) :: iconx1
        integer(kind=8) :: iconx2
        integer(kind=8) :: jtmp2
        integer(kind=8) :: lgtmp2
        integer(kind=8) :: jvalm(2)
        integer(kind=8) :: ilinu
        integer(kind=8) :: ellagr
        integer(kind=8) :: nbeltb
        integer(kind=8) :: ti1(*)
        integer(kind=8) :: ti2(*)
        integer(kind=8) :: v_crco(*)
        aster_logical :: lcontact_par
    end subroutine assma3
end interface
