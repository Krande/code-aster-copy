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
#include "MeshTypes_type.h"
!
interface
    subroutine lrcmve(ntvale, nmatyp, nbnoma, ntproa, lgproa,&
                      ncmprf, nomcmr, ntypel, npgmax, indpg,&
                      nbcmfi, nmcmfi, nbcmpv, ncmpvm, numcmp,&
                      jnumma, nochmd, nbma, npgma, npgmm,&
                      nspmm, typech, nutyma, adsl, adsv, adsd,&
                      lrenum, nuanom, codret)
        integer(kind=8) :: nbma
        integer(kind=8) :: npgmax
        integer(kind=8) :: ntypel
        character(len=*) :: ntvale
        integer(kind=8) :: nmatyp
        integer(kind=8) :: nbnoma
        character(len=*) :: ntproa
        integer(kind=8) :: lgproa
        integer(kind=8) :: ncmprf
        character(len=*) :: nomcmr(*)
        integer(kind=8) :: indpg(ntypel, npgmax)
        integer(kind=8) :: nbcmfi
        character(len=*) :: nmcmfi
        integer(kind=8) :: nbcmpv
        character(len=*) :: ncmpvm
        character(len=*) :: numcmp
        integer(kind=8) :: jnumma
        character(len=*) :: nochmd
        integer(kind=8) :: npgma(nbma)
        integer(kind=8) :: npgmm(nbma)
        integer(kind=8) :: nspmm(nbma)
        character(len=*) :: typech
        integer(kind=8) :: nutyma
        integer(kind=8) :: adsl
        integer(kind=8) :: adsv
        integer(kind=8) :: adsd
        integer(kind=8) :: nuanom(MT_NTYMAX, MT_NNOMAX)
        aster_logical :: lrenum
        integer(kind=8) :: codret
    end subroutine lrcmve
end interface
