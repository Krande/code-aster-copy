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

subroutine mtdsc2(matas, objet, eoul, adress)
    implicit none
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
    character(len=*) :: matas, objet, eoul
    integer(kind=8) :: adress
!
! -----------------------------------------------------------
!   BUT : RECUPERER L'ADRESSE D'UN OBJET D'UNE MATR_ASSE
!
!   MATAS K19  IN/JXIN : NOM DE LA MATR_ASSE
!   OBJET K4   IN      : NOM D'UN SUFFIXE : ABLO,ADIA,...
!   EOUL  K1   IN      : 'E' : EN ECRITURE
!                        'L' : EN LECTURE
!   ADRESS I   OUT     : ADRESSE DANS ZI, ZR, ... DE L'OBJET
!
!
!  ATTENTION : CETTE ROUTINE NE FAIT PAS JEMARQ/JEDEMA
!              POUR NE PAS INVALIDER L'ADRESSE "OUT"
! -----------------------------------------------------------
    character(len=19) :: mat
    character(len=14) :: nu
    character(len=4) :: obj
    integer(kind=8) :: i1, i2
    character(len=24), pointer :: refa(:) => null()
!
    mat = matas
    obj = objet
    ASSERT(obj(3:4) .eq. 'BL' .or. obj(3:4) .eq. 'DI' .or. obj(3:4) .eq. 'HC')
!
    call jeveuo(mat//'.REFA', eoul, vk24=refa)
    nu = refa(2)
!
    call jeexin(nu//'.SMOS.SMDI', i1)
    call jeexin(nu//'.SLCS.SCDI', i2)
!       SI LDLT, LES 2 OBJETS PEUVENT EXISTER
    ASSERT(i1 .gt. 0 .or. i2 .gt. 0)
!
    if (obj(2:2) .eq. 'X') then
!           -- ON PRIVILEGIE LE STOCKAGE MORSE :
        if (i1 .gt. 0) then
            obj = 'SM'//obj(3:4)
            ASSERT(obj .ne. 'SMBL')
            call jeveuo(nu//'.SMOS.'//obj, eoul, adress)
        else
            ASSERT(.false.)
!            -- LE STOCK. LIGNE_CIEL N'EXISTE QU'AVEC UN STOCK. MORSE
            obj = 'SC'//obj(3:4)
            call jeveuo(nu//'.SLCS.'//obj, eoul, adress)
        end if
!
    else if (obj(2:2) .eq. 'M') then
        ASSERT(i1 .gt. 0)
        ASSERT(obj .ne. 'SMBL')
        call jeveuo(nu//'.SMOS.'//obj, eoul, adress)
!
    else if (obj(2:2) .eq. 'C') then
        ASSERT(i2 .gt. 0)
        call jeveuo(nu//'.SLCS.'//obj, eoul, adress)
!
    else
        ASSERT(.false.)
    end if
!
end subroutine
