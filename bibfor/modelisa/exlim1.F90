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

subroutine exlim1(lismai, nbmail, modelz, basez, ligrez)
    implicit none
#include "jeveux.h"
#include "asterfort/exlim2.h"
#include "asterfort/dismoi.h"
!
    integer(kind=8) :: lismai(*), nbmail
    character(len=*) :: modelz, basez, ligrez

! But : Creer le ligrel "reduit" correspondant au ligrel d'un modele
!       mais en ne conservant que certaines mailles.
!
! in  : lismai : liste des numeros de mailles constituant le
!                ligrel a creer.
!                Remarque : les mailles de numero 0 sont ignorees.
! in  : nbmail : longueur de la liste des mailles
! in  : modelz : nom du modele referencant les mailles de lismai
!                des grels
! in  : basez  : base sur laquelle on cree le ligrel
! out : ligrez : ligrel a creer
!----------------------------------------------------------------------
!
!
    character(len=8) :: model
    character(len=19) :: ligrmo
!     ------------------------------------------------------------------
!

    model = modelz
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)

    call exlim2(lismai, nbmail, ligrmo, basez, ligrez)
end subroutine
