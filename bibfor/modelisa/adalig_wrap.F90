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

subroutine adalig_wrap(ligrz)
    implicit none
#include "asterf_types.h"
#include "asterfort/adalig.h"
!
    character(len=*), intent(in) :: ligrz
!----------------------------------------------------------------------
! But: Reorganiser la collection .LIEL de ligrz afin de regrouper
!      les elements de meme TYPE_ELEM dans un meme GREL.
!
! On veut :
!   * Limiter la taille des GRELS (pas plus de nelmx elements)
!   * Faire en sorte que l'equilibrage soit bon pour DISTRIBUTION / METHODE='GROUP_ELEM' :
!     * Pour chaque TYPE_ELEM :
!       On decoupe le paquet d'elements en un nombre de grels multiple de nbproc.
!     * Si partsd est n'est pas fourni :
!       L'equilibrage est presque parfait :
!       Les GRELS ont tous le meme nombre d'elements (a 1 pres)
!
!     * Si partsd est fourni:
!       * On ajoute une nouvelle contrainte pour les GRELS :
!         * le GREL kgrel ne contient que des elements des sous-domaines affectes au
!           processeur kproc [0, ..., nbproc-1] avec : mod(kgrel,nbproc)=kproc
!
!
! Arguments d'entree:
!     ligrz  (o) : nom du ligrel
!----------------------------------------------------------------------
    call adalig(ligrz)
end subroutine
