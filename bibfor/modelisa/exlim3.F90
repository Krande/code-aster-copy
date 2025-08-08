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

subroutine exlim3(motfaz, base, modelz, ligrel)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim4.h"
    character(len=*) :: motfaz, base, modelz, ligrel
! person_in_charge: jacques.pellet at edf.fr
! ======================================================================
! BUT  :  SCRUTER LES MOTS CLE TOUT/GROUP_MA/MAILLE POUR CREER
!         UN LIGREL "REDUIT" A PARTIR DU LIGREL DU MODELE MODELZ
!
! LA DIFFERENCE AVEC EXLIMA.F EST QUE CETTE ROUTINE SCRUTE TOUTES LES
! OCCURENCES DE MOTFAZ ET DETERMINE L'EVELOPPE DE LA LISTE DES MAILLES
!
! IN  : MODELZ : NOM DU MODELE
!
! OUT/JXOUT   : LIGREL  : LIGREL REDUIT
!     ATTENTION :
!          - LE NOM DE LIGREL EST TOUJOURS "OUT"
!          - PARFOIS ON REND LIGREL=LIGREL(MODELE) :
!             - ALORS ON NE TIENT DONC PAS COMPTE DE 'BASE'
!             - IL NE FAUT PAS LE DETRUIRE !
!          - PARFOIS ON EN CREE UN NOUVEAU SUR LA BASE 'BASE'
!             - LE NOM DU LIGREL EST OBTENU PAR GNOMSD
!     -----------------------------------------------------------------
!
    character(len=8) :: modele
    character(len=19) :: ligrmo
!     -----------------------------------------------------------------
!
!
    modele = modelz
    ASSERT(modele .ne. ' ')
!
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
    call exlim4(motfaz, base, ligrmo, ligrel)
!
!
end subroutine
