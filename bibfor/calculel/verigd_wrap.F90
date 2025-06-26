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
subroutine verigd_wrap(nomgdz, lcmp, ncmp, iret)
!
    implicit none
#include "asterfort/verigd.h"
!
    integer(kind=8) :: ncmp, iret
    character(len=*) :: nomgdz
    character(len=8) :: lcmp(ncmp)
! ---------------------------------------------------------------------
! BUT: VERIFIER LA COHERENCE D'UNE LISTE DE CMPS D'UNE GRANDEUR.
! ---------------------------------------------------------------------
!     ARGUMENTS:
! LCMP   IN   V(K8) : LISTE DES CMPS A VERIFIER
! NCMP   IN   I     : LONGUEUR DE LA LISTE LCMP
! IRET   OUT  I     : CODE RETOUR :
!                     /0 : OK
!                     /1 : NOMGDZ N'EST PAS UNE GRANDEUR.
!                     /2 : UNE CMP (AU MOINS) EST EN DOUBLE DANS LCMP
!                     /3 : UNE CMP (AU MOINS) DE LCMP N'EST PAS UNE
!                          CMP DE NOMGDZ
!
!  SI IRET>0 : ON EMET SYSTEMATIQUEMENT UNE ALARME.
!              C'EST A L'APPELANT D'ARRETER LE CODE SI NECESSAIRE.
!----------------------------------------------------------------------
!
    call verigd(nomgdz, lcmp, ncmp, iret)
!
end subroutine
