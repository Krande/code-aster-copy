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

subroutine irmhpc(idfimd, nomamd, nomast, nbnoeu)
! person_in_charge: nicolas.sellenet at edf.fr
!-----------------------------------------------------------------------
!     ECRITURE DU MAILLAGE -  FORMAT MED - LES NOEUDS
!        -  -     -                  -         --
!-----------------------------------------------------------------------
!     ENTREE:
!       IDFIMD  : IDENTIFIANT DU FICHIER MED
!       NOMAMD : NOM DU MAILLAGE MED
!       NOMAST : NOM UTILISATEUR DU MAILLAGE MED
!       NBNOEU : NOMBRE DE NOEUDS DU MAILLAGE
!-----------------------------------------------------------------------
!
implicit none
!
! 0.1. ==> ARGUMENTS
!
#include "asterfort/as_mmhgnw.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    med_idt :: idfimd
    integer :: nbnoeu
    character(len=*) :: nomamd, nomast
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
!
    integer, parameter :: ednoeu=3, typnoe=0
!
    integer :: codret
    integer :: jnumno
    integer :: ifm, nivinf
!
    character(len=24) :: nonulg
!
!====
! 1. PREALABLES
!====
!
    call jemarq()
!
    call infniv(ifm, nivinf)
!
    nonulg = nomast//'.NULOGL'
    call jeveuo(nonulg, 'L', jnumno)
!
    call as_mmhgnw(idfimd, nomamd, ednoeu, typnoe, zi(jnumno), nbnoeu, codret)
    print*, "COD: ", codret
!
    call jedema()
!
end subroutine
