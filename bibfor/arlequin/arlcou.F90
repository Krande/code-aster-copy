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
subroutine arlcou(mesh, iocc, model, typmai, &
                  nom1, nom2, partKine, &
                  dime, listRela, load)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/jemarq.h"
#include "asterfort/jelira.h"
#include "asterfort/arlcpl.h"
#include "asterfort/jedema.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), intent(in) :: iocc
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: typmai
    character(len=10), intent(in) :: nom1, nom2
    character(len=8), intent(in) :: partKine(2)
    integer(kind=8), intent(in) :: dime
    character(len=19), intent(in) :: listRela
    character(len=8), intent(in) :: load

! ----------------------------------------------------------------------

! ROUTINE ARLEQUIN

! CALCUL DES MATRICES DE COUPLAGE ELEMENTAIRE ARLEQUIN
! ROUTINE D'AIGUILLAGE SUIVANT GROUPE MEDIATEUR

! ----------------------------------------------------------------------

! IN  MAIL   : NOM DU MAILLAGE
! IN  NOMO   : NOM DU MODELE
! IN  TYPMAI : SD CONTENANT NOM DES TYPES ELEMENTS (&&CATA.NOMTM)
! IN  NOM1   : NOM DE LA SD DE STOCKAGE MAILLES GROUP_MA_1
! IN  NOM2   : NOM DE LA SD DE STOCKAGE MAILLES GROUP_MA_2
! IN  CINE   : CINEMATIQUES DES GROUPES DE MAILLE
! IN  DIME   : DIMENSION DE L'ESPACE GLOBAL (2 OU 3)

! ----------------------------------------------------------------------

    character(len=10) :: nomgr1, nomgr2
    character(len=8) :: cine1, cine2
    integer(kind=8) :: nbma1, nbma2
    character(len=19) :: ngrm1, ngrm2
    character(len=8) ::  k8bid

! ----------------------------------------------------------------------
    call jemarq()

! - INITIALISATIONS

    nomgr1 = nom1
    nomgr2 = nom2
    cine1 = partKine(1)
    cine2 = partKine(2)

! --- CALCUL DES MATRICES DE COUPLAGE ELEMENTAIRE ARLEQUIN

    ngrm1 = nom1(1:10)//'.GROUPEMA'
    ngrm2 = nom2(1:10)//'.GROUPEMA'
    call jelira(ngrm1, 'LONMAX', nbma1, k8bid)
    call jelira(ngrm2, 'LONMAX', nbma2, k8bid)

    call arlcpl(iocc, nbma1, nbma2, &
                mesh, model, typmai, &
                nomgr1, nomgr2, dime, listRela, load)

    call jedema()

end subroutine
