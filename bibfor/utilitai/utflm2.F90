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

subroutine utflm2(mailla, tabmai, nbma, dim, typmai, &
                  nbtrou, tatrou)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8), intent(in) :: mailla
    integer(kind=8), intent(in) :: nbma
    integer(kind=8), intent(in) :: tabmai(nbma)
    integer(kind=8), intent(in) :: dim
    character(len=*), intent(in) :: typmai
    integer(kind=8), intent(out) :: nbtrou
    integer(kind=8), intent(out) :: tatrou(nbma)

! person_in_charge: josselin.delmas at edf.fr
!
!     BUT:
!       FILTRER UNE LISTE DE MAILLE D'APRES LEUR DIMENSION VERSION 2
!       *           *        *                                     *
!       IDEM QUE UTFLMD MAIS AVEC UNE LISTE DE MAILLE
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   MAILLA    : NOM DU MAILLAGE
! IN   TABMAI    : LISTE DES MAILLES
! IN   NBMA      : LONGUEUR DE LA LISTE
! IN   NDIM      : DIMENSION DES MAILLES A TROUVER (-1,0,1,2,3)
! IN   TYPMAI    : SI DIM=-1, ON FILTRE SUR TYPMAI='QUAD4'/'TRIA3'/...
!                  SINON, ON NE SE SERT PAS DE TYPMAI
!      SORTIE :
!-------------
! OUT  NBTROU    : NOMBRE DE MAILLE TROUVEES
! OUT  TATROU    : LISTE DES MAILLES TROUVEES
!     REMARQUE : TATROU EST SUPPOSE DE LONGUEUR SUFFISANTE (>= NBTROU)
!
!.......................................................................
!
!
!
!
!
    integer(kind=8) :: nbtyp, i, ii, itrou, itych
    integer(kind=8), pointer :: dime_topo(:) => null()
    integer(kind=8), pointer :: liste_m_temp(:) => null()
    integer(kind=8), pointer :: liste_typmai(:) => null()
    character(len=8), pointer :: type_maille(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
!
!
! ----------------------------------------------------------------------
!
    call jemarq()
    ASSERT(nbma .gt. 0)
!
    call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbtyp)
!
!     -- SI DIM=-1, ON TRIE SUR TYPMAI :
    if (dim .eq. -1) then
        call jenonu(jexnom('&CATA.TM.NOMTM', typmai), itych)
        if (itych .eq. 0) then
            call utmess('F', 'CALCULEL2_67', sk=typmai)
        end if
!
    else
!
        AS_ALLOCATE(vk8=type_maille, size=nbtyp)
        AS_ALLOCATE(vi=dime_topo, size=nbtyp)
        AS_ALLOCATE(vi=liste_typmai, size=nbma)
!
! ------RECUPERATION DE TOUS LES TYPES DE MAILLE
!        ET DE LEUR DIMENSION TOPOLOGIQUE
!
        do i = 1, nbtyp
            call jenuno(jexnum('&CATA.TM.NOMTM', i), type_maille(i))
            call dismoi('DIM_TOPO', type_maille(i), 'TYPE_MAILLE', repi=dime_topo(i))
        end do
    end if
!
!
! ----RECUPERATION DE LA LISTE DES TYPES DE MAILLE DU MAILLAGE
!
    call jeveuo(mailla//'.TYPMAIL        ', 'L', vi=typmail)
!
    AS_ALLOCATE(vi=liste_m_temp, size=nbma)
!
    nbtrou = 0
    ii = 1
    do i = 1, nbma
!
! ------RECUPERATION DU TYPE DE LA MAILLE I
!
        if (dim .ne. -1) then
!
! --------SI LA DIMENSION TOPOLOGIQUE EST LA BONNE ON GARDE LA MAILLE
!
            liste_typmai(i) = typmail(tabmai(i))
            if (dime_topo(liste_typmai(i)) .eq. dim) then
                liste_m_temp(ii) = tabmai(i)
                nbtrou = nbtrou+1
                ii = ii+1
            end if
        else
!
! --------TRI SUR LE TYPE DE LA MAILLE :
!
            if (typmail(tabmai(i)) .eq. itych) then
                liste_m_temp(ii) = tabmai(i)
                nbtrou = nbtrou+1
                ii = ii+1
            end if
        end if
    end do
!
    if (nbtrou .eq. 0) goto 999
!
! ----SI LA LISTE N'EST PAS VIDE ON RECOPIE DANS TATROU DE TAILLE NBMA
!
    do itrou = 1, nbtrou
        tatrou(itrou) = liste_m_temp(itrou)
    end do
!
999 continue
!
    AS_DEALLOCATE(vk8=type_maille)
    AS_DEALLOCATE(vi=dime_topo)
    AS_DEALLOCATE(vi=liste_m_temp)
    AS_DEALLOCATE(vi=liste_typmai)
!
    call jedema()
!
end subroutine
