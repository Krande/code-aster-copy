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

subroutine xchdec(modelx, decou, chdec)
! person_in_charge: patrick.massin at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/cescel.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/typele.h"
!
    character(len=8) :: modelx
    character(len=8) :: decou
    character(len=19) :: chdec

! ROUTINE XFEM (METHODE XFEM - CREATION CHAM_ELEM)
!
! CREATION D'UN CHAM_ELEM_S VIERGE POUR ETENDRE LE CHAM_ELEM
! A PARTIR DE LA STRUCTURE D UN CHAMP EXISTANT
!
! ----------------------------------------------------------------------
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  NBMA   : NOMBRE DE MAILLES
! IN/OUT  TRAV   : TABLE TO BE FILLED ON FINDING EPSI
! OUT CHAM_ELEM :
!
    integer(kind=8) :: nncp, ibid
    character(len=19) :: ligrel, cham_elem_s
    integer(kind=8) :: nute, igr, nbgrel
    integer(kind=8) :: nel, jliel, jcesd, jcesl
    integer(kind=8) :: ima, iad, iel
    character(len=8) :: noma
    character(len=16) :: nomte
    character(len=8), pointer :: cesv(:) => null()
!
!     ------------------------------------------------------------------
!
    call jemarq()
    ligrel = modelx//'.MODELE'
    call dismoi('NOM_MAILLA', modelx, 'MODELE', repk=noma)
!
    call jelira(ligrel//'.LIEL', 'NMAXOC', nbgrel)
!
! allocation du CHAM_ELEM_S cham_elem_s
    cham_elem_s = '&&XCHAMELE.CHAMDIS'
    call cescre('V', cham_elem_s, 'ELEM', noma, 'NEUT_K8', &
                0, ' ', [-1], [-1], [-1])
    call jeveuo(cham_elem_s//'.CESD', 'L', jcesd)
    call jeveuo(cham_elem_s//'.CESL', 'E', jcesl)
    call jeveuo(cham_elem_s//'.CESV', 'E', vk8=cesv)

!   1. RECUPERATION DE TYPE DE FACETTES A GENERER
!     ------------------------------------------------------------------
    do igr = 1, nbgrel
        nel = nbelem(ligrel, igr)
        call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', jliel)
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
!
!       calcul du nombre de points de Gauss pour chaque element
!       du groupe d'elements
        do iel = 1, nel
            ima = zi(jliel-1+iel)
            if (ima .lt. 0) cycle
!
!           stockage du type de discontinuite
            call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
            iad = abs(iad)
            zl(jcesl-1+iad) = .true.
            cesv(iad) = decou
        end do
    end do
!
!---CONVERSION CHAM_ELEM_S -> CHAM_ELEM
    call cescel(cham_elem_s, ligrel, 'TOPOFA', 'PDECOU', 'NON', &
                nncp, 'V', chdec, 'F', ibid)
!
    call detrsd('CHAM_ELEM_S', cham_elem_s)
!
    call jedema()
end subroutine
