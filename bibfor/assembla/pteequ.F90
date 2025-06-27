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

subroutine pteequ(nume_equa, base, neq, igds, nb_cmp_field, &
                  field_to_cata)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/exisdg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/wkvect.h"
!
!
    character(len=19), intent(in) :: nume_equa
    integer(kind=8), intent(in) :: neq
    integer(kind=8), intent(in) :: igds
    integer(kind=8), intent(in) :: nb_cmp_field
    integer(kind=8), pointer :: field_to_cata(:)
    character(len=1), intent(in) :: base
!
! --------------------------------------------------------------------------------------------------
!
! Create DEEQ object for field
!
! --------------------------------------------------------------------------------------------------
!
! In  base           : JEVEUX base to create NUME_EQUA object
! In  nume_equa      : name of NUME_EQUA object
! In  neq            : number of equations
! In  igds           : index of GRANDEUR used to numbering
! In  nb_cmp_field   : number of components in field
! In  field_to_cata  : pointer to converter from local components (field) to global (catalog)
!
! Object   : NUME_EQUA.DEEQ
! Dimension: vector of size (2*neq)
! Contains : for ieq = 1,neq
!
!   DEEQ((ieq-1)*2+1)
!       SI LE NOEUD SUPPORT DE L'EQUA. IDDL EST PHYS.:
!           +NUMERO DU NOEUD
!       SI LE NOEUD SUPPORT DE L'EQUA. IDDL EST UN LAGRANGE DE BLOCAGE :
!           +NUMERO DU NOEUD PHYS. BLOQUE
!       SI LE NOEUD SUPPORT DE L'EQUA. IDDL EST UN LAGRANGE DE LIAISON :
!           0
!   DEEQ((ieq-1)*2+2)
!       SI LE NOEUD SUPPORT DE L'EQUA. IDDL EST PHYS.:
!           + NUM. DANS L'ORDRE DU CATAL. DES GRAND. DE LA CMP CORRESPONDANT A L'EQUATION IDDL.
!       SI LE NOEUD SUPPORT DE L'EQUA. IDDL EST UN LAGRANGE DE BLOCAGE :
!           - NUM. DANS L'ORDRE DU CATAL. DES GRAND. DE LA CMP CORRESPONDANT AU BLOCAGE.
!       SI LE NOEUD SUPPORT DE L'EQUA. IDDL EST UN LAGRANGE DE LIAISON :
!           0
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nb_cmp_fieldmx, nec, nb_ligr, l, jprno, i_cmp_field
    integer(kind=8) :: nb_node, i_node, iddl, iadg, i_cmp_cata, i_equa
    character(len=24) :: prno, nueq, deeq
    integer(kind=8), pointer :: p_nueq(:) => null()
    integer(kind=8), pointer :: p_deeq(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    prno = nume_equa(1:19)//'.PRNO'
    nueq = nume_equa(1:19)//'.NUEQ'
    deeq = nume_equa(1:19)//'.DEEQ'
!
! - Information about GRANDEUR
!
    call jelira(jexnum('&CATA.GD.NOMCMP', igds), 'LONMAX', nb_cmp_fieldmx)
    nec = nbec(igds)
    ASSERT(nb_cmp_fieldmx .ne. 0)
    ASSERT(nec .ne. 0)
!
! - Create .DEEQ object
!
    call jedetr(deeq)
    call wkvect(deeq, base//' V I', max(1, 2*neq), vi=p_deeq)
    call jeecra(deeq, "LONUTI", 2*neq)
!
! - Access to PRNO object
!
    call jelira(prno, 'NMAXOC', nb_ligr)
! - It's field: no "tardif" element/node, only mesh
    ASSERT(nb_ligr .eq. 1)
    call jelira(jexnum(prno, 1), 'LONMAX', l)
    ASSERT(l .gt. 0)
    call jeveuo(jexnum(prno, 1), 'L', jprno)
!
! - Access to NUEQ object
!
    call jeveuo(nueq, 'L', vi=p_nueq)
!
    nb_node = l/(nec+2)
    do i_node = 1, nb_node
        iddl = zi(jprno-1+(i_node-1)*(nec+2)+1)-1
        iadg = jprno-1+(i_node-1)*(nec+2)+3
        do i_cmp_field = 1, nb_cmp_field
            i_cmp_cata = field_to_cata(i_cmp_field)
            if (exisdg(zi(iadg), i_cmp_cata)) then
                iddl = iddl+1
                i_equa = p_nueq(iddl)
                p_deeq(2*(i_equa-1)+1) = i_node
                p_deeq(2*(i_equa-1)+2) = i_cmp_cata
            end if
        end do
    end do
!
    call jedema()
end subroutine
