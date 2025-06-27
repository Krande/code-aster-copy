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

subroutine nume_equa_gene_crsd(nume_equa_genez, base, nb_equa, nb_sstr, nb_link, &
                               model_genez, gran_namez)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/jecreo.h"
#include "asterfort/jeecra.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jenonu.h"
#include "asterfort/wkvect.h"
!
!
    character(len=*), intent(in) :: nume_equa_genez
    character(len=1), intent(in) :: base
    integer(kind=8), intent(in) :: nb_equa
    integer(kind=8), intent(in) :: nb_sstr
    integer(kind=8), intent(in) :: nb_link
    character(len=*), optional, intent(in) :: model_genez
    character(len=*), optional, intent(in) :: gran_namez
!
! --------------------------------------------------------------------------------------------------
!
! NUME_EQUA_GENE
!
! Create object
!
! --------------------------------------------------------------------------------------------------
!
! In  nume_equa_gene   : name of NUME_EQUA_GENE
! In  base        : JEVEUX base to create NUME_EQUA_GENE
! In  nb_equa     : number of equations
! In  nb_sstr     : number of sub_structures
! In  nb_link     : number of links
! In  model_gene  : name of model
! In  gran_name   : name of GRANDEUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: nume_equa_gene
    integer(kind=8) :: i_equa, i_ligr_sstr, i_ligr_link, nb_ligr
    logical :: l_ligr_sstr, l_ligr_link
    character(len=24) :: model_gene, gran_name
    integer(kind=8), pointer :: prgene_nueq(:) => null()
    integer(kind=8), pointer :: prgene_deeq(:) => null()
    integer(kind=8), pointer :: prgene_delg(:) => null()
    character(len=24), pointer :: prgene_refn(:) => null()
    integer(kind=8), pointer :: prgene_desc(:) => null()
    integer(kind=8), pointer :: prgene_nequ(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nume_equa_gene = nume_equa_genez
    model_gene = ' '
    if (present(model_genez)) then
        model_gene = model_genez
    end if
    gran_name = ' '
    if (present(gran_namez)) then
        gran_name = gran_namez
    end if
!
    call detrsd('NUME_EQUA_GENE', nume_equa_gene)
!
! - Create object NUEQ
!
    call wkvect(nume_equa_gene//'.NUEQ', base//' V I', nb_equa, vi=prgene_nueq)
!
! - Set to identity
!
    do i_equa = 1, nb_equa
        prgene_nueq(i_equa) = i_equa
    end do
!
! - Create object DEEQ
!
    call wkvect(nume_equa_gene//'.DEEQ', base//' V I', 2*nb_equa, vi=prgene_deeq)
!
! - Set
!
    do i_equa = 1, nb_equa
        prgene_deeq(2*(i_equa-1)+1) = i_equa
        prgene_deeq(2*(i_equa-1)+2) = 1
    end do
!
! - Number of LIGREL
!
    l_ligr_sstr = nb_sstr .gt. 0
    l_ligr_link = nb_link .gt. 0
    nb_ligr = 0
    if (l_ligr_sstr) then
        nb_ligr = nb_ligr+1
    end if
    if (l_ligr_link) then
        nb_ligr = nb_ligr+1
    end if
    ASSERT(nb_ligr .le. 2)
!
! - Create object LILI
!
    call jecreo(nume_equa_gene//'.LILI', base//' N K8')
    call jeecra(nume_equa_gene//'.LILI', 'NOMMAX', nb_ligr)
    if (l_ligr_sstr) then
        call jecroc(jexnom(nume_equa_gene//'.LILI', '&SOUSSTR'))
        call jenonu(jexnom(nume_equa_gene//'.LILI', '&SOUSSTR'), i_ligr_sstr)
        ASSERT(i_ligr_sstr .eq. 1)
    end if
    if (l_ligr_link) then
        call jecroc(jexnom(nume_equa_gene//'.LILI', 'LIAISONS'))
        call jenonu(jexnom(nume_equa_gene//'.LILI', 'LIAISONS'), i_ligr_link)
    end if
!
! - Create object PRNO
!
    call jecrec(nume_equa_gene//'.PRNO', base//' V I', 'NU', 'DISPERSE', 'VARIABLE', nb_ligr)
    if (l_ligr_sstr) then
        call jeecra(jexnum(nume_equa_gene//'.PRNO', i_ligr_sstr), 'LONMAX', 2*nb_sstr)
    end if
    if (l_ligr_link) then
        call jeecra(jexnum(nume_equa_gene//'.PRNO', i_ligr_link), 'LONMAX', 2*nb_link)
    end if
!
! - Create object ORIG
!
    call jecrec(nume_equa_gene//'.ORIG', base//' V I', 'NU', 'DISPERSE', 'VARIABLE', 2)
    if (l_ligr_sstr) then
        call jeecra(jexnum(nume_equa_gene//'.ORIG', i_ligr_sstr), 'LONMAX', nb_sstr)
    end if
    if (l_ligr_link) then
        call jeecra(jexnum(nume_equa_gene//'.ORIG', i_ligr_link), 'LONMAX', nb_link)
    end if
!
! - Create object NEQU
!
    call wkvect(nume_equa_gene//'.NEQU', base//' V I', 1, vi=prgene_nequ)
    prgene_nequ(1) = nb_equa
!
! - Create object DESC
!
    call wkvect(nume_equa_gene//'.DESC', base//' V I', 1, vi=prgene_desc)
    prgene_desc(1) = 2
!
! - Create object REFN
!
    call wkvect(nume_equa_gene//'.REFN', base//' V K24', 4, vk24=prgene_refn)
    prgene_refn(1) = model_gene
    prgene_refn(2) = gran_name
!
! - Create object DELG
!
    call wkvect(nume_equa_gene//'.DELG', base//' V I', nb_equa, vi=prgene_delg)

end subroutine
