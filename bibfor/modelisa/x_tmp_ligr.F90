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

subroutine x_tmp_ligr(mesh, ligrel, list_cells, n_list_cells)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/cormgi.h"
#include "asterfort/dismoi.h"
#include "asterfort/initel.h"
#include "asterfort/jeecra.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utflm2.h"
#include "asterfort/wkvect.h"
    character(len=8), intent(in) :: mesh
    character(len=19), intent(inout) :: ligrel
    character(len=19), optional, intent(in) :: list_cells
    integer(kind=8), optional, intent(in) :: n_list_cells
!
! person_in_charge: sam.cuvilliez at edf.fr
! ----------------------------------------------------------------------
!
! DEFI_FISS_XFEM / PROPA_FISS :
! -----------------------------
! - creer un ligrel temporaire (dans la base volatile) defini sur une
!   liste de mailles principales du maillage "mesh"
! - si les arguments optionnels "list_cells" et "n_list_cells" sont
!   absents, on prend toutes les mailles principales de "mesh"
!
! Version simplifiee de ce qui est fait dans op0018 (AFFE_MODELE)
!
! ----------------------------------------------------------------------
!
!    in     : mesh   -> nom du maillage
!
!    in/out : ligrel -> ligrel cree uniquement a partir des mailles
!                       principales de "mesh"
!
!    optionnel, in : list_cells   -> nom d'un vecteur jeveux contenant
!                                    liste de numero de mailles
!
!    optionnel, in : n_list_cells -> longueur du vecteur "list_cells"
!
! ----------------------------------------------------------------------
!
! Attention : "list_cells" doit etre une sous-liste de la liste
! de toutes les mailles principales de "mesh"
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nu_ma, nu_typ_ma, dim_ma, nu_typ_el1, nu_typ_el2
    integer(kind=8) :: i, nbma, ndim, nmaprin, idx_modeli, nb_grel, lont_liel
    integer(kind=8) :: ncells, nume_grel, long_grel, idx_in_grel
    integer(kind=8), pointer :: p_ligrel_nbno(:) => null()
    integer(kind=8), pointer :: lmatout(:) => null()
    integer(kind=8), pointer :: lmatmp(:) => null()
    integer(kind=8), pointer :: lmaprin(:) => null()
    integer(kind=8), pointer :: p_cata_typ_el(:) => null()
    integer(kind=8), pointer :: p_mesh_type_geom(:) => null()
    integer(kind=8), pointer :: p_dim_topo(:) => null()
    integer(kind=8), pointer :: p_liel(:) => null()
    integer(kind=8), pointer :: p_list_cells(:) => null()
    character(len=8), pointer :: p_ligrel_lgrf(:) => null()
    character(len=16) :: phenom, modeli(3)
    character(len=24) :: liel, lgrf, nbno
    aster_logical :: all_cells
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! - Check presence and coherence of optional arguments
!
    all_cells = .true.
!
    if (present(list_cells)) then
        ASSERT(present(n_list_cells))
        all_cells = .false.
    end if
    if (present(n_list_cells)) then
        ASSERT(present(list_cells))
    end if
!
! - Give a name for temporary LIGREL used to compute GRAD_NEUT_R option
!
    liel = ligrel//'.LIEL'
    lgrf = ligrel//'.LGRF'
    nbno = ligrel//'.NBNO'
!
! - Common definition for ligrel SD
!
    call wkvect(lgrf, 'V V K8', 2, vk8=p_ligrel_lgrf)
    call wkvect(nbno, 'V V I', 1, vi=p_ligrel_nbno)
    call jeecra(lgrf, 'DOCU', cval='MECA')
    p_ligrel_lgrf(1) = mesh
!   no sd_model => no sd_partition to get for this ligrel :
!     -> options will be computed sequentially
    p_ligrel_lgrf(2) = ''
    p_ligrel_nbno(1) = 0
!
! - Get mesh dimension "ndim" (2 or 3)
!
    call dismoi('DIM_GEOM', mesh, 'MAILLAGE', repi=ndim)
    ASSERT(ndim .eq. 2 .or. ndim .eq. 3)
!
! - Get list of cells to treat
!
    if (all_cells) then
!
! ----- Get list of all cells of dimension "ndim" that belong to mesh
!
        call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbma)
        ASSERT(nbma .gt. 0)
        AS_ALLOCATE(vi=lmatout, size=nbma)
        AS_ALLOCATE(vi=lmatmp, size=nbma)
        do i = 1, nbma
            lmatout(i) = i
            lmatmp(i) = 0
        end do
!
        call utflm2(mesh, lmatout, nbma, ndim, ' ', nmaprin, lmatmp)
        ASSERT(nmaprin .gt. 0)
        AS_ALLOCATE(vi=lmaprin, size=nmaprin)
        do i = 1, nmaprin
            lmaprin(i) = lmatmp(i)
        end do
        AS_DEALLOCATE(vi=lmatout)
        AS_DEALLOCATE(vi=lmatmp)
!
    else
!
! ----- Get list of all cells "p_list_cells" (optional argument)
!
        call jeveuo(list_cells, 'L', vi=p_list_cells)
!
    end if
!
! - Get number of cells to treat
!
    if (all_cells) then
        ncells = nmaprin
    else
        ncells = n_list_cells
    end if
    ASSERT(ncells .gt. 0)
!
! - Set modelisation type acording to dimension "ndim"
!
    phenom = 'PRESENTATION    '
    modeli(1) = '                '
    modeli(2) = '2D_GEOM         '
    modeli(3) = '3D_GEOM         '
!
! - Get cell type ids in mesh
!
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=p_mesh_type_geom)
!
! - Get element type ids available in CATA for the selected modelisation
!
    call jenonu(jexnom('&CATA.'//phenom(1:13)//'.MODL', modeli(ndim)), idx_modeli)
    call jeveuo(jexnum('&CATA.'//phenom, idx_modeli), 'L', vi=p_cata_typ_el)
!
! - Count number of GREL - Elements
!
    nb_grel = 0
    nu_typ_el2 = 0
    do i = 1, ncells
!
        if (all_cells) then
            nu_ma = lmaprin(i)
        else
            nu_ma = p_list_cells(i)
        end if
!
        nu_typ_ma = p_mesh_type_geom(nu_ma)
!       current cell dimension != mesh dimension => fatal error
        call jeveuo(jexnum('&CATA.TM.TMDIM', nu_typ_ma), 'L', vi=p_dim_topo)
        dim_ma = p_dim_topo(1)
        ASSERT(dim_ma .eq. ndim)
!       affectation of an element on current cell not possible => fatal error
        nu_typ_el1 = p_cata_typ_el(nu_typ_ma)
        ASSERT(nu_typ_el1 .gt. 0)
        if (nu_typ_el1 .ne. nu_typ_el2) then
            nu_typ_el2 = nu_typ_el1
            nb_grel = nb_grel+1
        end if
    end do
!
! - Create LIEL
!
    lont_liel = nb_grel+ncells
    call jecrec(liel, 'V V I', 'NU', 'CONTIG', 'VARIABLE', nb_grel)
    call jeecra(liel, 'LONT', lont_liel)
    call jeveuo(liel, 'E', vi=p_liel)
!
! - Store GREL in LIEL - Elements
!
    nu_typ_el2 = 0
    nume_grel = 0
    long_grel = 0
    idx_in_grel = 0
    do i = 1, ncells
!
        if (all_cells) then
            nu_ma = lmaprin(i)
        else
            nu_ma = p_list_cells(i)
        end if
!
        nu_typ_ma = p_mesh_type_geom(nu_ma)
        nu_typ_el1 = p_cata_typ_el(nu_typ_ma)
!
! ----- Create new GREL
!
        if (nu_typ_el1 .ne. nu_typ_el2 .and. nu_typ_el2 .ne. 0) then
            nume_grel = nume_grel+1
            long_grel = long_grel+1
            idx_in_grel = idx_in_grel+1
            p_liel(idx_in_grel) = nu_typ_el2
            call jecroc(jexnum(liel, nume_grel))
            call jeecra(jexnum(liel, nume_grel), 'LONMAX', long_grel)
            long_grel = 0
        end if
!
! ------ Add element in GREL
!
        long_grel = long_grel+1
        idx_in_grel = idx_in_grel+1
        p_liel(idx_in_grel) = nu_ma
        nu_typ_el2 = nu_typ_el1
!
! ----- Last element
!
        if (i .eq. ncells) then
            nume_grel = nume_grel+1
            long_grel = long_grel+1
            idx_in_grel = idx_in_grel+1
            p_liel(idx_in_grel) = nu_typ_el2
            call jecroc(jexnum(liel, nume_grel))
            call jeecra(jexnum(liel, nume_grel), 'LONMAX', long_grel)
        end if
    end do
!
! - Automatic GREL size adaptation
!
    call adalig(ligrel)
!
! - Set element/(IGREL,IM) object
!
    call cormgi('V', ligrel)
!
! - Init elements for this LIGREL
!
    call initel(ligrel)
!
! - De-allocation of tmp arrays
!
    if (all_cells) then
        AS_DEALLOCATE(vi=lmaprin)
    end if
!
    call jedema()
!
end subroutine
