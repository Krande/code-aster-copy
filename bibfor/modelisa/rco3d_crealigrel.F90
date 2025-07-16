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
subroutine rco3d_crealigrel(ligrel, noma, mod, list_pairs, nb_pairs, nt_nodes, &
                            list_total_no_co, nbnocot, map_noco_pair, map_noco_nbelem, &
                            map_noco_nbnoco)
    !
    implicit none
    !
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/initel.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jedema.h"
#include "asterfort/jedupo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/rco3d_elem.h"
!
    character(len=19), intent(in) :: ligrel
    character(len=8), intent(in) :: noma, mod
    integer(kind=8), intent(in):: nb_pairs, nt_nodes, nbnocot
    integer(kind=8), intent(out) :: map_noco_pair(:, :, :)
    integer(kind=8), intent(out) :: map_noco_nbelem(:, :)
    integer(kind=8), intent(out) :: map_noco_nbnoco(:, :, :)
    integer(kind=8), pointer :: list_total_no_co(:)
    integer(kind=8), pointer :: list_pairs(:)

! -------------------------------------------------------------------------------
!  SUBROUTINE: rco3d_crealigrel
!
!  DESCRIPTION:
!  This subroutine creates a LIGREL (relation link) for a 3D-shell connection.
!  It computes and maps node pairs between the shell and 3D interface elements.
!
!  INPUT PARAMETERS:
!  -----------------
!  ligrel           - IN(OUT)- K19  - Name of the LIGREL to create
!  noma             - IN     - K8   - Name of the mesh
!  mod              - IN     - K8   - Name of the model
!  nt_nodes         - IN     - I    - Total number of nodes associated with the link
!  nbnocot          - IN     - I    - Total number of shell nodes at the interface
!  list_total_no_co - IN     - PTR  - List of total shell nodes at the interface
!  list_pairs       - IN     - PTR  - Array of pairs representing the shell-3D link
!                                     Dimensions: (2 * nb_pairs)
!  nb_pairs         - IN     - I    - Number of pairs in the link
!
!  OUTPUT PARAMETERS:
!  ------------------
!  map_noco_pair    - OUT   - I(*) - Mapping of each shell node to the pairs that contain the node.
!                                     Dimensions: (9, nbnocot, nb_pairs)
!  map_noco_nbnoco  - OUT   - I(*) - Mapping of each shell node to the number
!                                    of shell nodes in the pair. Dimensions: (9, nbnocot, nb_pairs)
!  map_noco_nbelem  - OUT   - I(*) - Mapping of each shell node to the total number of pairs that
!                                     contain the node.
!                                     Dimensions: (9, nbnocot)
! -------------------------------------------------------------------------------
!
! --------- VARIABLES LOCALES ---------------------------

    character(len=8)  :: typg_co_name, typg_vo_name
    character(len=16) :: typg_racc_name, typf_racc_name
    integer(kind=8) :: typg_racc_nume, typf_racc_nume
    integer(kind=8) :: el_co, el_vo, typg_co_nume, typg_vo_nume
    integer(kind=8) :: nb_nodes_co, nb_nodes_vo
    integer(kind=8) :: i, j, k, index, liel_l, nb_grel, elem, deca, jv_liel
    integer(kind=8), pointer :: v_list_type(:) => null()
    integer(kind=8), pointer :: v_mesh_typmail(:) => null()
    integer(kind=8), pointer :: v_connex(:) => null()
    integer(kind=8), pointer :: v_connex_lcum(:) => null()
    integer(kind=8), pointer :: v_ligr_nbno(:) => null()
    integer(kind=8), pointer :: v_index_bool(:) => null()
    integer(kind=8), pointer :: v_ligrel_nema(:) => null()
    integer(kind=8), pointer :: v_ligrel_liel(:) => null()
    integer(kind=8), pointer :: v_list_no_pair(:) => null()
    character(len=16):: nomte
    integer(kind=8):: ndim, nddl, nnco, nn3d
    character(len=8):: typmaco, typma3d

    ! TABLEAUX DE DONNEES
    integer(kind=8), parameter :: nb_racc = 8
    character(len=8), parameter, dimension(nb_racc) :: coq_el = (/ &
                                                       'SEG2    ', 'SEG2    ', 'SEG3    ', &
                                                       'SEG2    ', 'SEG2    ', 'SEG3    ', &
                                                       'SEG3    ', 'SEG3    '/)

    character(len=8), parameter, dimension(nb_racc) :: vol_el = (/ &
                                                       'TRIA3   ', 'QUAD4   ', 'TRIA6   ', &
                                                       'TRIA6   ', 'QUAD8   ', 'TRIA3   ', &
                                                       'QUAD4   ', 'QUAD8   '/)

    character(len=8), parameter, dimension(nb_racc) :: mesh_type = (/ &
                                                       'SE2TR3  ', 'SE2QU4  ', 'SE3TR6  ', &
                                                       'SE2TR6  ', 'SE2QU8  ', 'SE3TR3  ', &
                                                       'SE3QU4  ', 'SE3QU8  '/)

    character(len=8), parameter, dimension(nb_racc) :: fe_type = (/ &
                                                       'RACS2T3 ', 'RACS2Q4 ', 'RACS3T6 ', &
                                                       'RACS2T6 ', 'RACS2Q8 ', 'RACS3T3 ', &
                                                       'RACS3Q4 ', 'RACS3Q8 '/)
    integer(kind=8), parameter, dimension(nb_racc) :: nb_nodes = (/ &
                                                      5, 6, 9, &
                                                      8, 10, 6, &
                                                      7, 11/)

    AS_ALLOCATE(vi=v_index_bool, size=nb_racc)

    call jemarq()

! --- INITIALIZATION
    map_noco_nbelem = 0
    map_noco_pair = 0
    map_noco_nbnoco = 0

! --- MAILLAGE INFOS

    call jeveuo(noma//'.TYPMAIL', 'L', vi=v_mesh_typmail)
    call jeveuo(noma//'.CONNEX', 'L', vi=v_connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', vi=v_connex_lcum)
!
!   CREATION LIGREL
!
!   ALLOCATION
    AS_ALLOCATE(vi=v_list_type, size=nb_pairs)
!
! - PAS DE NOEUDS RETARDES
!
    call wkvect(ligrel//'.NBNO', 'V V I', 1, vi=v_ligr_nbno)
    v_ligr_nbno(1) = 0

!
! - OBJET NEMA
!
    call jecrec(ligrel//'.NEMA', 'V V I', 'NU', 'CONTIG', 'VARIABLE', nb_pairs)
    call jeecra(ligrel//'.NEMA', 'LONT', nt_nodes+nb_pairs)
    call jeveuo(ligrel//'.NEMA', 'E', vi=v_ligrel_nema)
    deca = 0
    do i = 1, nb_pairs
        el_co = list_pairs(2*(i-1)+1)
        el_vo = list_pairs(2*(i-1)+2)
        typg_co_nume = v_mesh_typmail(el_co)
        typg_vo_nume = v_mesh_typmail(el_vo)
        call jenuno(jexnum('&CATA.TM.NOMTM', typg_co_nume), typg_co_name)
        call jenuno(jexnum('&CATA.TM.NOMTM', typg_vo_nume), typg_vo_name)
        do j = 1, nb_racc
            if ((typg_co_name .eq. coq_el(j)) .and. (typg_vo_name .eq. vol_el(j))) then
                index = j
                exit
            end if
        end do
        v_list_type(i) = index
        v_index_bool(index) = v_index_bool(index)+1
        typg_racc_name = mesh_type(index)
        call jenonu(jexnom('&CATA.TM.NOMTM', typg_racc_name), typg_racc_nume)
        typf_racc_name = fe_type(index)
        !
        !- NOMBRE DES NOEUDS DE CHAQUE PARTIE
        !
        nb_nodes_co = v_connex_lcum(el_co+1)-v_connex_lcum(el_co)
        nb_nodes_vo = v_connex_lcum(el_vo+1)-v_connex_lcum(el_vo)
        ASSERT(nb_nodes(index) .eq. (nb_nodes_co+nb_nodes_vo))
        !
        !- CREER L'ELEMENT
        !
        call jecroc(jexnum(ligrel//'.NEMA', i))
        call jeecra(jexnum(ligrel//'.NEMA', i), 'LONMAX', ival=nb_nodes(index)+1)
        call jeecra(jexnum(ligrel//'.NEMA', i), 'LONUTI', ival=nb_nodes(index)+1)
        v_ligrel_nema(deca+nb_nodes(index)+1) = typg_racc_nume
        !
        ! - NUMEROS DES NOEUDS DES ELEMENTS DE LA PARTIE COQUE
        !
        do j = 1, nb_nodes_co
            v_ligrel_nema(deca+j) = v_connex(v_connex_lcum(el_co)-1+j)
        end do
        !
        ! - NUMEROS DES NOEUDS PARTIE 3D
        !
        do j = 1, nb_nodes_vo
            v_ligrel_nema(deca+nb_nodes_co+j) = v_connex(v_connex_lcum(el_vo)-1+j)
        end do
        deca = deca+nb_nodes(index)+1

    end do
!
! - OBJET LIEL
!
    liel_l = 0
    nb_grel = 0
    do index = 1, nb_racc
        if (v_index_bool(index) .gt. 0) then
            liel_l = liel_l+v_index_bool(index)+1
            nb_grel = nb_grel+1
        end if
    end do
    call jecrec(ligrel//'.LIEL', 'V V I', 'NU', 'CONTIG', 'VARIABLE', nb_grel)
    call jeecra(ligrel//'.LIEL', 'LONT', liel_l)

    elem = 0
    do index = 1, nb_racc
        jv_liel = 0
        if (v_index_bool(index) .gt. 0) then
            elem = elem+1
            call jecroc(jexnum(ligrel//'.LIEL', elem))
            call jeecra(jexnum(ligrel//'.LIEL', elem), 'LONMAX', ival=v_index_bool(index)+1)
            call jeecra(jexnum(ligrel//'.LIEL', elem), 'LONUTI', ival=v_index_bool(index)+1)
            call jeveuo(jexnum(ligrel//'.LIEL', elem), 'E', vi=v_ligrel_liel)
            typf_racc_name = fe_type(index)
            call jenonu(jexnom('&CATA.TE.NOMTE', typf_racc_name), typf_racc_nume)
            v_ligrel_liel(v_index_bool(index)+1) = typf_racc_nume
            ! recupérer les infos sur les elements
            nomte = typf_racc_name
            !
            call rco3d_elem(nomte, ndim, nddl, typmaco, nnco, typma3d, nn3d)
            !
            do i = 1, nb_pairs
                call jeveuo(jexnum(ligrel//'.NEMA', i), 'L', vi=v_list_no_pair)
                !
                if (v_list_type(i) .eq. index) then
                    jv_liel = jv_liel+1
                    v_ligrel_liel(jv_liel) = -i
                    ! 2D et 3D arrays pour optimiser les accès aux données
                    ! dans la subroutine raco3d lors de l'assemblage
                    do j = 1, nnco
                        do k = 1, nbnocot
                            if (list_total_no_co(k) .eq. v_list_no_pair(j)) then
                                map_noco_nbelem(elem, k) = map_noco_nbelem(elem, k)+1
                                map_noco_pair(elem, k, map_noco_nbelem(elem, k)) = jv_liel
                                map_noco_nbnoco(elem, k, map_noco_nbelem(elem, k)) = nnco
                                exit
                            end if
                        end do
                    end do
                    !
                end if
            end do
            ASSERT(jv_liel .eq. v_index_bool(index))
        end if
    end do

!
!- OBJECT LGRF
!
    call jedupo(mod(1:8)//'.MODELE    .LGRF', 'V', ligrel//'.LGRF', .false._1)
    !call adalig(ligrel)
    call initel(ligrel)

!-FIN

    AS_DEALLOCATE(vi=v_index_bool)
    AS_DEALLOCATE(vi=v_list_type)

    call jedema()

end subroutine
