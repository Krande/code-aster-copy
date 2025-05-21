! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! aslint: disable=W1501
!
subroutine raco3d(numdlz, iocc, fonrez, lisrez, chargz)
    !
    implicit none
    !
    !    ATTENTION CETTE PROGRAMMATION SUPPOSE QUE L'OBJET NUEQ EST UN
    !    VECTEUR IDENTITE. A MODIFIER
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/apco3d.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mesomm.h"
#include "asterfort/nbelem.h"
!#include "asterfort/pacoor.h"
#include "asterfort/reliem.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"

#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
!
    integer :: iocc
    character(len=8) :: charge
    character(len=14) :: numddl
    character(len=19) :: lisrel
    character(len=*) :: numdlz, chargz, fonrez, lisrez
! -------------------------------------------------------
!     RACCO 3D-COQUE PAR DES RELATIONS LINEAIRES
!     ENTRE LES NOEUDS DES MAILLES DE BORD MODELISANT
!     LA TRACE DE LA SECTION DE LA POUTRE SUR LA COQUE
!     ET LE NOEUD DE LA POUTRE DONNE PAR L'UTILISATEUR
! -------------------------------------------------------
!  NUMDLZ        - IN    - K14  - : NOM DU NUMEDDL DU LIGREL DU MODELE
!                                     (IL A ETE CREE SUR LA VOLATILE)
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  FONREZ        - IN    - K4   - : 'REEL'
!  LISREZ        - IN    - K19  - : NOM DE LA SD LISTE_RELA
!  CHARGE        - IN    - K8   - : NOM DE LA SD CHARGE
!                - JXVAR -      -
! -------------------------------------------------------
!
!
! --------- VARIABLES LOCALES ---------------------------

    character(len=16) :: motfac, motcle(2), typmcl(2)
    character(len=19) :: ligrmo, ligrel
    character(len=24) :: lismaco, lismavo
    character(len=8)  :: mod, noma
    character(len=8)  :: typg_co_name, typg_vo_name
    character(len=16) :: typg_racc_name, typf_racc_name
    integer :: nbmavo, nbmaco, nt_nodes
    character(len=8), pointer :: lgrf(:) => null()
    integer :: nb_pairs, el_co, el_vo, typg_co_nume, typg_vo_nume
    integer :: i, j, index, deca
    integer, pointer :: list_pairs(:) => null()
    integer, pointer :: v_list_type(:)=> null()
    integer, pointer :: v_mesh_typmail(:) => null()
    integer, pointer :: v_ligr_nbno(:) => null()
    integer, pointer :: v_index_bool(:) => null()
    
    ! TABLEAUX DE DONNEES
    integer, parameter :: nb_racc = 2
    character(len=8), parameter, dimension(nb_racc) :: coq_el = (/&
                                'SEG2    ', 'SEG2    ' /)
    
    character(len=8), parameter, dimension(nb_racc) :: vol_el = (/&
                                'TRIA3   ', 'QUAD4   ' /)
    
    character(len=8), parameter, dimension(nb_racc) :: mesh_type = (/&
                                'SE2TR3  ', 'SE2QU4  ' /)
    
    character(len=8), parameter, dimension(nb_racc) :: fe_type = (/&
                                'RACS2T3 ', 'RACS2Q4 ' /)
    integer, parameter, dimension(nb_racc) :: nb_node = (/ &
                                5, 6 /)

    AS_ALLOCATE(vi=v_index_bool, size=nb_racc)


    !character(len=16) :: typf_cont_name
    !integer :: typf_cont_nume

    motcle(1) = 'GROUP_MA_1'
    motcle(2) = 'GROUP_MA_2'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'GROUP_MA'

    numddl = numdlz
    charge = chargz
    lisrel = lisrez

    lismavo = '&&RACO3D.LISTE_MAILLES.VOL'
    lismaco = '&&RACO3D.LISTE_MAILLES.COQ'
    motfac = 'LIAISON_ELEM'
    ligrel = '&&RACO3D'



! --- MODELE ASSOCIE AU LIGREL DE CHARGE :
!     ----------------------------------
    call dismoi('NOM_MODELE', charge(1:8), 'CHARGE', repk=mod)

! ---  LIGREL DU MODELE :
!      ----------------
    ligrmo = mod(1:8)//'.MODELE'
!
! --- MAILLAGE ASSOCIE AU MODELE :
!     --------------------------
    call jeveuo(ligrmo//'.LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)
! --- TYPE MAILLE
    call jeveuo(noma//'.TYPMAIL', 'L', vi=v_mesh_typmail)

!--- -----------------------------------
!-- RECUPERER LA LISTE DES MAILLES

    call reliem(' ', noma, 'NU_MAILLE', motfac, iocc, &
                1, motcle(1), typmcl(1), lismaco, nbmaco)

    call reliem(' ', noma, 'NU_MAILLE', motfac, iocc, &
                1, motcle(2), typmcl(2), lismavo, nbmavo)

!-- RECUPERER LA LISTE DES PAIRES
    nb_pairs = 0
    nt_nodes = 0
    call apco3d(noma, lismavo, lismaco, nbmavo, nbmaco, &
                list_pairs, nb_pairs, nt_nodes)
    
!-- CONSTRUCTION DU LIGREL
    
    
    !typf_cont_name = 'RACS2T3 '
    !call jenonu(jexnom('&CATA.TE.NOMTE', typf_cont_name), typf_cont_nume)
    !write(*,*) "#########################  ", typf_cont_nume

!   ALLOCATION
    AS_ALLOCATE(vi=v_list_type, size=nb_pairs)
!
! - NO LATE NODES
!
    call wkvect(ligrel//'.NBNO', 'V V I', 1, vi=v_ligr_nbno)
    v_ligr_nbno(1) = 0

!
! - OBJET NEMA
!
    call jecrec(ligrel//'.NEMA', 'V V I', 'NU', 'CONTIG', 'VARIABLE', nb_pairs)
    call jeecra(ligrel//'.NEMA', 'LONT', nt_nodes+nb_pairs)
    deca = 0
    do i = 1, nb_pairs
        el_co = list_pairs(2*(i-1) + 1)
        el_vo = list_pairs(2*(i-1) + 2)
        typg_co_nume = v_mesh_typmail(el_co)
        typg_vo_nume = v_mesh_typmail(el_vo)
        call jenuno(jexnum('&CATA.TM.NOMTM', typg_co_nume), typg_co_name)
        call jenuno(jexnum('&CATA.TM.NOMTM', typg_vo_nume), typg_vo_name)
        do j = 1, nb_racc
            if ( (typg_co_name .eq. coq_el(j)) .and. (typg_vo_name .eq. vol_el(j)) ) then
                    index = j
                    exit
            end if
        end do
        v_list_type(i) = index
        v_index_bool(index) = 1
        typg_racc_name = mesh_type(index)
        typf_racc_name = fe_type(index)
    end do

!FIN 
    AS_DEALLOCATE(vi=list_pairs)
    AS_DEALLOCATE(vi=v_index_bool)
    AS_DEALLOCATE(vi=v_list_type)

end subroutine