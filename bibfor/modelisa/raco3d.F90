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
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/digdel.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeexin.h"
#include "asterfort/mecact.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/jemarq.h"
#include "asterfort/jexnum.h"
#include "asterfort/rco3d_apco3d.h"
#include "asterfort/rco3d_crep.h"
#include "asterfort/rco3d_crealigrel.h"
#include "asterfort/rco3d_addrela.h"
!
    integer :: iocc
    character(len=8) :: charge
    character(len=14) :: numddl
    character(len=19) :: lisrel
    character(len=*) :: numdlz, chargz, fonrez, lisrez
! -------------------------------------------------------
!     LIAISON COQUE-3D PAR DES RELATIONS LINEAIRES
!     ENTRE LES MAILLES DE BORDS
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
    character(len=19) :: ligrmo, ligrel, chmlrac
    character(len=24) :: lismaco, lismavo, lisnoco
    character(len=8)  :: mod, noma, licmp(6)
    integer :: nbmavo, nbmaco, nt_nodes
    character(len=8), pointer :: lgrf(:) => null()
    integer :: nb_pairs, iret
    integer :: i, j, k, l, index, n1, elem
    real(kind=8) :: epai, icmp(6), crig
    integer, pointer :: list_pairs(:) => null()
    character(len=8) :: lpain(2), lpaout(1)
    character(len=24) :: lchin(2), lchout(1), valech
    integer :: nbnocot, jlisnoco
    integer, allocatable :: map_noco_pair(:, :, :)
    integer, allocatable :: map_noco_nbnoco(:, :, :)
    integer, allocatable :: map_noco_nbelem(:, :)
    integer, pointer :: list_total_no_co(:) => null()
    character(len=8) :: cara_elem

    call jemarq()

    motcle(1) = 'GROUP_MA_1'
    motcle(2) = 'GROUP_MA_2'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'GROUP_MA'

    numddl = numdlz
    charge = chargz
    lisrel = lisrez

    lismavo = '&&RACO3D.LMAILLES.VOL'
    lisnoco = '&&RACO3D.LNOEUDS.COQ'
    lismaco = '&&RACO3D.LMAILLES.COQ'
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

!--- RECUPERER EPAISSEUR

    call getvr8('LIAISON_ELEM', 'EPAISSEUR', iocc=iocc, scal=epai, nbret=n1)

    if (n1 .eq. 0) then
        ASSERT(.false.)
    end if

!--- RECUPERER COEF_RIGI_DRZ

    call getvr8('LIAISON_ELEM', 'COEF_RIGI_DRZ', iocc=iocc, scal=crig, nbret=n1)

    if (n1 .eq. 0) then
        crig = 1.d0-6
    end if
!--- -----------------------------------
!-- RECUPERER LA LISTE DES MAILLES

    call reliem(' ', noma, 'NU_MAILLE', motfac, iocc, &
                1, motcle(1), typmcl(1), lismaco, nbmaco)

    call reliem(' ', noma, 'NU_MAILLE', motfac, iocc, &
                1, motcle(2), typmcl(2), lismavo, nbmavo)

!- RECUPERER LA LISTE DES NOOEUDS DU BORD DE LA COQUE

    call reliem(' ', noma, 'NU_NOEUD', motfac, iocc, &
                1, motcle(1), typmcl(1), lisnoco, nbnocot)
    call jeveuo(lisnoco, 'L', jlisnoco)
    !
    AS_ALLOCATE(vi=list_total_no_co, size=nbnocot)
    !
    do i = 1, nbnocot
        list_total_no_co(i) = zi(jlisnoco-1+i)
    end do

!-- RECUPERER LA LISTE DES PAIRES
    nb_pairs = 0
    nt_nodes = 0

    call rco3d_apco3d(noma, lismavo, lismaco, nbmavo, nbmaco, epai, &
                list_pairs, nb_pairs, nt_nodes)
!

!-- CONSTRUCTION DU LIGREL
!
!   2D ET 3D ARRAYs POUR ACCELERER L ACCES AUX DONNEES AU MOMENT
!   DE L ASSEMBLAGE  DES MATRICES
                
    allocate (map_noco_pair(9, nbnocot, nb_pairs))
    allocate (map_noco_nbnoco(9, nbnocot, nb_pairs))
    allocate (map_noco_nbelem(9, nbnocot))
    !
    call rco3d_crealigrel(ligrel, noma, mod, list_pairs, &
                           nb_pairs, nt_nodes, &
                           list_total_no_co, nbnocot, map_noco_pair, &
                           map_noco_nbelem, map_noco_nbnoco)
!   LIRE EPAISSEUR
                           
    call getvid('LIAISON_ELEM', 'CARA_ELEM', iocc=iocc, scal=cara_elem, nbret=n1)
    chmlrac = '&&RACO3D.PCACOQU.CM'
    call alchml(ligrel, 'LIAI_CO_3D', 'PCACOQU', 'V', chmlrac, iret, ' ')
    call rco3d_crep(ligrel, noma, chmlrac, cara_elem, epai)

!    CREATION DE LA CARTE DES CARACTERISTIQUES DE LA PARTIE COQUE
!
    licmp(1) = 'EP'
    licmp(2) = 'ALPHA'
    licmp(3) = 'BETA'
    licmp(4) = 'CTOR'
    licmp(5) = 'EXCENT'
    licmp(6) = 'INERTIE'
    icmp(1) = epai
    icmp(2) = 0.0d0
    icmp(3) = 0.0d0
    icmp(4) = crig
    icmp(5) = 0.0d0
    icmp(6) = 0.0d0
    !
    call mecact('V', '&&RACO3D.PCACOQU', 'LIGREL', ligrel, 'CACOQU_R', &
                ncmp=6, lnomcmp=licmp, vr=icmp)

!--  Fields
!
    lpain(1) = 'PGEOMER'
    lpain(2) = 'PCACOQU'
    lchin(1) = noma//'.COORDO'
    lchin(2) = '&&RACO3D.PCACOQU.CM'
    lpaout(1) = 'PMATUNS'
    lchout(1) = '&&RACO3D.PMATUNS'

!-- Compute elementary matrices
    call calcul('S', 'LIAI_CO_3D', ligrel, 2, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')

!-- add the linear relations
    call rco3d_addrela(ligrel, noma, nb_pairs, nbnocot, &
                     list_total_no_co, map_noco_pair, map_noco_nbelem, &
                     map_noco_nbnoco, lchout(1) (1:19), fonrez, lisrel)


!FIN
    call detrsd('CHAM_ELEM', chmlrac)
    call detrsd('LIGREL', ligrel)
    call detrsd('CARTE', '&&RACO3D.PCACOQU')
    call detrsd('RESUELEM', '&&RACO3D.PMATUNS')
    call jedetr('&&RACO3D.LMAILLES.VOL')
    call jedetr('&&RACO3D.LMAILLES.COQ')
    call jedetr('&&RACO3D.LNOEUDS.COQ')
    AS_DEALLOCATE(vi=list_pairs)
    AS_DEALLOCATE(vi=list_total_no_co)
    deallocate (map_noco_pair)
    deallocate (map_noco_nbelem)
    deallocate (map_noco_nbnoco)

    call jedema()


end subroutine
