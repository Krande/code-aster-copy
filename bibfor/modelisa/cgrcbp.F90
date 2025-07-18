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

subroutine cgrcbp(mofaz, iocc, nomaz, l_write, nbgraj)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jecroc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/lxlgut.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: iocc, nbgraj
    character(len=*) :: mofaz, nomaz
    aster_logical :: l_write
!
!       CGRCBP -- TRAITEMENT DE L'OPTION RELA_CINE_BP
!                 DU MOT FACTEUR CREA_GROUP_NO DE
!                 LA COMMANDE DEFI_GROUP
!
!      CETTE FONCTIONNALITE PERMET DE CREER POUR CHAQUE TRIPLET DE
!      LIAISONS (3 DIRECTIONS DE L'ESPACE) UN GROUP_NO CONSTITUE DES
!      NOEUDS IMPLIQUES DANS CE DERNIER. CELA EST FAIT A PARTIR DE LA
!      LISTE DE RELATION CONSTRUITE PAR DEFI_CABLE_BP
!
! -------------------------------------------------------
!  MOFAZ         - IN    - K16  - : MOT FACTEUR 'CREA_GROUP_NO'
!  IOCC          - IN    - I    - : NUMERO D'OCCURENCE DU MOT-FACTEUR
!  NOMAZ         - IN    - K8   - : NOM DU MAILLAGE
!  L_WRITE       - IN    - L    - : .TRUE. POUR AJOUTER LES GROUPES
!  NBGRAJ        - OUT   -  I   - : NOMBRE DE GROUPES AJOUTES
! -------------------------------------------------------
!
!.========================= DEBUT DES DECLARATIONS ====================
!
!
! --------- VARIABLES LOCALES ---------------------------
    character(len=8) :: cabl_prec
    integer(kind=8) :: nbrela, irela, ibid, lgpref, lennom, iad2, j, iret
    integer(kind=8) :: nb_coef, adr, i_coef, nbno_liai, nbno_max, nuno
    character(len=8) :: noma, prefix, nom_ddl
    character(len=16) :: motfac
    character(len=19) :: list_rela
    character(len=24) :: grpno, nomgno
!
    integer(kind=8), pointer :: rlnr(:) => null()
    integer(kind=8), pointer :: v_nb_coef(:) => null()
    integer(kind=8), pointer :: pointeur(:) => null()
    character(len=8), pointer :: v_nomnoe(:) => null()
    character(len=8), pointer :: v_ddl(:) => null()
    character(len=8), pointer :: listno(:) => null()
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ================
    motfac = mofaz
    noma = nomaz
    nbno_max = 0
    grpno = noma//'.GROUPENO       '
    nbgraj = 0
!
    call getvid(motfac, 'CABLE_BP', iocc=iocc, scal=cabl_prec, nbret=ibid)
    call getvtx(motfac, 'PREF_GRNO', iocc=iocc, scal=prefix, nbret=ibid)
    lgpref = lxlgut(prefix)
!
    list_rela = cabl_prec//'.LIRELA'
    call jeexin(list_rela//'.RLNR', iret)
    if (iret .eq. 0) then
        call utmess('F', 'CHARGES2_48', sk=cabl_prec)
    end if
!
!   nombre de relations
!
    call jeveuo(list_rela//'.RLNR', 'L', vi=rlnr)
    nbrela = rlnr(1)
    call jelibe(list_rela//'.RLNR')
!
!   nombre de coefficients
    call jeveuo(list_rela//'.RLNT', 'L', vi=v_nb_coef)
!
!   pointeur sur la liste
    call jeveuo(list_rela//'.RLPO', 'L', vi=pointeur)
!
!   nom des noeuds
    call jeveuo(list_rela//'.RLNO', 'L', vk8=v_nomnoe)
!
!   nom des ddls
    call jeveuo(list_rela//'.RLDD', 'L', vk8=v_ddl)
!
!   calcul de nbno_max
    do irela = 1, nbrela, 3
        nb_coef = v_nb_coef(irela)
        nbno_liai = 0
!       position du dernier terme de la relation
        adr = pointeur(irela)
        adr = adr-nb_coef
        do i_coef = 1, nb_coef
            nom_ddl = v_ddl(adr+i_coef)
            if (nom_ddl(2:2) .eq. 'R') then
                cycle
            end if
            nbno_liai = nbno_liai+1
        end do
        if (nbno_liai .gt. nbno_max) nbno_max = nbno_liai
        nbgraj = nbgraj+1
    end do
!
    if (l_write) then
        AS_ALLOCATE(vk8=listno, size=nbno_max)
!
        do irela = 1, nbrela, 3
            nb_coef = v_nb_coef(irela)
            nbno_liai = 0
!           position du dernier terme de la relation
            adr = pointeur(irela)
            adr = adr-nb_coef
            do i_coef = 1, nb_coef
                nom_ddl = v_ddl(adr+i_coef)
                if (nom_ddl(2:2) .eq. 'R') then
                    cycle
                end if
                nbno_liai = nbno_liai+1
                listno(nbno_liai) = v_nomnoe(adr+i_coef)
            end do
            ASSERT(nbno_liai .gt. 1)
!
!           ajout du groupe
!
            lennom = lxlgut(listno(1))
            nomgno = prefix(1:lgpref)//listno(1) (1:lennom)
            call jecroc(jexnom(grpno, nomgno))
            call jeecra(jexnom(grpno, nomgno), 'LONMAX', nbno_liai)
            call jeecra(jexnom(grpno, nomgno), 'LONUTI', nbno_liai)
            call jeveuo(jexnom(grpno, nomgno), 'E', iad2)
            do j = 1, nbno_liai
                nuno = char8_to_int(listno(j))
                zi(iad2-1+j) = nuno
            end do
        end do
        AS_DEALLOCATE(vk8=listno)
    end if
!
    call jedema()
!.============================ FIN DE LA ROUTINE ======================
end subroutine
