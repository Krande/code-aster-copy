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

subroutine nugeel(nugene, modgen)
!    M.CORUS     DATE 23/02/10
!-----------------------------------------------------------------------
!  BUT:      < NUMEROTATION GENERALISEE AVEC ELIMINATION>
!
!  ON REMPLIT LE NUME_DDL COMME S'IL N'Y AVAIT QU'UNE SEULE SOUS
!  STRUCTURE. LES INFOS POUR LA RESTITUTION SONT DANS
!     SELIAI   : MATRICE DE PROJECTION
!     SIZLIA  : NB DE DDL POUR CHAQUE SOUS STRUCTURE
!     SST      : NOMS DES SOUS STRUCTURES, DANS L'ORDRE D'ASSEMBLAGE
!
    implicit none
!
!  DETERMINER LA NUMEROTATION DES DEGRES DE LIBERTE GENERALISES
!   A PARTIR D'UN MODELE GENERALISE
!
!-----------------------------------------------------------------------
!
! NUGENE   /I/: NOM K14 DU NUME_DDL_GENE
! MODGEN   /I/: NOM K8 DU MODELE GENERALISE
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/indlia.h"
#include "asterfort/iunifi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nume_equa_gene_crsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
!
!
    integer(kind=8) :: nindep, nbddl, nb_sstr
    character(len=6) :: pgc
    character(len=8) :: modgen, sst1, sst2
    character(len=14) :: nugene
    character(len=19) :: nume_equa_gene
    character(len=24) :: defli, nomsst, sizlia, sst
    character(len=24) :: valk, seliai
    aster_logical :: pbcone
!
!
!---------- VARIABLES PERSOS -------------------------------------------
!
!
    integer(kind=8) :: lsilia, lsst, i1, j1, imes, nblia, nbsst
    integer(kind=8) :: icomp, ltssnb, ltlia, nulia, i_ligr_sstr
    integer(kind=8) :: nusst1, nusst2, nusst, ltsst, lldefl
    integer(kind=8), pointer :: prgene_orig(:) => null()
    integer(kind=8), pointer :: prgene_prno(:) => null()
!
!
!-----------------------------------------------------------------------
    data pgc/'NUMGEN'/
!-----------------------------------------------------------------------
!
    call jemarq()
    imes = iunifi('MESSAGE')
!
!-----------------------------------------------------------------------
!
    defli = modgen//'      .MODG.LIDF'
    nomsst = modgen//'      .MODG.SSNO'
!
!----------------------RECUPERATION DES DIMENSIONS PRINCIPALES----------
!
    call jelira(defli, 'NMAXOC', nblia)
    call jelira(nomsst, 'NOMMAX', nbsst)

!  ON REMPLIT LE NUME_DDL COMME S'IL N'Y AVAIT QU'UNE SEULE SOUS
!  STRUCTURE.
    nb_sstr = 1

!
!
!
!----------------------------------------------------------------C
!--                                                            --C
!--  CONSTRUCTION D'UNE MATRICE REGROUPANT TOUTES LES LIAISONS --C
!--      ET TEST DE l'INDEPENDANCE DES RELATIONS LINEAIRES     --C
!--                                                            --C
!----------------------------------------------------------------C
!
    seliai = nugene(1:14)//'.ELIM.BASE'
    sizlia = nugene(1:14)//'.ELIM.TAIL'
    sst = nugene(1:14)//'.ELIM.NOMS'
!
    call indlia(modgen, seliai, nindep, nbddl, sst, &
                sizlia)
!
! - Create nume_equa_gene
!
    nume_equa_gene = nugene//'.NUME'
    call nume_equa_gene_crsd(nume_equa_gene, 'G', nindep, nb_sstr=nb_sstr, nb_link=0, &
                             model_genez=modgen, gran_namez='DEPL_R')
!
! - Set sub_structures
!
    call jenonu(jexnom(nume_equa_gene//'.LILI', '&SOUSSTR'), i_ligr_sstr)
    ASSERT(i_ligr_sstr .eq. 1)
    call jeveuo(jexnum(nume_equa_gene//'.PRNO', i_ligr_sstr), 'E', vi=prgene_prno)
    call jeveuo(jexnum(nume_equa_gene//'.ORIG', i_ligr_sstr), 'E', vi=prgene_orig)
    prgene_prno(1) = 1
    prgene_prno(2) = nindep
    prgene_orig(1) = 1

!
!----------------------BOUCLES DE COMPTAGE DES DDL----------------------
!
    icomp = 0

!
    call jeveuo(sizlia, 'L', lsilia)
    call jeveuo(sst, 'L', lsst)
!
    write (imes, *) '+++ NOMBRE DE SOUS-STRUSTURES: ', nbsst
    write (imes, *) '+++ NOMBRE DE LIAISONS: ', nblia
!
!
!
!------------------------ALLOCATIONS DIVERSES---------------------------
!
!
    call wkvect('&&'//pgc//'.SST.NBLIA', 'V V I', nbsst, ltssnb)
!
    call wkvect('&&'//pgc//'.LIA.SST', 'V V I', nblia*2, ltlia)
    call jecrec('&&'//pgc//'.SST.LIA', 'V V I', 'NU', 'DISPERSE', 'CONSTANT', &
                nbsst)
    call jeecra('&&'//pgc//'.SST.LIA', 'LONMAX', 2*nblia)
!
!
!   BOUCLE DE DETERMINATION DE LA RELATION
!   NUMERO TARDIF  LIAISON --> NUMERO SOUS-STRUCTURE DE PLUS PETIT
!                              NUMERO
!
!   ON CONSERVE POUR DETECTER LES SOUS STRUCTURES NON CONNECTEES
!
!
    do i1 = 1, nblia*2
        nulia = int((i1-1)/2)+1
        call jeveuo(jexnum(defli, nulia), 'L', lldefl)
        sst1 = zk8(lldefl)
        sst2 = zk8(lldefl+2)
        call jenonu(jexnom(nomsst, sst1), nusst1)
        call jenonu(jexnom(nomsst, sst2), nusst2)
!
        zi(ltssnb+nusst1-1) = 1
        zi(ltssnb+nusst2-1) = 1
        zi(ltlia+i1-1) = max(nusst1, nusst2)
    end do
!
!   BOUCLE PERMETTANT DE DETERMINER L'INVERSE
!   NUMERO TARDIF  SOUS-STRUCTURE --> NUMEROS TARDIF LIAISONS
!                     DONT ELLE EST LA STRUCTURE DE PLUS PETIT NUMERO
!
!   ET POUR DETECTER LES SOUS-STRUCTURES NON CONNECTEES
!
    pbcone = .false.
    do i1 = 1, nbsst
        icomp = 0
        call jenonu(jexnom(nomsst, zk8(lsst+i1-1)), nusst)
        if (zi(ltssnb+nusst-1) .eq. 0) then
            pbcone = .true.
            call jenuno(jexnum(nomsst, nusst), sst1)
            valk = sst1
            call utmess('E', 'ALGORITH13_75', sk=valk)
        end if
        call jecroc(jexnum('&&'//pgc//'.SST.LIA', i1))
        call jeveuo(jexnum('&&'//pgc//'.SST.LIA', i1), 'E', ltsst)
        do j1 = 1, nblia*2
            if (zi(ltlia+j1-1) .eq. nusst) then
                icomp = icomp+1
                zi(ltsst+icomp-1) = j1
            end if
        end do
    end do
!
    if (pbcone) then
        call utmess('F', 'ALGORITH13_76')
    end if
!
    call jedetr('&&'//pgc//'.LIA.SST')
    call jedetr('&&'//pgc//'.SST.NBLIA')
!
!--------------------REMPLISSAGE DES NUMERO D'EQUATION-----------------
!
!
!
    call jedema()
end subroutine
