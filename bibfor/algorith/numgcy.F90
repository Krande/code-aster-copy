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

subroutine numgcy(nugene, modgen)
!    O. NICOLAS
!-----------------------------------------------------------------------
!  BUT:      < NUMEROTATION GENERALISEE >
    implicit none
!
!  DETERMINER LA NUMEROTATION DES DEGRES DE LIBERTE GENERALISES
!   A PARTIR D'UN MODELE GENERALISE CAS CYCLIQUE POUR DESACORDAGE
!
!-----------------------------------------------------------------------
!
! NUGENE   /I/: NOM K14 DU NUME_DDL_GENE
! MODGEN   /I/: NOM K8 DU MODELE GENERALISE
!
!
!
#include "jeveux.h"
#include "asterfort/crsmos.h"
#include "asterfort/assert.h"
#include "asterfort/iunifi.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nume_equa_gene_crsd.h"
#include "asterfort/mgutdm.h"
!
!
    character(len=8) :: modgen, nomcou, kbid
    character(len=14) :: nugene
    character(len=19) :: nume_equa_gene, stomor
    character(len=24) :: defli, fprofl, nomsst
    integer(kind=8) :: ibid, i, i_ligr_link, nb_link, nb_sstr, i_ligr_sstr
    character(len=24) :: lili, prno, orig
    integer(kind=8), pointer :: prgene_orig(:) => null()
    integer(kind=8), pointer :: prgene_prno(:) => null()

!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icompl, icomps, ifimes, llprof, nblia
    integer(kind=8) :: nblig, nbmod, nbsst, neq
    integer(kind=8), pointer :: mael_raid_desc(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    ifimes = iunifi('MESSAGE')
!-----------------------------------------------------------------------
!
    kbid = ' '
    defli = modgen//'      .MODG.LIDF'
    fprofl = modgen//'      .MODG.LIPR'
    nomsst = modgen//'      .MODG.SSNO'
    nume_equa_gene = nugene//'.NUME'
    stomor = nugene//'.SMOS'
    lili = nume_equa_gene//'.LILI'
    prno = nume_equa_gene//'.PRNO'
    orig = nume_equa_gene//'.ORIG'

! ON RECUPERE LE NOMBRE DE LIAISON
    call jelira(defli, 'NMAXOC', nblia)
! ON RECUPERE LE NOMBRE DE SOUS-STRUCTURE
    call jelira(nomsst, 'NOMMAX', nbsst)
!
!----------------------BOUCLES DE COMPTAGE DES DDL----------------------
!
! ICOMPS EST LE NOMBRE TOTAL DE MODES DANS LES SOUS-STRUCTURES
    icomps = 0
! ICOMPS EST LE NOMBRE TOTAL DE MODES D'INTERFACE DANS LES
! SOUS-STRUCTURES
    icompl = 0
!
!   BOUCLE SUR LES SOUS-STRUCTURES
!
    do i = 1, nbsst
        call mgutdm(modgen, kbid, i, 'NOM_MACR_ELEM', ibid, &
                    nomcou)
        call jeveuo(nomcou//'.MAEL_RAID_DESC', 'L', vi=mael_raid_desc)
        nbmod = mael_raid_desc(2)
        icomps = icomps+nbmod
    end do
!
!   BOUCLE SUR LES LIAISONS
!
    call jeveuo(fprofl, 'L', llprof)
    do i = 1, nblia
        nblig = zi(llprof+(i-1)*9)
        icompl = icompl+nblig
    end do
!
    neq = icomps-icompl
!
    write (ifimes, *) '+++ NOMBRE DE SOUS-STRUCTURES: ', nbsst
    write (ifimes, *) '+++ NOMBRE DE LIAISONS: ', nblia
    write (ifimes, *) '+++ NOMBRE TOTAL D''EQUATIONS: ', neq
    write (ifimes, *) '+++ DONT NOMBRE D''EQUATIONS STRUCTURE: ', icomps
    write (ifimes, *) '+++ DONT NOMBRE D''EQUATIONS LIAISON: ', icompl
!  ON REMPLIT LE NUME_DDL COMME S'IL N'Y AVAIT QU'UNE SEULE SOUS
!  STRUCTURE.
    nb_sstr = 1
    nb_link = 1
!
! - Create nume_equa_gene
!
    call nume_equa_gene_crsd(nume_equa_gene, 'G', neq, nb_sstr=nb_sstr, nb_link=nb_link, &
                             model_genez=modgen, gran_namez='DEPL_R')
!
! - Set sub_structures
!
    call jenonu(jexnom(lili, '&SOUSSTR'), i_ligr_sstr)
    ASSERT(i_ligr_sstr .eq. 1)
    call jeveuo(jexnum(prno, i_ligr_sstr), 'E', vi=prgene_prno)
    call jeveuo(jexnum(orig, i_ligr_sstr), 'E', vi=prgene_orig)
    prgene_prno(1) = 1
    prgene_prno(2) = neq
    prgene_orig(1) = 1
!
! - Set links
!
    call jenonu(jexnom(lili, 'LIAISONS'), i_ligr_link)
    call jeveuo(jexnum(prno, i_ligr_link), 'E', vi=prgene_prno)
    call jeveuo(jexnum(orig, i_ligr_link), 'E', vi=prgene_orig)
    prgene_prno(1) = 0
    prgene_orig(1) = 1

!
!     CREATION DU STOCKAGES MORSE :
    call crsmos(stomor, 'PLEIN', neq)
!
!
    call jedema()
end subroutine
