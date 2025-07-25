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

subroutine numgen(nugene, modgen)
!    P. RICHARD     DATE 13/10/92
!-----------------------------------------------------------------------
!  BUT:      < NUMEROTATION GENERALISEE >
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
#include "asterfort/mgutdm.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/nume_equa_gene_crsd.h"
!
!
!
    character(len=6) :: pgc
    character(len=8) :: modgen, nomcou, sst1, sst2, kbid
    character(len=14) :: nugene
    character(len=19) :: nume_equa_gene
    character(len=24) :: defli, fprofl, nomsst
    character(len=24) :: valk
    aster_logical :: assok, pbcone
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid, icomp, icompl, icomps, ifimes, i_ligr_sstr, i_ligr_link
    integer(kind=8) :: j, k, lddeeq, ldnueq
    integer(kind=8) :: lldefl, lldesc, llprof, nb_sstr, nb_link
    integer(kind=8) :: ltlia, ltoras, ltssnb, ltsst, nbddl, nblia, nblig
    integer(kind=8) :: nbmod, neq, ntail, nuas, nulia, null2
    integer(kind=8) :: nult, nusst, nusst1, nusst2, nut, nutarl
    character(len=24) :: lili, prno, orig, deeq, nueq
    integer(kind=8), pointer :: nume_equa_gene_orig_s(:) => null()
    integer(kind=8), pointer :: nume_equa_gene_orig_l(:) => null()
    integer(kind=8), pointer :: nume_equa_gene_prno_s(:) => null()
    integer(kind=8), pointer :: nume_equa_gene_prno_l(:) => null()

!-----------------------------------------------------------------------
    data pgc/'NUMGEN'/
!-----------------------------------------------------------------------
!
    call jemarq()
    ifimes = iunifi('MESSAGE')
!
!-----------------------------------------------------------------------
!
    defli = modgen//'      .MODG.LIDF'
    fprofl = modgen//'      .MODG.LIPR'
    nomsst = modgen//'      .MODG.SSNO'
!
!--------------------CREATION DU .REFN----------------------------------
!                       ET DU DESC
    nume_equa_gene = nugene//'.NUME'
    lili = nume_equa_gene//'.LILI'
    prno = nume_equa_gene//'.PRNO'
    orig = nume_equa_gene//'.ORIG'
    deeq = nume_equa_gene//'.DEEQ'
    nueq = nume_equa_gene//'.NUEQ'
!
!----------------------RECUPERATION DES DIMENSIONS PRINCIPALES----------
!
    call jelira(defli, 'NMAXOC', nblia)
    call jelira(nomsst, 'NOMMAX', nb_sstr)

!
!----------------------BOUCLES DE COMPTAGE DES DDL----------------------
!
    icomp = 0
    icomps = 0
    icompl = 0
!
    call jeveuo(fprofl, 'L', llprof)
!
!   BOUCLE SUR LES SOUS-STRUCTURES
!
    do i = 1, nb_sstr
        kbid = '        '
        call mgutdm(modgen, kbid, i, 'NOM_MACR_ELEM', ibid, &
                    nomcou)
        call jeveuo(nomcou//'.MAEL_RAID_DESC', 'L', lldesc)
        nbmod = zi(lldesc+1)
        icomp = icomp+nbmod
        icomps = icomps+nbmod
    end do
!
!   BOUCLE SUR LES LIAISONS
!   (ON SUPPOSE QUE LES MATRICES DES LIAISONS 1 ET 2 ONT
!   MEME NOMBRE DE LIGNES = VERIF VERILI)
!
    do i = 1, nblia
        nblig = zi(llprof+(i-1)*9)
        icomp = icomp+2*nblig
        icompl = icompl+2*nblig
    end do
!
    neq = icomp
!
    write (ifimes, *) '+++ NOMBRE DE SOUS-STRUCTURES: ', nb_sstr
    write (ifimes, *) '+++ NOMBRE DE LIAISONS: ', nblia
    write (ifimes, *) '+++ NOMBRE TOTAL D''EQUATIONS: ', neq
    write (ifimes, *) '+++ DONT NOMBRE D''EQUATIONS STRUCTURE: ', icomps
    write (ifimes, *) '+++ DONT NOMBRE D''EQUATIONS LIAISON: ', icompl
!
! - Create NUME_EQUA_GENE
!
    nb_sstr = nb_sstr
    nb_link = 2*nblia
    call nume_equa_gene_crsd(nume_equa_gene, 'G', neq, nb_sstr=nb_sstr, nb_link=nb_link, &
                             model_genez=modgen, gran_namez='DEPL_R')
    call jeveuo(deeq, 'E', lddeeq)
    call jeveuo(nueq, 'E', ldnueq)
!
! - Set LIGREL for substructuring
!
    call jenonu(jexnom(lili, '&SOUSSTR'), i_ligr_sstr)
    call jeveuo(jexnum(prno, i_ligr_sstr), 'E', vi=nume_equa_gene_prno_s)
    call jeveuo(jexnum(orig, i_ligr_sstr), 'E', vi=nume_equa_gene_orig_s)
!
    do i = 1, nb_sstr
        kbid = '        '
        call mgutdm(modgen, kbid, i, 'NOM_MACR_ELEM', ibid, &
                    nomcou)
        call jeveuo(nomcou//'.MAEL_RAID_DESC', 'L', lldesc)
        nbmod = zi(lldesc+1)
        nume_equa_gene_orig_s(i) = i
        nume_equa_gene_prno_s((i-1)*2+2) = nbmod
    end do
!
! - Add LIGREL LIAISONS
!
    call jenonu(jexnom(lili, 'LIAISONS'), i_ligr_link)
    call jeveuo(jexnum(prno, i_ligr_link), 'E', vi=nume_equa_gene_prno_l)
    call jeveuo(jexnum(orig, i_ligr_link), 'E', vi=nume_equa_gene_orig_l)

!
    do i = 1, nblia
        nblig = zi(llprof+(i-1)*9)
        nume_equa_gene_orig_l((i-1)*2+1) = i
        nume_equa_gene_orig_l((i-1)*2+2) = i
        nume_equa_gene_prno_l((i-1)*4+2) = nblig
        nume_equa_gene_prno_l((i-1)*4+4) = nblig
    end do
!
    call wkvect('&&'//pgc//'.SST.NBLIA', 'V V I', nb_sstr, ltssnb)
    call wkvect('&&'//pgc//'.LIA.SST', 'V V I', nblia*2, ltlia)
    call jecrec('&&'//pgc//'.SST.LIA', 'V V I', 'NU', 'DISPERSE', 'CONSTANT', &
                nb_sstr)
    call jeecra('&&'//pgc//'.SST.LIA', 'LONMAX', 2*nblia)
!
!   BOUCLE DE DETERMINATION DE LA RELATION
!   NUMERO TARDIF  LIAISON --> NUMERO SOUS-STRUCTURE DE PLUS PETIT
!                              NUMERO
!
    do i = 1, nblia*2
        nulia = nume_equa_gene_orig_l(1+i-1)
        call jeveuo(jexnum(defli, nulia), 'L', lldefl)
        sst1 = zk8(lldefl)
        sst2 = zk8(lldefl+2)
        call jenonu(jexnom(nomsst, sst1), nusst1)
        call jenonu(jexnom(nomsst, sst2), nusst2)
        zi(ltssnb+nusst1-1) = 1
        zi(ltssnb+nusst2-1) = 1
        zi(ltlia+i-1) = max(nusst1, nusst2)
    end do
!
!   BOUCLE PERMETTANT DE DETERMINER L'INVERSE
!   NUMERO TARDIF  SOUS-STRUCTURE --> NUMEROS TARDIF LIAISONS
!                     DONT ELLE EST LA STRUCTURE DE PLUS PETIT NUMERO
!
!   ET POUR DETECTER LES SOUS-STRUCTURES NON CONNECTEES
!
    pbcone = .false.
    do i = 1, nb_sstr
        icomp = 0
        nusst = nume_equa_gene_orig_s(i)
        if (zi(ltssnb+nusst-1) .eq. 0) then
            pbcone = .true.
            call jenuno(jexnum(nomsst, nusst), sst1)
            valk = sst1
            call utmess('E', 'ALGORITH13_75', sk=valk)
        end if
        call jecroc(jexnum('&&'//pgc//'.SST.LIA', i))
        call jeveuo(jexnum('&&'//pgc//'.SST.LIA', i), 'E', ltsst)
        do j = 1, nblia*2
            if (zi(ltlia+j-1) .eq. nusst) then
                icomp = icomp+1
                zi(ltsst+icomp-1) = j
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
!--------------------DETERMINATION DE L'ORDRE D'ASSEMBLAGE--------------
!                            DES NOEUDS TARDIFS
!
    ntail = nb_sstr+2*nblia
    call wkvect('&&'//pgc//'.ORD.ASS', 'V V I', ntail, ltoras)
!
!   BOUCLE SUR LES SOUS-STRUCTURES
!
    icomp = 0
    do i = 1, nb_sstr
        call jeveuo(jexnum('&&'//pgc//'.SST.LIA', i), 'L', ltsst)
!
!  BOUCLE SUR LES LIAISONS POUR ASSEMBLAGES DES DUALISATION AVANT
!
        do j = 1, nblia*2
            assok = .true.
            nutarl = zi(ltsst+j-1)
            nulia = nume_equa_gene_orig_l(nutarl)
            if (nutarl .gt. 0) then
!
!   BOUCLE SUR LES NOEUDS TARDIFS DE LIAISON DE LA SOUS-STRUCTURE
!   COURANTE POUR EVITER LES DOUBLES ASSEMBLAGES
!   (NE PAS ASSEMBLER AVANT CE QUI DOIT L'ETRE APRES)
!
                if (j .ne. 1) then
                    do k = 1, j-1
                        nult = zi(ltsst+k-1)
                        null2 = nume_equa_gene_orig_l(nult)
                        if (null2 .eq. nulia .and. nult .ne. 0) assok = .false.
                    end do
                end if
                if (assok) then
                    icomp = icomp+1
                    zi(ltoras+icomp-1) = -nutarl
                end if
            end if
        end do
!
!  ASSEMBLAGE DE LA SOUS-STRUCTURE COURANTE
!
        icomp = icomp+1
        zi(ltoras+icomp-1) = i
!
!   ASSEMBLAGE DES DUALISATIONS APRES LA SOUS-STRUCTURE COURANTE
!
        do j = 1, nblia*2
            assok = .true.
            nutarl = zi(ltsst+j-1)
            if (nutarl .gt. 0) then
                do k = 1, icomp
                    nut = -zi(ltoras+k-1)
                    if (nut .eq. nutarl) assok = .false.
                end do
                if (assok) then
                    icomp = icomp+1
                    zi(ltoras+icomp-1) = -nutarl
                end if
            end if
        end do
    end do
!
    call jedetr('&&'//pgc//'.SST.LIA')
!
!--------------------REMPLISSAGE DES NUMERO D'EQUATION-----------------
!
    icomp = 1
!
    do i = 1, ntail
        nuas = zi(ltoras+i-1)
!
!  CAS DE LA SOUS-STRUCTURE
        if (nuas .gt. 0) then
            nbddl = nume_equa_gene_prno_s((nuas-1)*2+2)
            nume_equa_gene_prno_s(1+(nuas-1)*2) = icomp
!
! CAS DE LA LIAISON
        else
            nbddl = nume_equa_gene_prno_l(-(nuas+1)*2+2)
            nume_equa_gene_prno_l(1-(nuas+1)*2) = icomp
        end if
!
        do j = icomp, icomp+nbddl-1
            zi(ldnueq+j-1) = j
            zi(lddeeq+2*j-1) = nuas
            zi(lddeeq+2*j-2) = j-icomp+1
        end do
        icomp = icomp+nbddl
    end do
!
!----------------------SAUVEGARDES DIVERSES-----------------------------
!
    call jedetr('&&'//pgc//'.ORD.ASS')
!
    call jedema()
end subroutine
