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

subroutine cynupl(nume_equa, indirf, modcyc, mailsk, nbsec)
    implicit none
!
!***********************************************************************
!    O. NICOLAS     DATE 10/01/06
!-----------------------------------------------------------------------
!  BUT: < CYCLIQUE NUMEROTATION GLOBALE >
!
!    CREER A PARTIR D'UN RESULTAT CYCLIQUE ET UN MAILLAGE GLOBAL
!  SQUELLETTE LE PROFIL CHAMNO ET UNE FAMILLE NUMEROTEE
! ,DONT CHAQUE OBJET CORRESPOND A UN SECTEUR, ET EST DIMENSIONNE
!  A 2*NBDDL ENGENDRE PAR LE SECTEUR:C
!
!            2*(I-1)+1 --> NUMERO EQUATION DANS PFCHNO SECTEUR
!            2*(I-1)+2 --> NUMERO EQUATION DANS PROFNO GLOBAL
!
!-----------------------------------------------------------------------
!
! NOM----- / /:
!
! NUME_EQUA   /I/: NOM K19 DU NUME_EQUA A CREER
! INDIRF   /I/: NOM K24 DE LA FAMILLE DES INDIRECTIONS A CREER
! MODCYC   /I/: NOM DU RESULTAT CYCLIQUE EN AMONT
! MAILSK   /I/: NOM DU MAILLAGE SKELETTE
! NBSEC    /I/: NBRE DE SECTEUR
!
!
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/isdeco.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/profchno_crsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/wkvect.h"
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, icomp, iec, ieq, ier, i_ligr_mesh, i_ligr_link
    integer(kind=8) :: ipoint, j, lddeeq, ldnueq, ldprno, linueq
    integer(kind=8) ::   llprno, ltinse, lttds, nbcmp, nbcpmx
    integer(kind=8) :: nbddl, nbnot, nbsec, nddlt, neqsec, nsecpr
    integer(kind=8) :: ntail, nugd, numnos, numsec
!-----------------------------------------------------------------------
    parameter(nbcpmx=320)
    character(len=6) :: pgc
    character(len=8) :: modcyc, mailsk, nomgd
    character(len=19) :: nume_equa_sec, nume_equa, chamno
    character(len=24) :: indirf, lili, prno, deeq, nueq
    integer(kind=8) :: idec(nbcpmx), nec
    integer(kind=8), pointer :: skeleton(:) => null()
    integer(kind=8), pointer :: vnueq(:) => null()
!
!-----------------------------------------------------------------------
!
    call jemarq()
    pgc = 'CYNUPL'
!
!
!----------------RECUPERATION DU NUME_EQUA:
    call rsexch('F', modcyc, 'DEPL', 1, chamno, &
                ier)
    call dismoi('NUME_EQUA', chamno, 'CHAM_NO', repk=nume_equa_sec)
!
!---------------RECUPERATION DU NOMBRE DE COMPOSANTES-------------------
!
!     -- QUESTION "POURRIE" :
    call dismoi('NOM_GD', nume_equa_sec, 'NUME_EQUA', repk=nomgd)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbcmp)
    call jenonu(jexnom('&CATA.GD.NOMGD', nomgd), nugd)
    nec = nbec(nugd)
    ASSERT(nec .le. 11)
!
!
!-------------RECUPERATION DIMENSION MAILLAGE SQUELETTE-----------------
!
    call dismoi('NB_NO_MAILLA', mailsk, 'MAILLAGE', repi=nbnot)
!
!------------RECUPERATION DU .INV.SKELETON------------------------------
!
    call jeveuo(mailsk//'.INV.SKELETON', 'L', vi=skeleton)
!
!--------------RECUPERATION DU PRNO DU SECTEUR--------------------------
!
    call jenonu(jexnom(nume_equa_sec//'.LILI', '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(nume_equa_sec//'.PRNO', i_ligr_mesh), 'L', llprno)
    call jeveuo(nume_equa_sec//'.NUEQ', 'L', vi=vnueq)
    call dismoi('NB_EQUA', nume_equa_sec, 'NUME_EQUA', repi=neqsec)
!
!--------------------ALLOCATION DU VECTEUR DE TRAVAIL-------------------
!     POUR STOCKAGE NOMBRE DE DDL GLOBAUX ENGENDRE PAR SECTEUR
!
    call wkvect('&&'//pgc//'.TAIL.DDL.SECT', 'V V IS', nbsec, lttds)
!
!--------------BOUCLE DE COMPTAGE DES DDL FINAUX------------------------
!
    nddlt = 0
    do i = 1, nbnot
        numsec = skeleton(i)
        numnos = skeleton(1+nbnot+i-1)
        nddlt = nddlt+zi(llprno+(numnos-1)*(2+nec)+1)
        zi(lttds+numsec-1) = zi(lttds+numsec-1)+zi(llprno+(numnos-1)*( &
                                                   2+nec)+1)
    end do
!
!-----------------ALLOCATION DES DIVERS OBJETS--------------------------
!
    lili = nume_equa//'.LILI'
    prno = nume_equa//'.PRNO'
    deeq = nume_equa//'.DEEQ'
    nueq = nume_equa//'.NUEQ'
!
! - Create NUME_EQUA
!
    call profchno_crsd(nume_equa, 'G', nb_equa=nddlt, nb_ligrz=2, &
                       prno_lengthz=nbnot*(2+nec))
    call jeveuo(deeq, 'E', lddeeq)
    call jeveuo(nueq, 'E', ldnueq)
!
! - Create object LIAISON
!
    call jecroc(jexnom(lili, 'LIAISONS'))
    call jenonu(jexnom(lili, 'LIAISONS'), i_ligr_link)
    call jeecra(jexnum(prno, i_ligr_link), 'LONMAX', 1)
!
!
    call jecrec(indirf, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbsec)
!
!
    do i = 1, nbsec
        call jecroc(jexnum(indirf, i))
        ntail = 2*zi(lttds+i-1)
        call jeecra(jexnum(indirf, i), 'LONMAX', ntail)
        zi(lttds+i-1) = 0
    end do
!
!-------------------------REMPLISSAGE DES OBJETS------------------------
!
    call jenonu(jexnom(lili, '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(prno, i_ligr_mesh), 'E', ldprno)
!
    nsecpr = 1
    call jeveuo(jexnum(indirf, nsecpr), 'E', ltinse)
    icomp = 0
    do i = 1, nbnot
!
        numsec = skeleton(i)
        numnos = skeleton(1+nbnot+i-1)
        ieq = zi(llprno+(numnos-1)*(2+nec))
        nbddl = zi(llprno+(numnos-1)*(2+nec)+1)
        call isdeco(zi(llprno+(numnos-1)*(2+nec)+2), idec, nbcmp)
!
        zi(ldprno+(i-1)*(2+nec)) = icomp+1
        zi(ldprno+(i-1)*(2+nec)+1) = nbddl
        do iec = 1, nec
            zi(ldprno+(i-1)*(2+nec)+1+iec) = zi(llprno+(numnos-1)*(2+ &
                                                                   nec)+1+iec)
        end do
        if (numsec .ne. nsecpr) then
            call jelibe(jexnum(indirf, nsecpr))
            nsecpr = numsec
            call jeveuo(jexnum(indirf, nsecpr), 'E', ltinse)
        end if
        iad = 0
        do j = 1, nbcmp
            if (idec(j) .gt. 0) then
                iad = iad+1
                icomp = icomp+1
                zi(lddeeq+(icomp-1)*2) = i
                zi(lddeeq+(icomp-1)*2+1) = j
                zi(ldnueq+icomp-1) = icomp
                linueq = vnueq(1+ieq+iad-2)
                ipoint = ltinse-1+2*zi(lttds+numsec-1)
                zi(ipoint+1) = linueq
                zi(ipoint+2) = icomp
                zi(lttds+numsec-1) = zi(lttds+numsec-1)+1
            end if
        end do
    end do
!
    call jelibe(jexnum(indirf, nsecpr))
!
!------------------SAUVEGARDE DES OBJETS--------------------------------
!
    call jedetr('&&'//pgc//'.TAIL.DDL.SECT')
!
    call jedema()
end subroutine
