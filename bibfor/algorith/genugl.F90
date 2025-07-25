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

subroutine genugl(nume_equa, indirf, modgen, mailsk)
    implicit none
!
!***********************************************************************
!  P. RICHARD     DATE 16/11/92
!-----------------------------------------------------------------------
!  BUT: < GENERALE NUMEROTATION GLOBALE >
!
!  CREER, A PARTIR D'UN MODELE GENERALISE ET D'UN MAILLAGE GLOBAL
!  SQUELETTE, LE PROFIL CHAMNO ET UNE FAMILLE NUMEROTEE DONT CHAQUE
!  OBJET CORRESPOND A UNE SOUS-STRUCTURE, ET EST DIMENSIONNE A 2*NBDDL
!  ENGENDRE PAR LES SOUS-STRUCTURES :
!            2*(I-1)+1 --> NUMERO EQUATION DANS NUMDDL SST
!            2*(I-1)+2 --> NUMERO EQUATION DANS PROFNO GLOBAL
!
!-----------------------------------------------------------------------
!
! NUME_EQUA   /I/: NOM K19 DU NUME_EQUA A CREER
! INDIRF  /I/ : NOM K24 DE LA FAMILLE DES INDIRECTIONS A CREER
! MODGEN  /I/ : NOM DU MODELE GENERALISE EN AMONT
! MAILSK  /I/ : NOM DU MAILLAGE SKELETTE
!
!
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/isdeco.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mgutdm.h"
#include "asterfort/utmess.h"
#include "asterfort/profchno_crsd.h"
#include "asterfort/wkvect.h"
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, icomp, iddn, idds, iec, ieq, i_ligr_mesh, ibid, i_ligr_link
    integer(kind=8) :: ipoint, j, k, lddeeq, ldinse, ldnueq
    integer(kind=8) :: ldprno, linueq, llprno, lttds, nbcmp
    integer(kind=8) :: nbcou, nbcpmx, nbddl, nbec, nbnot, nbsst, nddlt
    integer(kind=8) :: ntail, numno, nusst, jrefn
!-----------------------------------------------------------------------
    parameter(nbcpmx=300)
    character(len=6) :: pgc
    character(len=8) :: mailsk, modgen, kbid
    character(len=8) :: k8bid
    character(len=19) :: numddl, nume_equa
    character(len=24) :: indirf, lili, prno, deeq, nueq
    integer(kind=8) :: idec(nbcpmx)
    integer(kind=8), pointer :: vnueq(:) => null()
    integer(kind=8), pointer :: skeleton(:) => null()
!
!-----------------------------------------------------------------------
    data pgc/'GENUGL'/
!-----------------------------------------------------------------------
!
    call jemarq()
    kbid = '  '
    call mgutdm(modgen, kbid, 1, 'NB_CMP_MAX', nbcmp, &
                k8bid)
!
!-----RECUPERATION DU NOMBRE D'ENTIERS CODES DE LA GRANDEUR DEPL_R------
!
    call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nbec)
    if (nbec .gt. 10) then
        call utmess('F', 'MODELISA_94')
    end if
!
!-----RECUPERATION DU NOMBRE DE SOUS-STRUCTURES-------------------------
!
    call jelira(modgen//'      .MODG.SSNO', 'NOMMAX', nbsst)
!
!-----RECUPERATION DIMENSION MAILLAGE SQUELETTE-------------------------
!
    call dismoi('NB_NO_MAILLA', mailsk, 'MAILLAGE', repi=nbnot)
!
!-----RECUPERATION DU .INV.SKELETON-------------------------------------
!
    call jeveuo(mailsk//'.INV.SKELETON', 'L', vi=skeleton)
!
!-----ALLOCATION DU VECTEUR DE TRAVAIL POUR STOCKAGE NOMBRE DE DDL
!     GLOBAUX ENGENDRE PAR SOUS STRUCTURE
!
    call wkvect('&&'//pgc//'.TAIL.DDL.SST', 'V V I', nbsst, lttds)
!
!-----BOUCLE DE COMPTAGE DES DDL FINAUX---------------------------------
!
    nddlt = 0
    do i = 1, nbsst
        kbid = '  '
        call mgutdm(modgen, kbid, i, 'NOM_NUME_DDL', ibid, &
                    numddl)
        numddl(15:19) = '.NUME'
        call jenonu(jexnom(numddl//'.LILI', '&MAILLA'), ibid)
        call jeveuo(jexnum(numddl//'.PRNO', ibid), 'L', llprno)
        do j = 1, nbnot
            nusst = skeleton(j)
            numno = skeleton(1+nbnot+j-1)
            if (nusst .eq. i) then
                nddlt = nddlt+zi(llprno+(numno-1)*(2+nbec)+1)
                zi(lttds+i-1) = zi(lttds+i-1)+zi(llprno+(numno-1)*(2+ &
                                                                   nbec)+1)
            end if
        end do
    end do
!
!-----ALLOCATION DES DIVERS OBJETS-------------------------------------
!
    lili = nume_equa//'.LILI'
    prno = nume_equa//'.PRNO'
    deeq = nume_equa//'.DEEQ'
    nueq = nume_equa//'.NUEQ'
!
! - Create NUME_EQUA
!
    call profchno_crsd(nume_equa, 'G', nb_equa=nddlt, nb_ligrz=2, &
                       prno_lengthz=nbnot*(2+nbec))
    call jeveuo(deeq, 'E', lddeeq)
    call jeveuo(nueq, 'E', ldnueq)
    call jeveuo(nume_equa//'.REFN', 'E', jrefn)
    zk24(jrefn) = mailsk
    zk24(jrefn+1) = 'DEPL_R'
    call wkvect(nume_equa//'.NEQU', 'G V I', 2, jrefn)
    zi(jrefn) = nddlt
    zi(jrefn+1) = nddlt
    call wkvect(nume_equa//'.DELG', 'G V I', nddlt, jrefn)

!
! - Create object LIAISON
!
    call jecroc(jexnom(lili, 'LIAISONS'))
    call jenonu(jexnom(lili, 'LIAISONS'), i_ligr_link)
    call jeecra(jexnum(prno, i_ligr_link), 'LONMAX', 1)

    call jecrec(indirf, 'V V I', 'NU', 'DISPERSE', 'VARIABLE', &
                nbsst)
!
    do i = 1, nbsst
        ntail = 2*zi(lttds+i-1)
        if (ntail .gt. 0) then
            call jeecra(jexnum(indirf, i), 'LONMAX', ntail)
            call jecroc(jexnum(indirf, i))
        end if
    end do
!
!-----REMPLISSAGE DES OBJETS--------------------------------------------
!
    call jenonu(jexnom(lili, '&MAILLA'), i_ligr_mesh)
    call jeveuo(jexnum(prno, i_ligr_mesh), 'E', ldprno)
!
    icomp = 0
!
!  BOUCLE SUR LES SST DU MODELE GENERALISE
!
    do i = 1, nbsst
        nbcou = zi(lttds+i-1)
        idds = 0
!
!  TEST SI LA SST COURANTE ENGENDRE DES DDL GLOBAUX
!
        if (nbcou .gt. 0) then
            kbid = '  '
            call mgutdm(modgen, kbid, i, 'NOM_NUME_DDL', ibid, &
                        numddl)
            numddl(15:19) = '.NUME'
            call jenonu(jexnom(numddl//'.LILI', '&MAILLA'), ibid)
            call jeveuo(jexnum(numddl//'.PRNO', ibid), 'L', llprno)
            call jeveuo(numddl//'.NUEQ', 'L', vi=vnueq)
            call jeveuo(jexnum(indirf, i), 'E', ldinse)
!
!  BOUCLE SUR LES DDL GLOBAUX
!
            do j = 1, nbnot
                nusst = skeleton(j)
!
!  TEST SI LE NOEUD GLOBAL EST ENGENDRE PAR LA SST
!
                if (nusst .eq. i) then
                    numno = skeleton(1+nbnot+j-1)
                    ieq = zi(llprno+(numno-1)*(2+nbec))
                    nbddl = zi(llprno+(numno-1)*(2+nbec)+1)
                    call isdeco(zi(llprno+(numno-1)*(2+nbec)+2), idec, nbcmp)
                    zi(ldprno+(j-1)*(2+nbec)) = icomp+1
                    zi(ldprno+(j-1)*(2+nbec)+1) = nbddl
                    do iec = 1, nbec
                        zi(ldprno+(j-1)*(2+nbec)+1+iec) = zi(llprno+( &
                                                             numno-1)*(2+nbec)+1+iec)
                    end do
                    iddn = 0
                    do k = 1, nbcmp
                        if (idec(k) .gt. 0) then
                            iddn = iddn+1
                            icomp = icomp+1
                            zi(lddeeq+(icomp-1)*2) = j
                            zi(lddeeq+(icomp-1)*2+1) = k
                            zi(ldnueq+icomp-1) = icomp
                            linueq = vnueq(1+ieq+iddn-2)
                            ipoint = ldinse+2*idds-1
                            zi(ipoint+1) = linueq
                            zi(ipoint+2) = icomp
                            idds = idds+1
                        end if
                    end do
                end if
            end do
            call jelibe(numddl//'.NUEQ')
            call jelibe(jexnum(indirf, i))
        end if
    end do
!
!-----SAUVEGARDE DES OBJETS---------------------------------------------
!
    call jedetr('&&'//pgc//'.TAIL.DDL.SST')
!
    call jedema()
end subroutine
