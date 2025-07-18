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

subroutine gesdef(nomres, numddl)
!    P. RICHARD     DATE 19/02/91
!-----------------------------------------------------------------------
!  BUT: < GESTION DES DEFORMEES A CALCULER >
    implicit none
!
!     EN FONCTION DES DDL ACTIFS
!  ET DES MASQUE DE DDL AUX NOEUDS
!   FINIR DE REMPLIR LA TABLEAU DE DESCRIPTION DES DEFORMEES
!
! UTILISATION D'UN TABLEAU VOLATIL DE DESCRIPTION DES DDL PORTES PAR
!   LES NOEUDS:
!
!  ATTENTION: RESTRICTION A 12 DDL PAR NOEUD 6 PHYSIQUES ET LES 6
!            LAGRANGES CORRESPONDANTS
!
!    DX=1    (1)    LAG SUR UN DX=-1   (7)        () ASTER
!    DY=2    (2)    LAG SUR UN DY=-2   (7)
!    DZ=3    (3)    LAG SUR UN DZ=-3   (7)
!    DRX=4    (4)    LAG SUR UN DRX=-4   (7)
!    DRY=5    (5)    LAG SUR UN DRY=-5   (7)
!    DRZ=6    (6)    LAG SUR UN DRZ=-6   (7)
!
!   REMARQUE:LES DDL DU IPRNO NE PORTE QUE DE DDL A CODE POSITIF CAR LE
!     IPRNO CORRESPOND AU PRNO DU LIGREL MAILLAGE -> DECODAGE SUR 6
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM UTILISATEUR DES RESULTATS DE L'OPERATEUR
! NUMDDL   /I/: NOM DE LA NUMEROTATION CORRESPONDANT AU PB
! NBDEF   /O/: NOMBRE TOTAL DE DEFORMEES A CALCULER
!
!
!
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/maskau.h"
#include "asterfort/maskcb.h"
#include "asterfort/maskmn.h"
#include "asterfort/recddl.h"
#include "asterfort/wkvect.h"
    character(len=6) :: pgc
    character(len=8) :: nomres
    character(len=19) :: numddl
    character(len=24) :: desdef, deeq, temmat, temidc
    integer(kind=8) :: ikyp(4)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad1, iad2, ik, lldeeq
    integer(kind=8) :: lldes, ltidec, ltmat, nbcmp, nbdef, nbec, nbno
    integer(kind=8) :: nbnot, nbtem, neq, nomax
    integer(kind=8), pointer :: idc_desc(:) => null()
!-----------------------------------------------------------------------
    data pgc/'GESDEF'/
!-----------------------------------------------------------------------
!
!------------RECUPERATION DONNEES GRANDEURS SOUS-JACENTE----------------
!             ET CREATION VECEUR DE TRAVAIL DECODAGE
!
    call jemarq()
    call dismoi('NB_CMP_MAX', nomres, 'INTERF_DYNA', repi=nbcmp)
    call dismoi('NB_EC', nomres, 'INTERF_DYNA', repi=nbec)
    temidc = '&&'//pgc//'.IDEC'
    call wkvect(temidc, 'V V I', nbcmp*nbec*2, ltidec)
!
!-----------REQUETTE ADRESSE DE LA TABLE DESCRIPTION DES DEFORMEES------
!
    desdef = nomres//'.IDC_DEFO'
    call jeveuo(desdef, 'E', lldes)
    call jelira(desdef, 'LONMAX', nbnot)
    nbnot = nbnot/(2+nbec)
!
!---------------COMPTAGE DES NOEUDS DES DIVERS TYPES INTERFACE----------
!
    do i = 1, 4
        ikyp(i) = 0
    end do
!
    do i = 1, nbnot
        ik = zi(lldes+nbnot+i-1)
        ik = -ik
        ikyp(ik) = ikyp(ik)+1
    end do
!
    nomax = max(ikyp(1), ikyp(2))
    nomax = max(nomax, ikyp(3))
    nomax = max(nomax, ikyp(4))
!
!-----------CREATION MATRICE DES ENTIER CODES DDL ASSEMBLES------------
!    COLONNES SEPARES POUR LES DDL PHYSIQUES ET LES LAGRANGES
!
    nomax = 2*nomax*nbec
    temmat = '&&'//pgc//'.MATDDL'
    call wkvect(temmat, 'V V I', nomax, ltmat)
!
!--------------------REQUETE SUR LE DEEQ DU NUMDDL----------------------
!
!
!
    deeq = numddl//'.DEEQ'
    call jeveuo(deeq, 'L', lldeeq)
    call jelira(deeq, 'LONMAX', neq)
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!
!--------------TRAITEMENT DES MODES D'ATTACHE (MAC NEAL)----------------
!
!
!   DECALAGE EVENTUEL DE LA LISTE DES NOEUDS MN DANS LA LISTE GLOBALE
!
    nbtem = 0
    nbdef = 0
!
    nbno = ikyp(1)
!
    call recddl(nbcmp, zi(lldes+nbtem), nbno, nbec, zi(lldeeq), &
                neq, zi(ltmat), zi(ltidec))
!
    iad1 = lldes+nbnot*2+nbtem*nbec
    iad2 = lldes+nbnot+nbtem
    call maskmn(nbcmp, nbno, nbec, zi(ltmat), zi(iad1), &
                zi(iad2), nbdef)
!
!----------TRAITEMENT DES MODES CONTRAINTS (CRAIG BAMPTON)--------------
!
!   DECALAGE EVENTUEL DE LA LISTE DES NOEUDS MN DANS LA LISTE GLOBALE
!
    nbtem = ikyp(1)
!
    nbno = ikyp(2)
!
    call recddl(nbcmp, zi(lldes+nbtem), nbno, nbec, zi(lldeeq), &
                neq, zi(ltmat), zi(ltidec))
!
    iad1 = lldes+nbnot*2+nbtem*nbec
    iad2 = lldes+nbnot+nbtem
    call maskcb(nbcmp, nbno, nbec, zi(ltmat), zi(iad1), &
                zi(iad2), nbdef)
!
!-------TRAITEMENT DES MODES CONTRAINTS HARMONIQUES(CB-HARMO)-----------
!
!   DECALAGE EVENTUEL DE LA LISTE DES NOEUDS MN DANS LA LISTE GLOBALE
!
    nbtem = ikyp(1)+ikyp(2)
!
    nbno = ikyp(3)
!
    call recddl(nbcmp, zi(lldes+nbtem), nbno, nbec, zi(lldeeq), &
                neq, zi(ltmat), zi(ltidec))
!
    iad1 = lldes+nbnot*2+nbtem*nbec
    iad2 = lldes+nbnot+nbtem
    call maskcb(nbcmp, nbno, nbec, zi(ltmat), zi(iad1), &
                zi(iad2), nbdef)
!
!-----------------TRAITEMENT DES NOEUDS D'INTERFACE AUCUN---------------
!
!
!   DECALAGE EVENTUEL DE LA LISTE DES NOEUDS AU DANS LA LISTE GLOBALE
!
    nbtem = ikyp(1)+ikyp(2)+ikyp(3)
!
    nbno = ikyp(4)
!
    call recddl(nbcmp, zi(lldes+nbtem), nbno, nbec, zi(lldeeq), &
                neq, zi(ltmat), zi(ltidec))
!
    iad1 = lldes+nbnot*2+nbtem*nbec
    iad2 = lldes+nbnot+nbtem
    call maskau(nbno, nbec, zi(iad1))
!
!------------------------FINITION DU .DESC------------------------------
!
    call jeveuo(nomres//'.IDC_DESC', 'E', vi=idc_desc)
    idc_desc(5) = nbdef
!
!---------------------LIBERATION DES OBJETS-----------------------------
!
!
!-------------------DESTRUCTION DES OBJETS VOLATILES--------------------
!
    call jedetr(temmat)
    call jedetr(temidc)
!
    call jedema()
end subroutine
