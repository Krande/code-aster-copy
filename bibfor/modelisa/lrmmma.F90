! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
! person_in_charge: nicolas.sellenet at edf.fr
!
subroutine lrmmma(fid, nomamd, nbmail, nbnoma, nbtyp,&
                  typgeo, nomtyp, nnotyp, renumd, nmatyp,&
                  nommai, connex, typmai, prefix, infmed,&
                  modnum, numnoa)
!
implicit none
!
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "asterfort/as_mmhcyr.h"
#include "asterfort/as_mmhear.h"
#include "asterfort/codent.h"
#include "asterfort/codlet.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
med_idt :: fid
integer :: nbmail, nbnoma
integer :: nbtyp
integer :: nmatyp(MT_NTYMAX), numnoa(MT_NTYMAX, *), modnum(MT_NTYMAX)
integer :: nnotyp(MT_NTYMAX), typgeo(MT_NTYMAX)
integer :: renumd(*)
integer :: infmed
character(len=6) :: prefix
character(len=8) :: nomtyp(*)
character(len=24) :: connex, nommai, typmai
character(len=*) :: nomamd
!
! --------------------------------------------------------------------------------------------------
!
!     LECTURE DU MAILLAGE -  FORMAT MED - LES MAILLES
!
! --------------------------------------------------------------------------------------------------
!
!     ENTREES :
!       FID    : IDENTIFIANT DU FICHIER MED
!       NOMAMD : NOM DU MAILLAGE MED
!       NBMAIL : NOMBRE DE MAILLES DU MAILLAGE
!       NBNOMA : NOMBRE CUMULE DE NOEUDS PAR MAILLE
!       NBTYP  : NOMBRE DE TYPES POSSIBLES POUR MED
!       TYPGEO : TYPE MED POUR CHAQUE MAILLE
!       NNOTYP : NOMBRE DE NOEUDS POUR CHAQUE TYPE DE MAILLES
!       NOMTYP : NOM DES TYPES POUR CHAQUE MAILLE
!       RENUMD : RENUMEROTATION DES TYPES ENTRE MED ET ASTER
!       NMATYP : NOMBRE DE MAILLES PAR TYPE
!       PREFIX : PREFIXE POUR LES TABLEAUX DES RENUMEROTATIONS
!                A UTILISER PLUS TARD
!       MODNUM : INDICATEUR SI LA SPECIFICATION DE NUMEROTATION DES
!                NOEUDS DES MAILLES EST DIFFERENTES ENTRE ASTER ET MED:
!                     MODNUM = 0 : NUMEROTATION IDENTIQUE
!                     MODNUM = 1 : NUMEROTATION DIFFERENTE
!       NUMNOA : TABLEAU DE CORRESPONDANCE DES NOEUDS MED/ASTER.
!                NUMNOA(ITYP,J) : NUMERO DANS MED DU J IEME NOEUD
!                DE LA MAILLE DE TYPE ITYP DE ASTER
!
!     SORTIES:
!       NOMMAI : NOM DES MAILLES
!       CONNEX : CONNECTIVITE
!       TYPMAI : TYPE DES MAILLES
!
! --------------------------------------------------------------------------------------------------
!
    character(len=6), parameter :: nompro = 'LRMMMA'
    integer, parameter :: edmail=0,edfuin=0,ednoda=0
    integer :: codret
    integer :: ima, idec
    integer :: iaux, jaux
    integer :: imatyp, ityp, letype
    integer :: jcnxma, jmatyp
    integer :: jnomma, jtypma
    integer :: code
    integer :: jnomty(MT_NTYMAX), jnumty(MT_NTYMAX), jcxtyp(MT_NTYMAX)
    integer :: nromai
    integer :: ifm, niv
    character(len=15) :: saux15
    character(len=8) :: saux08
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
!====
! 1. ALLOCATIONS
!====
! 1.1. ==> BASE GLOBALE: OBJETS NOMS / CONNECTIVITE / TYPE DES MAILLES
!
    call wkvect(typmai, 'G V I', nbmail, jtypma)
!
! 1.2. ==> CREATION DE PLUSIEURS VECTEURS PAR TYPE DE MAILLE PRESENT :
!              UN VECTEUR CONTENANT LES NOMS DES MAILLES/TYPE
!           +  UN VECTEUR CONTENANT LES NUMEROS DES MAILLES/TYPE
!           +  UN VECTEUR CONTENANT LA CONNECTIVITE DES MAILLE/TYPE
!              (CONNECTIVITE = NOEUDS + UNE VALEUR BIDON(0) SI BESOIN)
!
    call wkvect('&&'//nompro//'.NOMMAI', 'V V K8', nbmail, jnomma)
    call wkvect('&&'//nompro//'.IMATYP', 'V V I', nbmail*2, jmatyp)
!
    do ityp = 1 , MT_NTYMAX
        if (nmatyp(ityp) .ne. 0) then
            call wkvect('&&'//nompro//'.NOM.'//nomtyp(ityp), 'V V K16', nmatyp(ityp), jnomty(ityp))
            call wkvect('&&'//prefix//'.NUM.'//nomtyp(ityp), 'V V I', nmatyp(ityp), jnumty(ityp))
            iaux = nmatyp(ityp) * nnotyp(ityp)
            call wkvect('&&'//nompro//'.CNX.'//nomtyp(ityp), 'V V I', iaux, jcxtyp(ityp))
        endif
    end do
!
!====
! 2. LECTURE
!    ON PARCOURT TOUS LES TYPES POSSIBLES POUR MED ET ON DECLENCHE LES
!    LECTURES SI DES MAILLES DE CE TYPE SONT PRESENTES DANS LE MAILLAGE
!    REMARQUE : GRACE A LA RENUMEROTATION, ON PARCOURT LES TYPES DE
!    MAILLES DANS L'ORDRE CROISSANT DE LEUR TYPE MED. CE N'EST PAS
!    OBLIGATOIRE SI ON DISPOSE DES TABLEAUX DE NUMEROTATION DES MAILLES.
!    MAIS QUAND CES TABLEAUX SONT ABSENTS, LES CONVENTIONS MED PRECISENT
!    QUE LA NUMEROTATION IMPLICITE SE FAIT DANS CET ORDRE. DONC ON
!    LE FAIT !
!====
!
    nromai = 1
    do letype = 1 , nbtyp
!
! 2.0. ==> PASSAGE DU NUMERO DE TYPE MED AU NUMERO DE TYPE ASTER
!
        ityp = renumd(letype)
!
        if (infmed .ge. 3) then
            write (ifm,201) nomtyp(ityp), nmatyp(ityp)
        endif
201     format('TYPE ',a8,' : ',i10,' MAILLES')
!
        if (nmatyp(ityp) .ne. 0) then
!
! 2.1. ==> LE NUMERO DES MAILLES
!          ON NE TIENT PAS COMPTE DE LA NUMEROTATION DES MAILLES
!          PRESENTE DANS LE FICHIER MED: ON NUMEROTE LES MAILLES
!          SELON LEUR ORDRE D'APPARITION.
!
            do jaux = 1 , nmatyp(ityp)
                zi(jnumty(ityp)+jaux-1) = nromai
                nromai = nromai + 1
            end do
!
!
! 2.2. ==> LE NOM DES MAILLES
!          SI LE FICHIER NE CONTIENT PAS DE NOMMAGE DES MAILLES, ON LEUR
!          DONNE UN NOM PAR DEFAUT FORME AVEC LE PREFIXE 'M' SUIVI DE
!          LEUR NUMERO
!
            call as_mmhear(fid, nomamd, zk16(jnomty(ityp)), edmail, typgeo(ityp), codret)
!
            if (codret .ne. 0) then
                if (infmed .ge. 3) then
                    call utmess('I', 'MED_19', sk=nomtyp(ityp))
                endif
                if (nbmail .ge. 10000000) then
!           + DE 10 MILLIONS DE MAILLES (AU TOTAL), ON PASSE EN BASE 36
                    do iaux = 1, nmatyp(ityp)
                        code = zi(jnumty(ityp)+iaux-1)
                        call codlet(code, 'G', saux15)
                        zk16(jnomty(ityp)+iaux-1) = 'M'//saux15
                    end do
                else
!           MOINS DE 10 MILLIONS DE MAILLES, ON RESTE EN BASE 10
                    do iaux = 1, nmatyp(ityp)
                        code = zi(jnumty(ityp)+iaux-1)
                        call codent(code, 'G', saux15)
                        zk16(jnomty(ityp)+iaux-1) = 'M'//saux15
                    end do
                endif
                codret = 0
            endif
!
! 2.3. ==> LES CONNECTIVITES
!          LA CONNECTIVITE EST FOURNIE EN STOCKANT TOUS LES NOEUDS A
!          LA SUITE POUR UNE MAILLE DONNEE.
!          C'EST CE QUE MED APPELLE LE MODE ENTRELACE
!
            call as_mmhcyr(fid, nomamd, zi(jcxtyp(ityp)), nmatyp(ityp) * nnotyp(ityp), edfuin,&
                           edmail, typgeo(ityp), ednoda, codret)
            if (codret .ne. 0) then
                saux08='mmhcyr'
                call utmess('F', 'DVP_97', sk=saux08, si=codret)
            endif
            do imatyp = 1 , nmatyp(ityp)
                ima = zi(jnumty(ityp)+imatyp-1)
                zk8(jnomma+ima-1) = zk16(jnomty(ityp)+imatyp-1)(1:8)
                zi (jtypma+ima-1) = ityp
                idec = (ima-1)*2
                zi(jmatyp+idec) = ityp
                zi(jmatyp+idec+1) = imatyp
            end do
        endif
    end do
!
!====
! 3. CREATION OBJET CONNECTIVITE DES MAILLES DANS LE BON ORDRE
!====
!
    call jecrec(connex, 'G V I', 'NU', 'CONTIG', 'VARIABLE', nbmail)
    call jeecra(connex, 'LONT', nbnoma)
    call jeecra(connex, 'NUTIOC', nbmail)
!
    do ima = 1 , nbmail
        idec = (ima-1)*2
        ityp = zi(jmatyp+idec)
        call jeecra(jexnum(connex, ima), 'LONMAX', nnotyp(ityp))
        call jeecra(jexnum(connex, ima), 'LONUTI', nnotyp(ityp))
        call jeveuo(jexnum(connex, ima), 'E', jcnxma)
        imatyp = zi(jmatyp+idec+1)
        idec = (imatyp-1)*nnotyp(ityp)
!       POUR LES TYPES DE MAILLE DONT LA NUMEROTATION DES NOEUDS
!       ENTRE ASTER ET MED EST IDENTIQUE:
        if (modnum(ityp) .eq. 0) then
            do jaux = 1 , nnotyp(ityp)
                zi(jcnxma+jaux-1) = zi(jcxtyp(ityp)+jaux-1+idec)
            end do
!       SINON (CF LRMTYP POUR CONNAITRE LA CORRESPONDANCE DES
!       NOEUDS LOCAUX ENTRE ASTER ET MED)
        else
            do jaux = 1 , nnotyp(ityp)
                zi(jcnxma+jaux-1) = zi(jcxtyp(ityp)+numnoa(ityp,jaux)- 1+idec)
            end do
        endif
    end do
!
!====
! 4. CREATION OBJET NOMS DES MAILLES DANS LE BON ORDRE
!====
!
    call jecreo(nommai, 'G N K8')
    call jeecra(nommai, 'NOMMAX', nbmail)
!
    do iaux = 1 , nbmail
        call jecroc(jexnom(nommai, zk8(jnomma+iaux-1)))
    end do
!
!====
! 5. LA FIN
!====
!
    call jedetr('&&'//nompro//'.NOMMAI')
    call jedetr('&&'//nompro//'.IMATYP')
    do ityp = 1 , MT_NTYMAX
        if (nmatyp(ityp) .ne. 0) then
            call jedetr('&&'//nompro//'.NOM.'//nomtyp(ityp))
            call jedetr('&&'//nompro//'.CNX.'//nomtyp(ityp))
        endif
    end do
!
    call jedema()
!
end subroutine
