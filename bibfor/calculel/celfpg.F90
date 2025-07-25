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

subroutine celfpg(celz, nomobj, iret)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: celz, nomobj
    integer(kind=8) :: iret
! person_in_charge: jacques.pellet at edf.fr
! ------------------------------------------------------------------
! BUT : EXTRAIRE DU CHAM_ELEM (ELGA) CELZ UN OBJET JEVEUX CONTENANT
!       LE SCHEMA DE POINT DE GAUSS DES MAILLES
! ------------------------------------------------------------------
!     ARGUMENTS:
! CELZ    IN/JXIN  K19 : SD CHAM_ELEM A EXAMINER
!         REMARQUE:   SI CELZ N'EST PAS "ELGA", NOMOBJ N'EST PAS CREE
!
! NOMOBJ  IN/JXVAR K24 : OBJET JEVEUX A CREER (SUR LA BASE 'V')
!    EN SORTIE, L'OBJET NOMOBJ EST UN VECTEUR DE K16 DIMENSIONNE AU
!    NOMBRE DE MAILLES DU MAILLAGE.
!    V(IMA)(1: 8) : NOM DE L'ELREFE POUR LA MAILLE IMA (OU ' ')
!    V(IMA)(9:16) : NOM DE LA FAMILLE DE PG POUR LA MAILLE IMA (OU ' ')
!
!    SI LA FAMILLE DE PG EST UNE FAMILLE "LISTE" (PAR EXEMPLE MATER),
!    V(IMA) CONTIENT LE NUMERO (<0) DE LA FAMILLE. EX : "-185"
!
!    REMARQUE :
!      SI L'OBJET NOMOBJ EXISTE AVANT L'APPEL A CELFPG :
!      * ON VERFIE QUE LE LIGREL EST LE MEME QUE CELUI QUI A SERVI A
!      CREER NOMOBJ
!      * ON NE LE RECALCULE PAS MAIS ON VA VERIFIER
!      QUE LA FAMILLE DES ELEMENTS DE CELZ EST COHERENTE AVEC CELLE
!      DECRITE DANS NOMOBJ. (CAS D'UNE BOUCLE SUR LES NOMSYM)
!      * ON PEUT NE FAIRE LA VERIF QUE SUR LE 1ER ELEMENT DE CHAQUE GREL
!      PUISQUE LE LIGREL EST LE MEME.
!
!
! IRET    CODE RETOUR
!         0 : PAS D'ERREUR
!         1 : CELZ NE CORRESPOND PAS A NOMOBJ (S'IL EXISTE)
    aster_logical :: lexi
!     ------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: ma, nomgd
    character(len=16) :: nofpg
    character(len=19) :: cel, ligrel, ligrsv
    character(len=24) :: nomosv
    integer(kind=8) :: jobj, nbma, jcelv, nec
    integer(kind=8) :: igr, iel, illiel
    integer(kind=8) :: nbgr, imolo, jmolo, numa, nbel, kfpg
    integer(kind=8) :: iexi
    character(len=24), pointer :: celk(:) => null()
    integer(kind=8), pointer :: liel(:) => null()
    integer(kind=8), pointer :: celd(:) => null()
    save ligrsv, nomosv
!
#define numail(igr,iel) liel(zi(illiel+igr-1)+iel-1)
!     ------------------------------------------------------------------
!
    call jemarq()
    cel = celz
    iret = 0
!
!     -- SI CE N'EST PAS UN CHAMP ELGA, IL N'Y A RIEN A FAIRE :
    call jeveuo(cel//'.CELK', 'L', vk24=celk)
    if (celk(3) .ne. 'ELGA') goto 30
!
!
!     1 CALCUL DE LIGREL,NBMA,NEC :
!     --------------------------------------------------------
    call dismoi('NOM_MAILLA', cel, 'CHAM_ELEM', repk=ma)
    call dismoi('NOM_LIGREL', cel, 'CHAM_ELEM', repk=ligrel)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('NOM_GD', cel, 'CHAM_ELEM', repk=nomgd)
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
!
!
!     2 RECUPERATION DES OBJETS DU CHAM_ELEM ET DU LIGREL :
!     -------------------------------------------------------
    call jeveuo(cel//'.CELV', 'L', jcelv)
    call jeveuo(cel//'.CELD', 'L', vi=celd)
    call jeveuo(ligrel//'.LIEL', 'L', vi=liel)
    call jeveuo(jexatr(ligrel//'.LIEL', 'LONCUM'), 'L', illiel)
    nbgr = celd(2)
!
!
!     3 REPLISSAGE DE NOMOBJ :
!     -------------------------------------------------------
!
    call jeexin(nomobj, iexi)
    lexi = (iexi .gt. 0)
!
    if (.not. lexi) then
        call wkvect(nomobj, 'V V K16', nbma, jobj)
        ligrsv = ligrel
        nomosv = nomobj
    else
        ASSERT(nomobj .eq. nomosv)
        ASSERT(ligrel .eq. ligrsv)
        call jeveuo(nomobj, 'L', jobj)
    end if
!
    do igr = 1, nbgr
        nbel = nbelem(ligrel, igr)
        imolo = celd(celd(4+igr)+2)
        if (imolo .eq. 0) goto 20
        call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', jmolo)
        kfpg = zi(jmolo-1+4+nec+1)
        if (kfpg .gt. 0) then
            call jenuno(jexnum('&CATA.TM.NOFPG', kfpg), nofpg)
        else
            call codent(kfpg, 'G', nofpg)
        end if
!
        if (.not. lexi) then
            do iel = 1, nbel
                numa = numail(igr, iel)
                if (numa .gt. 0) zk16(jobj-1+numa) = nofpg
            end do
!
        else
            iel = 1
            numa = numail(igr, iel)
            if (zk16(jobj-1+numa) .ne. nofpg) then
                iret = 1
                goto 30
            end if
        end if
20      continue
    end do
!
!
30  continue
    call jedema()
end subroutine
