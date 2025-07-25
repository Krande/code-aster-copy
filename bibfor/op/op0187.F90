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

subroutine op0187()
    implicit none
! person_in_charge: samuel.geniaut at edf.fr
!     =================================================================
!                      OPERATEUR POST_MAIL_XFEM
!                      ------------------------
!     BUT : GENERER UN MAILLAGE DESTINE UNIQUEMENT AU POST-TRAITEMENT
!           DES ELEMENTS XFEM, ET METTANT EN EVIDENCE LES SOUS-ELEMENTS
!           DES MAILLES FISSUREES
!     =================================================================
!     ------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/cargeo.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/titre.h"
#include "asterfort/wkvect.h"
#include "asterfort/xpocrf.h"
#include "asterfort/xpodim.h"
#include "asterfort/xpofon.h"
#include "asterfort/xpoini.h"
#include "asterfort/xpomac.h"
#include "asterfort/xpomax.h"
#include "asterfort/xposep.h"
    integer(kind=8) :: ibid, iret, nsetot, nnntot, ncotot, nbnoc, nbmac, ifm, niv
    integer(kind=8) :: nbgma2, jnivgr, nftot, mftot, nfcomf, ngfon
    character(len=2) :: prefno(4)
    character(len=8) :: maxfem, mo, malini, k8b, nomres
    character(len=16) :: k16b
    character(len=19) :: k19b
    character(len=24) :: mailx, mailc, listno, k24b, logrma, dirgrm, listgr
    character(len=24) :: nivgrm
    character(len=24) :: nogrfi
    aster_logical :: pre1
!
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
    pre1 = .false.
!
!     ------------------------------------------------------------------
!     1. RECUPERATION DES CONCEPTS UTILISATEURS
!        ET DEFINITION DES NOMBRES DE MAILLES ET DE NOEUDS EN FOND DE
!        FISSURE
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '1. XPOINI'
    call xpoini(maxfem, mo, malini, k8b, k24b, &
                k8b, k8b, prefno, nogrfi)
    call xpofon(mo, mftot, nftot, nfcomf, ngfon)
!
!     ------------------------------------------------------------------
!     2. SEPARATION DES MAILLES DE MALINI EN 2 GROUPES
!              - MAILC : MAILLES NON AFFECTEES D'UN MODELE
!                        OU NON SOUS-DECOUPEES (CLASSIQUE)
!              - MAILX : MAILLES SOUS-DECOUPEES (X-FEM)
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '2. XPOSEP'
    mailc = '&&OP0187.MAILC'
    mailx = '&&OP0187.MAILX'
    logrma = '&&OP0187.LOGRMA'
    listgr = '&&OP0187.LISTGR'
    call xposep(mo, malini, mailc, mailx, nsetot, &
                nnntot, ncotot, logrma, listgr)
!
    if (niv .gt. 1) then
        write (ifm, *) 'NOMBRE DE NOUVELLES MAILLES A CREER', nsetot+ &
            mftot
        write (ifm, *) 'NOMBRE DE NOUVEAUX NOEUDS A CREER', nnntot+nftot
    end if
!
!     ------------------------------------------------------------------
!     3. DIMENSIONNEMENT DES OBJETS DE MAXFEM
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '3. XPODIM'
    listno = '&&OP0187.LISTNO'
    dirgrm = '&&OP0187.DIRGRM'
    call xpodim(malini, mailc, k8b, k24b, nsetot+mftot, &
                nnntot+nftot, ncotot+nfcomf, listno, k19b, k19b, &
                k19b, k19b, k19b, k19b, k19b, &
                ibid, k8b, nbnoc, nbmac, logrma, &
                dirgrm, maxfem, ngfon, k19b, k19b, &
                pre1, mo)
!
!     ------------------------------------------------------------------
!     4. TRAITEMENT DES MAILLES DE MAILC
!            LES MAILLES DE MAILC ET LES NOEUDS ASSOCIÉS SONT COPIES
!            DANS MAXFEM A L'IDENTIQUE
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '4. XPOMAC'
!
!     CREATION DU VECTEUR DE REMPLISSAGE DES GROUP_MA
    nivgrm = '&&OP0187.NIVGRM'
    call jelira(maxfem//'.GROUPEMA', 'NUTIOC', nbgma2)
    if (nbgma2 .gt. 0) call wkvect(nivgrm, 'V V I', nbgma2, jnivgr)
!
    call xpomac(malini, mailc, listno, nbnoc, nbmac, &
                maxfem, nivgrm, k19b, k19b, k19b, &
                k19b, k19b, k19b, k8b, k19b, &
                k19b, pre1)
!
!     ------------------------------------------------------------------
!     5. TRAITEMENT DES MAILLES DE MAILX
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '5. XPOMAX'
    call xpomax(mo, malini, mailx, nbnoc, nbmac, &
                prefno, nogrfi, maxfem, k19b, k19b, &
                k19b, k19b, k19b, k19b, listgr, &
                dirgrm, nivgrm, k8b, ngfon, k19b, &
                k19b, pre1, 0)
!
!     ------------------------------------------------------------------
!     6. TRAITEMENT DES FONDS DE FISSURE
!     ------------------------------------------------------------------
!
    if (niv .gt. 1) write (ifm, *) ' '
    if (niv .gt. 1) write (ifm, *) '6. XPOCRF'
    call xpocrf(mo, maxfem, mftot, nftot)
!
    if (niv .gt. 1) write (ifm, *) 'FIN DE POST_MAIL_XFEM'
!
    call titre()
!
! --- CARACTERISTIQUES GEOMETRIQUES :
!     -----------------------------
    call getres(nomres, k16b, k16b)
    call cargeo(nomres)
!
    call jeexin(dirgrm, iret)
    if (iret .ne. 0) call jedetr(dirgrm)
!
    call jeexin(listgr, iret)
    if (iret .ne. 0) call jedetr(listgr)
!
    call jeexin(nivgrm, iret)
    if (iret .ne. 0) call jedetr(nivgrm)
!
    call jedema()
end subroutine
