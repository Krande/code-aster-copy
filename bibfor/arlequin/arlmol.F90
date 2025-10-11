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
!
subroutine arlmol(nomo, mailar, modarl, tabcor)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/adalig.h"
#include "asterfort/assert.h"
#include "asterfort/cormgi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/initel.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
!     ARGUMENTS:
!     ----------
    character(len=8) :: mailar, modarl
    character(len=8) :: nomo
    character(len=24) :: tabcor
!
    integer(kind=8) :: ima, nbma, ibid
    integer(kind=8) :: jnbno, jad, jlgrf, jdime, jtyel, jtabco
    character(len=8) :: k8bid
    character(len=19) :: ligarl, ligrel
    integer(kind=8) :: numori, ityel, iret
!
! ----------------------------------------------------------------------
!
! ROUTINE ARLEQUIN
!
! CREATION DU LIGREL DU PSEUDO-MODELE
!
! ----------------------------------------------------------------------
!
!
! IN  NOMO   : NOM DU MODELE
! IN  MAILAR : NOM DU PSEUDO-MAILLAGE
! IN  TABCOR : TABLEAU DE CORRESPONDANCE
!            POUR CHAQUE NOUVEAU NUMERO ABSOLU DANS MAILAR
!             -> ANCIEN NUMERO ABSOLU DANS MAIL
!             -> SI NEGATIF, LA NOUVELLE MAILLE EST ISSUE D'UNE
!                DECOUPE DE LA MAILLE DE NUMERO ABSOLU ABS(NUM) DANS
!                MAIL
! IN  MODARL : NOM DU PSEUDO-MODELE
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ligarl = ' '
    call exisd('MODELE', modarl, iret)
    if (iret .ne. 0) then
        call dismoi('NOM_LIGREL', modarl, 'MODELE', repk=ligarl)
    end if
!
! --- DESTRUCTION DU LIGREL S'IL EXISTE
!
    call exisd('LIGREL', ligarl, iret)
    if (iret .ne. 0) then
        call detrsd('LIGREL', ligarl)
    else
        ligarl = modarl(1:8)//'.MODELE'
    end if
!
! --- INFORMATIONS SUR LE MODELE ORIGINAL
!
    call dismoi('NOM_LIGREL', nomo, 'MODELE', repk=ligrel)
    call jeveuo(ligrel//'.TYFE', 'L', jtyel)
!
! --- ACCES AU TABLEAU DE CORRESPONDANCE
!
    call jeveuo(tabcor, 'L', jtabco)
!
! --- INFORMATIONS SUR LE MAILLAGE
!
    call jeveuo(mailar(1:8)//'.DIME', 'L', jdime)
    nbma = zi(jdime-1+3)
!
! --- CREATION DE .NOMA + ATTRIBUT DOCU
!
    call wkvect(ligarl//'.LGRF', 'V V K8', 4, jlgrf)
    zk8(jlgrf-1+1) = mailar
    call jeecra(ligarl//'.LGRF', 'DOCU', ibid, 'MECA')
!
! --- CREATION DE L'OBJET .LIEL: ON LE CREE AU MAX. AUTANT DE LIEL
! --- QUE DE MAILLES
!
    call jecrec(ligarl//'.LIEL', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbma)
    call jeecra(ligarl//'.LIEL', 'LONT', 2*nbma, k8bid)
!
    do ima = 1, nbma
!
! --- CREATION OBJET DE LA COLLECTION
!
        call jecroc(jexnum(ligarl//'.LIEL', ima))
        call jeecra(jexnum(ligarl//'.LIEL', ima), 'LONMAX', 2, k8bid)
        call jeveuo(jexnum(ligarl//'.LIEL', ima), 'E', jad)
!
! --- NUMERO DANS LE MAILLAGE ORIGINAL
!
        numori = zi(jtabco+ima-1)
!
! --- TYPE DE LA MAILLE DANS LE PSEUDO-MODELE
!
        if (numori < 0) then
            numori = abs(numori)
        end if
        ityel = zi(jtyel-1+numori)
!
! --- PAS D'EF AFFECTE SUR LA MAILLE !
!
        if (ityel == 0) then
            ASSERT(.false.)
        else
            zi(jad) = ima
            zi(jad+1) = ityel
        end if
    end do
!
! --- PAS DE NOEUDS TARDIFS
!
    call wkvect(ligarl//'.NBNO', 'V V I', 1, jnbno)
    zi(jnbno-1+1) = 0
!
! --- ADAPTATION DE .LIEL
!
    call adalig(ligarl)
!
! --- CREATION DE L'OBJET .REPE
!
    call cormgi('V', ligarl)
!
! --- INITIALISATION DU LIGREL (OBJETS PRNM/PRNS)
!
    call initel(ligarl)
!
    call jedema()
!
end subroutine
