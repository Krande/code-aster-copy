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
subroutine nbnode(noma, motfac, nzocu, nopono, nnocu)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/reliem.h"
#include "asterfort/wkvect.h"
    character(len=8) :: noma
    character(len=16) :: motfac
    integer(kind=8) :: nzocu
    character(len=24) :: nopono
    integer(kind=8) :: nnocu
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (LIAISON_UNILATERALE - LECTURE)
!
! DECOMPTE DES NOEUDS AFFECTES PAR ZONE
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  MOTFAC : MOT_CLEF FACTEUR POUR LIAISON UNILATERALE
! IN  NZOCU  : NOMBRE DE ZONES DE LIAISON_UNILATERALE
! OUT NOPONO : NOM DE L'OBJET JEVEUX CONTENANT LE VECTEUR D'INDIRECTION
! OUT NNOCU  : NOMBRE DE TOTAL DE NOEUDS POUR TOUTES LES OCCURRENCES
!
!
!
!
    character(len=8) :: k8bla
    integer(kind=8) :: izone
    integer(kind=8) :: jnp
    integer(kind=8) :: nbmocl
    character(len=16) :: limocl(2), tymocl(2)
    character(len=24) :: listmn, listnn
    integer(kind=8) :: nbmano, nbnono, nbno
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    call wkvect(nopono, 'V V I', nzocu+1, jnp)
!
    zi(jnp) = 1
    nnocu = 0
    nbmocl = 2
    k8bla = ' '
!
! --- NOM DES SD TEMPORAIRES
!
    listmn = '&&NBNODE.MAIL.NOEU'
    listnn = '&&NBNODE.NOEU.NOEU'
!
! --- ON COMPTE LES NOEUDS DES ZONES
!
    do izone = 1, nzocu
        tymocl(1) = 'GROUP_MA'
        tymocl(2) = 'MAILLE'
        limocl(1) = 'GROUP_MA'
        limocl(2) = 'MAILLE'
        call reliem(k8bla, noma, 'NU_NOEUD', motfac, izone, &
                    nbmocl, limocl, tymocl, listmn, nbmano)
!
        tymocl(1) = 'GROUP_NO'
        tymocl(2) = 'NOEUD'
        limocl(1) = 'GROUP_NO'
        limocl(2) = 'NOEUD'
        call reliem(k8bla, noma, 'NU_NOEUD', motfac, izone, &
                    nbmocl, limocl, tymocl, listnn, nbnono)
        nbno = nbmano+nbnono
        nnocu = nnocu+nbno
        zi(jnp+izone) = zi(jnp+izone-1)+nbno
    end do
!
    call jedetr(listmn)
    call jedetr(listnn)
!
    call jedema()
!
end subroutine
