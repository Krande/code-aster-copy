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
subroutine exnode(noma, motfac, nzocu, nnocu, nolino)
!
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/reliem.h"
#include "asterfort/wkvect.h"
    character(len=8) :: noma
    character(len=16) :: motfac
    integer(kind=8) :: nzocu
    character(len=24) :: nolino
    integer(kind=8) :: nnocu
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (LIAISON_UNILATERALE - LECTURE)
!
! AFFECTATION DES NOEUDS AFFECTES PAR ZONE
!
! ----------------------------------------------------------------------
!
!
! IN  NOMA   : NOM DU MAILLAGE
! IN  MOTFAC : MOT_CLEF FACTEUR POUR LIAISON UNILATERALE
! IN  NZOCU  : NOMBRE DE ZONES DE LIAISON_UNILATERALE
! OUT NOLINO : NOM DE L'OBJET JEVEUX CONTENANT LA LISTE DES NOEUDS
! IN  NNOCU  : NOMBRE DE TOTAL DE NOEUDS POUR TOUTES LES OCCURRENCES
!
!
!
!
    character(len=8) :: k8bla
    integer(kind=8) :: izone, ino, jdecal, jlist
    integer(kind=8) :: jnl
    integer(kind=8) :: nbmocl
    character(len=16) :: limocl(2), tymocl(2)
    character(len=24) :: listmn, listnn
    integer(kind=8) :: nbmano, nbnono
!
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    jdecal = 1
    nbmocl = 2
    k8bla = ' '
!
! --- NOM DES SD TEMPORAIRES
!
    listmn = '&&EXNODE.MAIL.NOEU'
    listnn = '&&EXNODE.NOEU.NOEU'
!
! --- CREATION DU VECTEUR
!
    call wkvect(nolino, 'V V I', nnocu, jnl)
!
! --- ON LIT LES NOEUDS DES ZONES
!
    do izone = 1, nzocu
        tymocl(1) = 'GROUP_MA'
        tymocl(2) = 'MAILLE'
        limocl(1) = 'GROUP_MA'
        limocl(2) = 'MAILLE'
        call reliem(k8bla, noma, 'NU_NOEUD', motfac, izone, &
                    nbmocl, limocl, tymocl, listmn, nbmano)
        if (nbmano .ne. 0) then
            call jeveuo(listmn, 'L', jlist)
            do ino = 1, nbmano
                zi(jnl+jdecal-1) = zi(jlist+ino-1)
                jdecal = jdecal+1
            end do
        end if
        call jedetr(listmn)
!
        tymocl(1) = 'GROUP_NO'
        tymocl(2) = 'NOEUD'
        limocl(1) = 'GROUP_NO'
        limocl(2) = 'NOEUD'
        call reliem(k8bla, noma, 'NU_NOEUD', motfac, izone, &
                    nbmocl, limocl, tymocl, listnn, nbnono)
!
        if (nbnono .ne. 0) then
            call jeveuo(listnn, 'L', jlist)
            do ino = 1, nbnono
                zi(jnl+jdecal-1) = zi(jlist+ino-1)
                jdecal = jdecal+1
            end do
        end if
        call jedetr(listnn)
!
    end do
!
    call jedema()
!
end subroutine
