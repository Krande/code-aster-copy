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
subroutine discax(noma, nbn, iaxe, nuno, diax)
    implicit none
!     CREATION D'UNE LISTE ORDONNEE DE NOEUDS SUR UNE STRUCTURE POUTRE
!     DROITE : ORDRE CROISSANT DU PARAMETRE LE LONG DE L'AXE DIRECTEUR
!     DE LA POUTRE
!     APPELANT : SPECFF
!-----------------------------------------------------------------------
! IN  : NOMA   : NOM DU CONCEPT MAILLAGE
! IN  : NBN    : NOMBRE DE NOEUDS DU MAILLAGE
! IN  : IAXE   : ENTIER DEFINISSANT L'AXE DIRECTEUR
!       IAXE = 1 L'AXE DIRECTEUR EST L'AXE DES X DU REPERE GLOBAL
!       IAXE = 2 L'AXE DIRECTEUR EST L'AXE DES Y DU REPERE GLOBAL
!       IAXE = 3 L'AXE DIRECTEUR EST L'AXE DES Z DU REPERE GLOBAL
! OUT : NUNO   : LISTE DES NUMEROS DES NOEUDS DU MAILLAGE, REORDONNEE
!                PAR VALEURS CROISSANTES DU PARAMETRE LE LONG DE L'AXE
! OUT : DIAX   : LISTE DES VALEURS DU PARAMETRE LE LONG DE L'AXE
!                ORDRE CROISSANT
!
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: noma
    integer(kind=8) :: nbn, iaxe, nuno(nbn)
    real(kind=8) :: diax(nbn)
!
    character(len=8) :: nomnoe
    character(len=24) :: coorma
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: icoma, imin, innoe, ino, jno
    real(kind=8) :: xmin
!-----------------------------------------------------------------------
    call jemarq()
!
! --- 1.ACCES AUX OBJETS DU CONCEPT MAILLAGE
!
    coorma = noma//'.COORDO    .VALE'
    call jeveuo(coorma, 'L', icoma)
!
! --- 2.ON RECOPIE LA DISCRETISATION (NON ORDONNEE) LUE DANS L'OBJET
! ---   .COORDO    .VALE DU CONCEPT MAILLAGE
! ---   ON RECOPIE SIMULTANEMENT LA LISTE DES NOMS DES NOEUDS
!
    call wkvect('&&DISCAX.TEMP.NNOE', 'V V K8', nbn, innoe)
    do ino = 1, nbn
        diax(ino) = zr(icoma+3*(ino-1)+iaxe-1)
        zk8(innoe+ino-1) = int_to_char8(ino)
    end do
!
! --- 3.ON REORDONNE LA DISCRETISATION PAR VALEURS CROISSANTES
! ---   ON REORDONNE SIMULTANEMENT LA LISTE DES NOMS DES NOEUDS
! ---   ON EN DEDUIT LA LISTE ORDONNEE DES NUMEROS DES NOEUDS
!
    do ino = 1, nbn-1
        xmin = diax(ino)
        nomnoe = zk8(innoe+ino-1)
        imin = ino
        do jno = ino+1, nbn
            if (diax(jno) .lt. xmin) then
                xmin = diax(jno)
                nomnoe = zk8(innoe+jno-1)
                imin = jno
            end if
        end do
        diax(imin) = diax(ino)
        zk8(innoe+imin-1) = zk8(innoe+ino-1)
        diax(ino) = xmin
        zk8(innoe+ino-1) = nomnoe
        nuno(ino) = char8_to_int(zk8(innoe+ino-1))
    end do
    nuno(nbn) = char8_to_int(zk8(innoe+nbn-1))
!
    call jedetr('&&DISCAX.TEMP.NNOE')
    call jedema()
end subroutine
