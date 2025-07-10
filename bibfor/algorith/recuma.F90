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
subroutine recuma(mailla, nbma, nbgr, nomma, nomgr, &
                  nbto, numnot)
!    P. RICHARD     DATE 13/07/90
!-----------------------------------------------------------------------
!  BUT: RASSEMBLER LES MAILLES DE NOMMA ET DES GROUPNO DE NOMGR
    implicit none
!          ET TRANSCODER DANS NUMNOT
!
!-----------------------------------------------------------------------
!
! MAILLA   /I/: NOM UTILISATEUR DU MAILLAGE
! NBMA     /I/: NOMBRE DE MAILLE EN ARGUMENT DE LA COMMANDE
! NBGR     /I/: NOMBRE DE GROUPES DE MAILLES EN ARGUMENTS
! NOMMA    /I/: NOMS DES MAILLES DONNES EN ARGUMENTS
! NOMGR    /I/: NOMS DES GROUPES DE MAILLES EN ARGUMENTS
! NBTO     /O/: NOMBRE TOTAL DE MAILLES ASSOCIES A L'INTERFACE
! NUMNOT   /O/: VECTEUR DES NUMERO DES MAILLES D'INTERFACE
!
!
!
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: nbma, nbgr, nbto, numnot(nbto)
    character(len=8) :: mailla, nomma(nbma)
    character(len=24) :: valk(2), nomgr(nbgr), nomcou
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadg, icomp, j, nb, numa
!-----------------------------------------------------------------------
    call jemarq()
    icomp = 0
!
!-------RECUPERATION ET TRANSCODAGE DES MAILLES DES GROUPES-------------
!
    if (nbgr .gt. 0) then
        do i = 1, nbgr
            nomcou = nomgr(i)
            call jelira(jexnom(mailla//'.GROUPEMA', nomcou), 'LONUTI', nb)
            call jeveuo(jexnom(mailla//'.GROUPEMA', nomcou), 'L', iadg)
            do j = 1, nb
                icomp = icomp+1
                numnot(icomp) = zi(iadg+j-1)
            end do
        end do
    end if
!
!
!-------RECUPERATION ET TRANSCODAGE DES MAILLES-------------------------
!
!
!
    if (nbma .gt. 0) then
        do i = 1, nbma
            nomcou = nomma(i)
            numa = char8_to_int(nomcou)
!
            if (numa .eq. 0) then
                valk(1) = mailla
                valk(2) = nomcou
                call utmess('E', 'ALGORITH14_10', nk=2, valk=valk)
            end if
!
            icomp = icomp+1
            numnot(icomp) = numa
!
        end do
    end if
    nbto = icomp
!
    call jedema()
end subroutine
