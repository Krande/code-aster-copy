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
subroutine gilig2(nfic, nbnono, niv)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nfic, nbnono, niv
!
!     BUT: LIRE LES N LIGNES DES POINTS DU MAILLAGE GIBI :
!                 ( PROCEDURE SAUVER)
!
!     IN: NFIC   : UNITE DE LECTURE
!         NBNONO : NOMBRE D'OBJETS NOMMES
!         NIV    : NIVEAU GIBI
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbnom, nbnum, iaptno, iaptnu, nbfois, nbrest, icoj, i, j
!     ------------------------------------------------------------------
!
    call jemarq()
!
!     -- ON LIT LES NOEUDS NOMMES:
!
!
    if (nbnono .gt. 0) then
!
        if (niv .eq. 3) then
            nbnom = 8
            nbnum = 16
        else
            nbnom = 8
            nbnum = 10
        end if
!
!     -- ON CREE LES 2 OBJETS QUI CONTIENDRONT LES NOMS ET NUMEROS
!        DES POINTS NOMMES:
!
        call wkvect('&&GILIRE.POINT_NOM', 'V V K8', nbnono, iaptno)
        call wkvect('&&GILIRE.POINT_NUM', 'V V I', nbnono, iaptnu)
!
        nbfois = nbnono/nbnom
        nbrest = nbnono-nbnom*nbfois
        icoj = 0
        do i = 1, nbfois
            read (nfic, 1007) (zk8(iaptno-1+j), j=icoj+1, icoj+nbnom)
            icoj = icoj+nbnom
        end do
        if (nbrest .gt. 0) then
            read (nfic, 1007) (zk8(iaptno-1+j), j=icoj+1, icoj+nbrest)
        end if
!
        nbfois = nbnono/nbnum
        nbrest = nbnono-nbnum*nbfois
        icoj = 0
        do i = 1, nbfois
            if (niv .eq. 3) then
                read (nfic, 1009) (zi(iaptnu-1+j), j=icoj+1, icoj+nbnum)
                icoj = icoj+nbnum
            else
                read (nfic, 1008) (zi(iaptnu-1+j), j=icoj+1, icoj+nbnum)
                icoj = icoj+nbnum
            end if
        end do
        if (nbrest .gt. 0) then
            if (niv .eq. 3) then
                read (nfic, 1009) (zi(iaptnu-1+j), j=icoj+1, icoj+nbrest)
            else
                read (nfic, 1008) (zi(iaptnu-1+j), j=icoj+1, icoj+nbrest)
            end if
        end if
    end if
!
1007 format(8(1x, a8))
1008 format(10(i8))
1009 format(16(i5))
!
    call jedema()
!
end subroutine
