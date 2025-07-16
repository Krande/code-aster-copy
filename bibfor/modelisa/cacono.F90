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
subroutine cacono(noma, ndim, llist1, llist2, no1, &
                  no2, norm1, norm2, inoma)
    implicit none
#include "jeveux.h"
#include "asterfort/canorm.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/pacoor.h"
#include "asterfort/panbno.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: ndim, no1, no2
    character(len=8) :: noma
    character(len=24) :: llist1, llist2
    real(kind=8) :: norm1(*), norm2(*)
!     BUT : CALCUL DES NORMALES AUX NOEUDS NO1 ET NO2
!     APPARTENANT AUX LISTES DE MAILLES LLIST1 ET LLIST2
!     CUMUL DE CES NORMALES
!     CALCUL DU JEU : SOMME_SUR_I(X1(I)*N1(I)+X2(I)*N2(I))
!
! IN  NOMA    K8  : NOM DU MAILLAGE
! IN  NDIM    I   : DIMENSION DU MODELE (2 SI COORD_2D OU 3 SI COORD_3D)
! IN  LLIST1  K8  : NOM DE LA LISTE DE MAILLE ASSOCIEE AU NOEUD NO1
! IN  LLIST2  K8  : NOM DE LA LISTE DE MAILLE ASSOCIEE AU NOEUD NO2
! IN  NO1     I   : NOEUD NUMERO 1 (NUMERO ABSOLU)
! IN  NO2     I   : NOEUD NUMERO 2 (NUMERO ABSOLU)
!
! OUT NORM1   R8  : NORMALE AU NOEUD 1
! OUT NORM2   R8  : NORMALE AU NOEUD 2
! OUT INOMA   I   :  = 0    SI LE NOEUD 1 N'APPARTIENT PAS A LLIST1
!                        ET SI LE NOEUD 2 N'APPARTIENT PAS A LLIST2
!                    = 1    SINON
!                    =-1    SI LE NOEUD 1 APPARTIENT A UNE MAILLE POI1
!                    =-2    SI LE NOEUD 2 APPARTIENT A UNE MAILLE POI1
! ROUTINES APPELEES :
!     CANORM      PACOOR      PANBNO
!
!
!
    real(kind=8) :: vecnor(3), coor(27)
    integer(kind=8) :: nbma, numma, ipoi1, inoma, i, ilist, ipoi2, ima, ityp, imad, ino
    integer(kind=8) :: nbno, nbnott(3), inorm, j
! DEBUT ----------------------------------------------------------------
!
!     INORM = 0 POUR VECTEUR NORMAL DE NORME "SURFACE"
!     INORM = 1 POUR VECTEUR NORMAL UNITAIRE
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call jemarq()
    inorm = 1
    do i = 1, 3
        norm1(i) = 0.0d0
        norm2(i) = 0.0d0
    end do
    call jeveuo(llist1, 'L', ilist)
    nbma = zi(ilist)
    inoma = 0
    ipoi1 = 0
    ipoi2 = 0
!
    do ima = 1, nbma
        ityp = zi(ilist+2*(ima-1)+2)
        call panbno(ityp, nbnott)
        nbno = nbnott(1)+nbnott(2)+nbnott(3)
        numma = zi(ilist+2*(ima-1)+1)
        call jeveuo(jexnum(noma//'.CONNEX', numma), 'L', imad)
        do ino = 1, nbno
            if (zi(imad-1+ino) .eq. no1) then
!           CAS D'UNE MAILLE POI1
                if (nbno .eq. 1) then
                    ipoi1 = 1
                    goto 100
                end if
                inoma = 1
                call pacoor(noma, numma, nbno, coor)
                call canorm(coor, vecnor, ndim, ityp, inorm)
                do j = 1, 3
                    norm1(j) = norm1(j)+vecnor(j)
                end do
                goto 100
!
            end if
!
        end do
100     continue
    end do
!
    call jeveuo(llist2, 'L', ilist)
    nbma = zi(ilist)
    do ima = 1, nbma
        ityp = zi(ilist+2*(ima-1)+2)
        call panbno(ityp, nbnott)
        nbno = nbnott(1)+nbnott(2)+nbnott(3)
        numma = zi(ilist+2*(ima-1)+1)
        call jeveuo(jexnum(noma//'.CONNEX', numma), 'L', imad)
        do ino = 1, nbno
            if (zi(imad-1+ino) .eq. no2) then
!           CAS D'UNE MAILLE POI1
                if (nbno .eq. 1) then
                    ipoi2 = 1
                    goto 200
                end if
                inoma = 1
                call pacoor(noma, numma, nbno, coor)
                call canorm(coor, vecnor, ndim, ityp, inorm)
                do j = 1, 3
                    norm2(j) = norm2(j)+vecnor(j)
                end do
                goto 200
!
            end if
!
        end do
200     continue
    end do
!
!     ON CHERCHE LE CAS OU NO2 APPARTIENT SEULEMENT A UNE MAILLE POI1
!     ON VERIFIE DANS CE CAS QUE NO1 N'APPARTIENT PAS A UNE MAILLE POI1
!
    if (ipoi1 .eq. 1) then
        inoma = -1
    else if (ipoi2 .eq. 1) then
        inoma = -2
    end if
    if ((ipoi1 .eq. 1) .and. (ipoi2 .eq. 1)) then
        call utmess('F', 'MODELISA2_43')
    end if
!
! FIN ------------------------------------------------------------------
    call jedema()
end subroutine
