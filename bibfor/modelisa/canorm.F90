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
subroutine canorm(coor, normal, ndim, ityp, inorm)
    implicit none
#include "jeveux.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: ndim, ityp, inorm
    real(kind=8) :: coor(*), normal(3)
!
!     BUT : CALCUL DE LA NORMALE A UNE MAILLE  EN UN NOEUD
!     AVEC OU SANS NORMALISATION DE CE VECTEUR
!
! IN  COOR    R8 : TABLEAU DES COORDONNEES DES NBNO NOEUDS DE LA MAILLE
!                  DE DIMENSION (3*NBNO)
! IN  NDIM    I  : DIMENSION DU MODELE (2 SI COORD_2D OU 3 SI COORD_3D)
! IN  ITYP    I  : TYPE DE LA MAILLE
! IN  INORM   I  : INDICATEUR DE NORMALISATION
!                  INORM = 0 PAS DE NORMALISATION
!                  INORM = 1 NORMALISATION
! OUT NORMALE R8 : NORMALE CALCULEE
!
! ROUTINES APPELLEES :
!     NORMEV     PROVEC
!     JENUNO     JEXNUM
!
!
!
    real(kind=8) :: xx(3), yy(3), norme, surf, vect(3)
    character(len=8) :: nomtm
!
! DEBUT ----------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: j
!-----------------------------------------------------------------------
    call jenuno(jexnum('&CATA.TM.NBNO', ityp), nomtm)
    if (nomtm(1:3) .eq. 'SEG') then
        if (ndim .eq. 3) then
            normal(1) = coor(4)-coor(1)
            normal(2) = coor(5)-coor(2)
            normal(3) = coor(6)-coor(3)
            if (inorm .eq. 1) then
                call normev(normal, norme)
            else
                call utmess('F', 'MODELISA3_20')
            end if
!          ELSE
!            NORMAL(1) = COOR(5) - COOR(2)
!            NORMAL(2) = COOR(1) - COOR(4)
!            NORMAL(3) = 0.0D0
!            IF (INORM.EQ.1) THEN
!              CALL NORMEV(NORMAL,NORME)
!            ENDIF
!
!          ENDIF
!
        else
            normal(1) = coor(5)-coor(2)
            normal(2) = coor(1)-coor(4)
            normal(3) = 0.0d0
            if (inorm .eq. 1) then
                call normev(normal, norme)
            end if
!
        end if
!
    else if (nomtm(1:4) .eq. 'TRIA') then
        if (ndim .eq. 2) then
            call utmess('F', 'MODELISA3_21')
!
        else
            do j = 1, 3
                xx(j) = coor(3+j)-coor(j)
                yy(j) = coor(6+j)-coor(3+j)
            end do
            call provec(xx, yy, normal)
            do j = 1, 3
                normal(j) = normal(j)/2.0d0
            end do
            if (inorm .eq. 1) then
                call normev(normal, norme)
            end if
!
        end if
!
    else if (nomtm(1:4) .eq. 'QUAD') then
        if (ndim .eq. 2) then
            call utmess('F', 'MODELISA3_22')
!
        else
!
!     PRODUIT VECTORIEL (N1N3) * (N2N4) POUR CALCULER LE VECTEUR NORMAL
!
            do j = 1, 3
                xx(j) = coor(6+j)-coor(j)
                yy(j) = coor(9+j)-coor(3+j)
            end do
            call provec(xx, yy, normal)
            call normev(normal, norme)
            if (inorm .eq. 0) then
!
!     ON CALCULE UNE APPROXIMATION DE LA SURFACE
!     DANS L'ORDRE
!     (N1N2) * (N1N3)
!     (N1N3) * (N1N4)
!     (N2N3) * (N2N4)
!     (N2N4) * (N2N1)
!
                surf = 0.0d0
                do j = 1, 3
                    xx(j) = coor(3+j)-coor(j)
                    yy(j) = coor(6+j)-coor(j)
                end do
                call provec(xx, yy, vect)
                call normev(vect, norme)
                surf = surf+norme
                do j = 1, 3
                    xx(j) = coor(6+j)-coor(j)
                    yy(j) = coor(9+j)-coor(j)
                end do
                call provec(xx, yy, vect)
                call normev(vect, norme)
                surf = surf+norme
                do j = 1, 3
                    xx(j) = coor(6+j)-coor(3+j)
                    yy(j) = coor(9+j)-coor(3+j)
                end do
                call provec(xx, yy, vect)
                call normev(vect, norme)
                surf = surf+norme
                do j = 1, 3
                    xx(j) = coor(9+j)-coor(3+j)
                    yy(j) = coor(j)-coor(3+j)
                end do
                call provec(xx, yy, vect)
                call normev(vect, norme)
                surf = surf+norme
                surf = surf/4.0d0
                do j = 1, 3
                    normal(j) = normal(j)*surf
                end do
            end if
!
        end if
!
    end if
!
!
! FIN ------------------------------------------------------------------
end subroutine
