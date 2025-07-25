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
subroutine immepy(nbcnx, xyzma, x3dca, itetra, xbar, &
                  immer)
    implicit none
!  DESCRIPTION : TENTATIVE D'IMMERSION D'UN NOEUD CABLE X3DCA(3) DANS
!  -----------   UNE MAILLE PYRAMIDE APPARTENANT A LA STRUCTURE BETON
!                APPELANT : IMMENO
!
!  IN     : NBCNX  : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DE LA MAILLE PYRAMIDE DANS
!                    LAQUELLE EST TENTEE L'IMMERSION
!  IN     : XYZMA  : REAL*8 , TABLEAU DE DIMENSIONS (3,NNOMAX)
!                    TABLEAU DES COORDONNEES DES NOEUDS DE LA MAILLE
!                    PYRAMIDE DANS LAQUELLE EST TENTEE L'IMMERSION
!  IN     : X3DCA  : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU NOEUD CABLE
!  OUT    : ITETRA : INTEGER , SCALAIRE
!                    SI IMMERSION REUSSIE : INDICATEUR DU SOUS-DOMAINE
!                    TETRAEDRE AUQUEL APPARTIENT LE NOEUD CABLE
!                    ITETRA = 1 OU 2
!  OUT    : XBAR   : REAL*8 , VECTEUR DE DIMENSION 4
!                    SI IMMERSION REUSSIE : COORDONNEES BARYCENTRIQUES
!                    DU NOEUD CABLE DANS LE SOUS-DOMAINE TETRAEDRE
!                    AUQUEL IL APPARTIENT
!  OUT    : IMMER  : INTEGER , SCALAIRE
!                    INDICE D'IMMERSION
!                    IMMER = -1  IMMERSION NON REUSSIE
!                    IMMER =  0  LE NOEUD CABLE EST A L'INTERIEUR
!                                DE LA MAILLE
!                    IMMER = 100 + 10 * NUMERO DE FACE
!                                LE NOEUD CABLE EST SUR UNE FACE
!                                DE LA MAILLE
!                    IMMER = 100 + 10 * NUMERO DE FACE + NUMERO D'ARETE
!                                LE NOEUD CABLE EST SUR UNE ARETE
!                                DE LA MAILLE
!                    IMMER =  2  LE NOEUD CABLE COINCIDE AVEC UN DES
!                                NOEUDS DE LA MAILLE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/cotfac.h"
#include "asterfort/tstbar.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbcnx, itetra, immer
    real(kind=8) :: xyzma(3, *), x3dca(*), xbar(*)
!
! VARIABLES LOCALES
! -----------------
    real(kind=8) :: d, dx, dy, dz
    integer(kind=8) :: ktest, f1(4), idc, id(6), ii, j
    aster_logical :: facnp1
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    facnp1 = .false.
!
!CCC    ORIENTATION FACE QUADRANGULAIRE NOEUDS 1-2-3-4
!
    f1(1) = 1
    f1(2) = 2
    f1(3) = 3
    f1(4) = 4
!
    call cotfac(xyzma, f1(1), f1(2), f1(3), 5, &
                xyzma(1, f1(4)), idc)
    if (idc .lt. 0) then
        ii = f1(4)
        f1(4) = f1(1)
        f1(1) = f1(2)
        f1(2) = f1(3)
        f1(3) = ii
        facnp1 = .true.
    end if
!
!
    ii = 0
!CCC    POSITION COTE INTERNE PREMIERE FACE (2 PLANS)
    call cotfac(xyzma, f1(1), f1(2), f1(3), 5, &
                x3dca(1), id(1))
    if (id(1) .ge. 0) then
        ii = ii+1
        call cotfac(xyzma, f1(3), f1(4), f1(1), 5, &
                    x3dca(1), id(2))
!
!CCC    POSITION COTE INTERNE DEUXIEME FACE (1 PLAN)
        if (id(2) .ge. 0) then
            ii = ii+1
            call cotfac(xyzma, 1, 2, 5, 3, &
                        x3dca(1), id(3))
!
!CCC    POSITION COTE INTERNE TROISIEME FACE (1 PLAN)
            if (id(3) .ge. 0) then
                ii = ii+1
                call cotfac(xyzma, 2, 3, 5, 4, &
                            x3dca(1), id(4))
!
!CCC    POSITION COTE INTERNE QUATRIEME FACE (1 PLAN)
                if (id(4) .ge. 0) then
                    ii = ii+1
                    call cotfac(xyzma, 3, 4, 5, 1, &
                                x3dca(1), id(5))
!
!CCC    POSITION COTE INTERNE CINQUIEME FACE (1 PLAN)
                    if (id(5) .ge. 0) then
                        ii = ii+1
                        call cotfac(xyzma, 4, 1, 5, 2, &
                                    x3dca(1), id(6))
                        if (id(6) .ge. 0) ii = ii+1
                    end if
                end if
            end if
        end if
    end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      NOEUD A L EXTERIEUR DU VOLUME DE LA MAILLE
!      ON A TROUVE : ON RESSORT COMPLETEMENT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (ii .lt. 6) then
!
        immer = -1
        goto 999
!
    else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ktest = 0
        do j = 1, 6
            ktest = ktest+id(j)
        end do
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      NOEUD IMMERGE DANS LE VOLUME DE LA MAILLE
!      CALCUL DES COORDONNES BARYCENTRIQUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
        if (ktest .gt. 3) then
!
            if (nbcnx .eq. 13) then
                do j = 6, 13, 1
                    dx = xyzma(1, j)-x3dca(1)
                    dy = xyzma(2, j)-x3dca(2)
                    dz = xyzma(3, j)-x3dca(3)
                    d = dx*dx+dy*dy+dz*dz
                    if (d .lt. r8prem()) then
                        immer = 2
                        goto 999
                    end if
                end do
            end if
!
!     TEST D'APPARTENANCE A UN SOUS-DOMAINE TETRAEDRE, PAR DETERMINATION
!     DES COORDONNEES BARYCENTRIQUES (DECOUPAGE DU PYRA EN DEUX TETRAS)
!
!.....TETRAEDRE 1-2-3-5
!
            itetra = 1
            call tstbar(4, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), xyzma(1, 5), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.... TETRAEDRE 1-3-4-5
!
            itetra = 2
            call tstbar(4, xyzma(1, 1), xyzma(1, 3), xyzma(1, 4), xyzma(1, 5), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!  DANS LE CAS DE FACES REORIENTEE (FACNP1 VRAI) ON TESTE LA PRESENCE
!  DANS LES PETITS TETRAEDRES DEFINIS PAR CHAQUE FACE.
!
!.... TETRAEDRE 1-2-3-4
!
            if (facnp1) then
                itetra = 3
                call tstbar(4, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), xyzma(1, 4), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
            if (immer .lt. 0) then
                call utmess('F', 'MODELISA4_72')
            end if
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      NOEUD COINCIDANT AVEC UN NOEUD SOMMET -APPARTIENT A + DE 3 PLANS
!      ON A TROUVE : ON RESSORT COMPLETEMENT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
        else
!
            immer = 2
            goto 999
!
        end if
    end if
!
999 continue
!
! --- FIN DE IMMEPY.
end subroutine
