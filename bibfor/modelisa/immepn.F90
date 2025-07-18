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
subroutine immepn(nbcnx, xyzma, x3dca, itetra, xbar, &
                  immer)
    implicit none
!  DESCRIPTION : TENTATIVE D'IMMERSION D'UN NOEUD CABLE X3DCA(3) DANS
!  -----------   UNE MAILLE PENTAEDRE APPARTENANT A LA STRUCTURE BETON
!                APPELANT : IMMENO
!
!  IN     : NBCNX  : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DE LA MAILLE PENTAEDRE DANS
!                    LAQUELLE EST TENTEE L'IMMERSION
!  IN     : XYZMA  : REAL*8 , TABLEAU DE DIMENSIONS (3,NNOMAX)
!                    TABLEAU DES COORDONNEES DES NOEUDS DE LA MAILLE
!                    PENTAEDRE DANS LAQUELLE EST TENTEE L'IMMERSION
!  IN     : X3DCA  : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU NOEUD CABLE
!  OUT    : ITETRA : INTEGER , SCALAIRE
!                    SI IMMERSION REUSSIE : INDICATEUR DU SOUS-DOMAINE
!                    TETRAEDRE AUQUEL APPARTIENT LE NOEUD CABLE
!                    ITETRA = 1 OU 2 OU 3
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
    integer(kind=8) :: idc, id(8), ii, j, ktest
    real(kind=8) :: d, dx, dy, dz
    integer(kind=8) :: f1(4), f2(4), f3(4)
!
    aster_logical :: facnp1, facnp2, facnp3
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
!     INDICATEUR SI DES FACES NON PLANES A 4 NOEUDS ONT DU ETRE
!     RENUMEROTEES AFIN QUE LA SURFACE DECRITE DEVIENNE ENVELOPPE
!     CONVEXE DU VOLUME TETRAEDRE
!
    facnp1 = .false.
    facnp2 = .false.
    facnp3 = .false.
!
!CCC    ORIENTATION PREMIERE FACE QUADRANGULAIRE NOEUDS 4-1-2-5
!
    f1(1) = 4
    f1(2) = 1
    f1(3) = 2
    f1(4) = 5
!
    call cotfac(xyzma, f1(1), f1(2), f1(3), 3, &
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
!CCC    ORIENTATION DEUXIEME FACE QUADRANGULAIRE NOEUDS 5-2-3-6
!
    f2(1) = 5
    f2(2) = 2
    f2(3) = 3
    f2(4) = 6
!
    call cotfac(xyzma, f2(1), f2(2), f2(3), 1, &
                xyzma(1, f2(4)), idc)
    if (idc .lt. 0) then
        ii = f2(4)
        f2(4) = f2(1)
        f2(1) = f2(2)
        f2(2) = f2(3)
        f2(3) = ii
        facnp2 = .true.
    end if
!
!CCC    ORIENTATION TROISIEME FACE QUADRANGULAIRE NOEUDS 3-1-4-6
!
    f3(1) = 3
    f3(2) = 1
    f3(3) = 4
    f3(4) = 6
!
    call cotfac(xyzma, f3(1), f3(2), f3(3), 2, &
                xyzma(1, f3(4)), idc)
    if (idc .lt. 0) then
        ii = f3(4)
        f3(4) = f3(1)
        f3(1) = f3(2)
        f3(2) = f3(3)
        f3(3) = ii
        facnp3 = .true.
    end if
!
!
!
    ii = 0
!CCC    POSITION COTE INTERNE PREMIERE FACE QUAD (2 PLANS)
    call cotfac(xyzma, f1(1), f1(2), f1(3), 3, &
                x3dca(1), id(1))
    if (id(1) .ge. 0) then
        ii = ii+1
        call cotfac(xyzma, f1(3), f1(4), f1(1), 3, &
                    x3dca(1), id(2))
!CCC    POSITION COTE INTERNE DEUXIEME FACE QUAD (2 PLANS)
        if (id(2) .ge. 0) then
            ii = ii+1
            call cotfac(xyzma, f2(1), f2(2), f2(3), 1, &
                        x3dca(1), id(3))
            if (id(3) .ge. 0) then
                ii = ii+1
                call cotfac(xyzma, f2(3), f2(4), f2(1), 1, &
                            x3dca(1), id(4))
!CCC    POSITION COTE INTERNE TROISIEME FACE QUAD (2 PLANS)
                if (id(4) .ge. 0) then
                    ii = ii+1
                    call cotfac(xyzma, f3(1), f3(2), f3(3), 2, &
                                x3dca(1), id(5))
                    if (id(5) .ge. 0) then
                        ii = ii+1
                        call cotfac(xyzma, f3(3), f3(4), f3(1), 2, &
                                    x3dca(1), id(6))
!CCC    POSITION COTE INTERNE QUATRIEME FACE TRIA (1 PLAN)
                        if (id(6) .ge. 0) then
                            ii = ii+1
                            call cotfac(xyzma, 4, 5, 6, 3, &
                                        x3dca(1), id(7))
!CCC    POSITION COTE INTERNE CINQUIEME FACE TRIA (1 PLAN)
                            if (id(7) .ge. 0) then
                                ii = ii+1
                                call cotfac(xyzma, 1, 2, 3, 4, &
                                            x3dca(1), id(8))
                                if (id(8) .ge. 0) ii = ii+1
                            end if
                        end if
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
    if (ii .lt. 8) then
!
        immer = -1
        goto 999
!
    else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ktest = 0
        do j = 1, 8
            ktest = ktest+id(j)
        end do
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      NOEUD IMMERGE DANS LE VOLUME DE LA MAILLE
!      CALCUL DES COORDONNES BARYCENTRIQUES SAUF SI
!      COINCIDENCE AVEC NOEUD MILIEU SI MAILLE HEXA20 OU HEXA27
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
        if (ktest .gt. 5) then
!
            if (nbcnx .eq. 15) then
                do j = 7, 15, 1
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
!   TEST D'APPARTENANCE A UN SOUS-DOMAINE TETRAEDRE, PAR DETERMINATION
!   DES COORDONNEES BARYCENTRIQUES (DECOUPAGE DU PENTA EN TROIS TETRAS)
!
!.....TETRAEDRE 2-3-4-5
!
            itetra = 1
            call tstbar(4, xyzma(1, 2), xyzma(1, 3), xyzma(1, 4), xyzma(1, 5), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.... TETRAEDRE 3-4-5-6
!
            itetra = 2
            call tstbar(4, xyzma(1, 3), xyzma(1, 4), xyzma(1, 5), xyzma(1, 6), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.... TETRAEDRE 1-2-3-4
!
            itetra = 3
            call tstbar(4, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), xyzma(1, 4), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!  DANS LE CAS DE FACES REORIENTEE (FACNP VRAI) ON TESTE LA PRESENCE
!  DANS LES PETITS TETRAEDRES DEFINIS PAR CHAQUE FACE.
!
!.... TETRAEDRE 4-1-2-5
!
            if (facnp1) then
                itetra = 4
                call tstbar(4, xyzma(1, 4), xyzma(1, 1), xyzma(1, 2), xyzma(1, 5), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 5-2-3-6
!
            if (facnp2) then
                itetra = 5
                call tstbar(4, xyzma(1, 5), xyzma(1, 2), xyzma(1, 3), xyzma(1, 6), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 3-1-4-6
!
            if (facnp3) then
                itetra = 6
                call tstbar(4, xyzma(1, 3), xyzma(1, 1), xyzma(1, 4), xyzma(1, 6), &
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
!
999 continue
!
! --- FIN DE IMMEPN.
end subroutine
