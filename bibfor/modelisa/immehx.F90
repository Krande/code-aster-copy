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
subroutine immehx(nbcnx, xyzma, x3dca, itetra, xbar, &
                  immer)
    implicit none
!  DESCRIPTION : TENTATIVE D'IMMERSION D'UN NOEUD CABLE X3DCA(3) DANS
!  -----------   UNE MAILLE HEXAEDRE APPARTENANT A LA STRUCTURE BETON
!                APPELANT : IMMENO
!
!  IN     : NBCNX  : INTEGER , SCALAIRE
!                    NOMBRE DE NOEUDS DE LA MAILLE HEXAEDRE DANS
!                    LAQUELLE EST TENTEE L'IMMERSION
!  IN     : XYZMA  : REAL*8 , TABLEAU DE DIMENSIONS (3,NNOMAX)
!                    TABLEAU DES COORDONNEES DES NOEUDS DE LA MAILLE
!                    HEXAEDRE DANS LAQUELLE EST TENTEE L'IMMERSION
!  IN     : X3DCA  : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU NOEUD CABLE
!  OUT    : ITETRA : INTEGER , SCALAIRE
!                    SI IMMERSION REUSSIE : INDICATEUR DU SOUS-DOMAINE
!                    TETRAEDRE AUQUEL APPARTIENT LE NOEUD CABLE
!                    ITETRA = 1 OU 2 OU 3 OU 4 OU 5 OU 6
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
    integer(kind=8) :: idc, id(12), ii, j, ktest
    real(kind=8) :: d, dx, dy, dz
    integer(kind=8) :: f1(4), f2(4), f3(4), f4(4), f5(4), f6(4)
    aster_logical :: facnp1, facnp2, facnp3, facnp4, facnp5, facnp6
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
    facnp4 = .false.
    facnp5 = .false.
    facnp6 = .false.
!
!CCC    ORIENTATION PREMIERE FACE NOEUDS 1-2-3-4
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
!CCC    ORIENTATION DEUXIEME FACE NOEUDS 3-4-8-7
!
    f2(1) = 3
    f2(2) = 4
    f2(3) = 8
    f2(4) = 7
!
    call cotfac(xyzma, f2(1), f2(2), f2(3), 5, &
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
!CCC    ORIENTATION TROISIEME FACE NOEUDS 6-7-8-5
!
    f3(1) = 6
    f3(2) = 7
    f3(3) = 8
    f3(4) = 5
!
    call cotfac(xyzma, f3(1), f3(2), f3(3), 1, &
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
!CCC    ORIENTATION QUATRIEME FACE NOEUDS 6-5-1-2
!
    f4(1) = 6
    f4(2) = 5
    f4(3) = 1
    f4(4) = 2
!
    call cotfac(xyzma, f4(1), f4(2), f4(3), 4, &
                xyzma(1, f4(4)), idc)
    if (idc .lt. 0) then
        ii = f4(4)
        f4(4) = f4(1)
        f4(1) = f4(2)
        f4(2) = f4(3)
        f4(3) = ii
        facnp4 = .true.
    end if
!
!CCC    ORIENTATION CINQUIEME FACE NOEUDS 6-2-3-7
!
    f5(1) = 6
    f5(2) = 2
    f5(3) = 3
    f5(4) = 7
!
    call cotfac(xyzma, f5(1), f5(2), f5(3), 5, &
                xyzma(1, f5(4)), idc)
    if (idc .lt. 0) then
        ii = f5(4)
        f5(4) = f5(1)
        f5(1) = f5(2)
        f5(2) = f5(3)
        f5(3) = ii
        facnp5 = .true.
    end if
!
!CCC    ORIENTATION SIXIEME FACE NOEUDS 1-5-8-4
!
    f6(1) = 1
    f6(2) = 5
    f6(3) = 8
    f6(4) = 4
!
    call cotfac(xyzma, f6(1), f6(2), f6(3), 2, &
                xyzma(1, f6(4)), idc)
    if (idc .lt. 0) then
        ii = f6(4)
        f6(4) = f6(1)
        f6(1) = f6(2)
        f6(2) = f6(3)
        f6(3) = ii
        facnp6 = .true.
    end if
!
    ii = 0
!CCC    POSITION COTE INTERNE PREMIERE FACE (2 PLANS)
    call cotfac(xyzma, f1(1), f1(2), f1(3), 5, &
                x3dca(1), id(1))
    if (id(1) .ge. 0) then
        ii = ii+1
        call cotfac(xyzma, f1(3), f1(4), f1(1), 5, &
                    x3dca(1), id(2))
!CCC    POSITION COTE INTERNE DEUXIEME FACE (2 PLANS)
        if (id(2) .ge. 0) then
            ii = ii+1
            call cotfac(xyzma, f2(1), f2(2), f2(3), 5, &
                        x3dca(1), id(3))
            if (id(3) .ge. 0) then
                ii = ii+1
                call cotfac(xyzma, f2(3), f2(4), f2(1), 5, &
                            x3dca(1), id(4))
!CCC    POSITION COTE INTERNE TROISIEME FACE (2 PLANS)
                if (id(4) .ge. 0) then
                    ii = ii+1
                    call cotfac(xyzma, f3(1), f3(2), f3(3), 1, &
                                x3dca(1), id(5))
                    if (id(5) .ge. 0) then
                        ii = ii+1
                        call cotfac(xyzma, f3(3), f3(4), f3(1), 1, &
                                    x3dca(1), id(6))
!CCC    POSITION COTE INTERNE QUATRIEME FACE (2 PLANS)
                        if (id(6) .ge. 0) then
                            ii = ii+1
                            call cotfac(xyzma, f4(1), f4(2), f4(3), 4, &
                                        x3dca(1), id(7))
                            if (id(7) .ge. 0) then
                                ii = ii+1
                                call cotfac(xyzma, f4(3), f4(4), f4(1), 4, &
                                            x3dca(1), id(8))
!CCC    POSITION COTE INTERNE CINQUIEME FACE (2 PLANS)
                                if (id(8) .ge. 0) then
                                    ii = ii+1
                                    call cotfac(xyzma, f5(1), f5(2), f5(3), 5, &
                                                x3dca(1), id(9))
                                    if (id(9) .ge. 0) then
                                        ii = ii+1
                                        call cotfac(xyzma, f5(3), f5(4), f5(1), 5, &
                                                    x3dca(1), id(10))
!CCC    POSITION COTE INTERNE SIXIEME FACE (2 PLANS)
                                        if (id(10) .ge. 0) then
                                            ii = ii+1
                                            call cotfac(xyzma, f6(1), f6(2), f6(3), 2, &
                                                        x3dca(1), id(11))
                                            if (id(11) .ge. 0) then
                                                ii = ii+1
                                                call cotfac(xyzma, f6(3), f6(4), f6(1), 2, &
                                                            x3dca(1), id(12))
                                                if (id(12) .ge. 0) ii = ii+1
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end if
    end if
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      NOEUD A L EXTERIEUR DU VOLUME DE LA MAILLE
!      ON A TROUVE : ON RESSORT COMPLETEMENT
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    if (ii .lt. 12) then
!
        immer = -1
        goto 999
!
    else
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ktest = 0
        do j = 1, 12
            ktest = ktest+id(j)
        end do
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      NOEUD IMMERGE DANS LE VOLUME DE LA MAILLE
!      CALCUL DES COORDONNES BARYCENTRIQUES SAUF SI
!      COINCIDENCE AVEC NOEUD MILIEU SI MAILLE HEXA20 OU HEXA27
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
        if (ktest .gt. 6) then
!
            if (nbcnx .eq. 20) then
                do j = 9, 20, 1
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
            if (nbcnx .eq. 27) then
                do j = 21, 27, 1
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
!     TEST D'APPARTENANCE A UN SOUS-DOMAINE TETRAEDRE PAR DETERMINATION
!     DES COORDONNEES BARYCENTRIQUES (DECOUPAGE D HEXA EN CINQ TETRAS)
!
!.....TETRAEDRE 6-3-8-1
!
            itetra = 1
            call tstbar(4, xyzma(1, 6), xyzma(1, 3), xyzma(1, 8), xyzma(1, 1), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.... TETRAEDRE 1-3-8-4
!
            itetra = 2
            call tstbar(4, xyzma(1, 1), xyzma(1, 3), xyzma(1, 8), xyzma(1, 4), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.... TETRAEDRE 6-8-1-5
!
            itetra = 3
            call tstbar(4, xyzma(1, 6), xyzma(1, 8), xyzma(1, 1), xyzma(1, 5), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.....TETRAEDRE 1-3-6-2
!
            itetra = 4
            call tstbar(4, xyzma(1, 1), xyzma(1, 3), xyzma(1, 6), xyzma(1, 2), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!.... TETRAEDRE 6-8-3-7
!
            itetra = 5
            call tstbar(4, xyzma(1, 6), xyzma(1, 8), xyzma(1, 3), xyzma(1, 7), &
                        x3dca(1), xbar(1), immer)
            if (immer .ge. 0) goto 999
!
!  DANS LE CAS DE FACES REORIENTEE (FACNP VRAI) ON TESTE LA PRESENCE
!  DANS LES PETITS TETRAEDRES DEFINIS PAR CHAQUE FACE.
!
!.... TETRAEDRE 1-2-3-4
!
            if (facnp1) then
                itetra = 6
                call tstbar(4, xyzma(1, 1), xyzma(1, 2), xyzma(1, 3), xyzma(1, 4), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 3-4-8-7
!
            if (facnp2) then
                itetra = 7
                call tstbar(4, xyzma(1, 3), xyzma(1, 4), xyzma(1, 8), xyzma(1, 7), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 6-7-8-5
!
            if (facnp3) then
                itetra = 8
                call tstbar(4, xyzma(1, 6), xyzma(1, 7), xyzma(1, 8), xyzma(1, 5), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 6-5-1-2
!
            if (facnp4) then
                itetra = 9
                call tstbar(4, xyzma(1, 6), xyzma(1, 5), xyzma(1, 1), xyzma(1, 2), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 6-2-3-7
!
            if (facnp5) then
                itetra = 10
                call tstbar(4, xyzma(1, 6), xyzma(1, 2), xyzma(1, 3), xyzma(1, 7), &
                            x3dca(1), xbar(1), immer)
                if (immer .ge. 0) goto 999
            end if
!
!.... TETRAEDRE 1-5-8-4
!
            if (facnp6) then
                itetra = 11
                call tstbar(4, xyzma(1, 1), xyzma(1, 5), xyzma(1, 8), xyzma(1, 4), &
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
! --- FIN DE IMMEHX.
end subroutine
