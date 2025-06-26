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
subroutine fonno5(noma, indic, noe, na, nb, &
                  ndim, nbnoel, indr, vnor, vdir)
    implicit none
#include "jeveux.h"
#include "asterfort/gdire3.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
    character(len=8) :: noma
    integer(kind=8) :: indic(4), noe(4, 4), na, nb, ndim, nbnoel, indr(2)
    real(kind=8) :: vnor(2, 3), vdir(2, 3)
!
!
!     ------------------------------------------------------------------
!     BUT : CALCUL DES VECTEURS DE LA BASE LOCALE :
!             - VNOR : VECTEUR NORMAL A LA SURFACE DE LA FISSURE
!             - VDIR : VECTEUR DANS LA DIRECTION DE PROPAGATION
!           CAS OU CONFIG_INIT = COLLEE
!
!           RQ : CHACUN CONTIENT EN FAIT 2 VECTEURS (UN PAR LEVRE)
!     ------------------------------------------------------------------
!
! ENTREES
!     NOMA   : NOM DU MAILLAGE
!     INDIC  : INDICE DES FACES INTERNES
!     NBNOFF : NOMBRE DE NOEUD EN FOND DE FISSURE
!     NOE    : NOEUDS DES FACES CONTENANT NA et NB ET APPARTENANT AUX
!              MAILLES CONNECTEES AU NOEUD SOMMET COURANT
!              ET AUX LEVRES
!     NA     : NUMERO DU NOEUD SOMMET COURANT
!     NB     : NUMERO DU NOEUD SOMMET SUIVANT
!
! SORTIES
!     NBNOEL : NOMBRE DE NOEUDS SOMMETS PAR ELEMENTS
!     INDR   : INDICES DES FACES LIBRES
!     VNOR   : VECTEUR NORMAL A LA SURFACE DE LA FISSURE
!     VDIR   : VECTEUR DANS LA DIRECTION DE PROPAGATION
!
!     ----------------------------------------------------
!
    integer(kind=8) :: compte
    integer(kind=8) :: indice, inp, ino1, ino2, ico
    integer(kind=8) :: m(8)
    real(kind=8) :: vect1(3), vect2(3), vect3(3), vect4(3), norm1
    real(kind=8) :: coord(3, 4)
    real(kind=8), pointer :: vale(:) => null()
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
!     RECUPERATION DE L'ADRESSE DES COORDONNEES
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
    compte = 0
    indice = 0
    do inp = 1, 4
        if ((inp .ne. indic(1)) .and. (inp .ne. indic(2)) .and. (inp .ne. indic(3)) .and. &
            (inp .ne. indic(4))) then
            compte = compte+1
!         RENUMEROTATION LOCALE POUR AVOIR DANS LES DEUX PREMIERS
!         NOEUDS LES NOEUDS DU FOND DE FISSURE
!         CELA PERMET DEFINIR LE VECTEUR NORMAL A L'AIDE DE GDIRE3
            if (ndim .eq. 3) then
                if (noe(inp, 4) .eq. 0) then
                    nbnoel = 3
                    indice = 1
                else
                    nbnoel = 4
                    indice = 2
                end if
                do ino1 = 1, nbnoel
                    m(ino1) = noe(inp, ino1)
                    m(ino1+nbnoel) = noe(inp, ino1)
                end do
                do ino1 = 1, nbnoel
                    if (m(ino1) .eq. na) then
                        if (m(ino1+1) .eq. nb) then
                            do ino2 = 1, nbnoel
                                noe(inp, ino2) = m(ino1-1+ino2)
                            end do
                        else
                            do ino2 = 1, nbnoel
                                noe(inp, ino2) = m(ino1+1+nbnoel-ino2)
                            end do
                        end if
                    end if
                end do
                do ico = 1, 3
                    do ino1 = 1, nbnoel
                        coord(ico, ino1) = vale((noe(inp, ino1)-1)*3+ico)
                    end do
                    vect1(ico) = coord(ico, 2)-coord(ico, 1)
                end do
                call normev(vect1, norm1)
!           CALCUL DU VECTEUR DIRECTION DE PROPAGATION
                call gdire3(coord, vect4(1), vect4(2), vect4(3), indice)
!
!           CALCUL DU VECTEUR NORMAL A LA FACE
                call provec(vect4, vect1, vect3)
                call normev(vect3, norm1)
                do ico = 1, 3
                    vnor(compte, ico) = vect3(ico)
                    vdir(compte, ico) = vect4(ico)
                end do
            else
!           LE NOEUD DU FOND DOIT ETRE LE PREMIER
                nbnoel = 2
                m(1) = noe(inp, 1)
                m(2) = noe(inp, 2)
                if (m(1) .ne. na) then
                    noe(inp, 1) = m(2)
                    noe(inp, 2) = m(1)
                end if
!           CALCUL DU VECTEUR DIRECTION DE PROPAGATION
                do ico = 1, 3
                    do ino1 = 1, 2
                        coord(ico, ino1) = vale((noe(inp, ino1)-1)*3+ico)
                    end do
                    vect1(ico) = coord(ico, 1)-coord(ico, 2)
                end do
                call normev(vect1, norm1)
!           CALCUL DU VECTEUR NORMAL A L'ARETE
                vect3(:) = 0.d0
                vect3(3) = 1.d0
                call provec(vect3, vect1, vect2)
                do ico = 1, 3
                    vnor(compte, ico) = vect2(ico)
                    vdir(compte, ico) = vect1(ico)
                end do
            end if
            indr(compte) = inp
        end if
    end do
    call jedema()
end subroutine
