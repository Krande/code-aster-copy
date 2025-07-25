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
subroutine xvoise(nnotot, nse, nnop, nno, jcnset, &
                  cninv, cvoise)
    implicit none
#include "jeveux.h"
    integer(kind=8) :: nvois
    parameter(nvois=3)
    integer(kind=8) :: nnotot, nse, nnop, nno, jcnset
    integer(kind=8) :: cninv(nnotot, nse+1), cvoise(nvois, nse)
!
!     BUT:
!         RECHERCHE DES VOISINS DES SOUS ELEMENTS DE
!         L'ELEMENT XFEM PARENT (EN 2D), PUIS ECRITURE DANS ZI()
!
!
!     ARGUMENTS:
!     ----------
!
!      ENTREE :
!-------------
! IN   NDIM   : DIMENSION
! IN   NNOTOT : NOMBRE TOTAL DE NOEUDS (POINTS D'INTERSECTION INCLUS)
! IN   NSE    : NOMBRE TOTAL DE SOUS ELEMENT DE L'ELEMENT PARENT
! IN   NNOP   : NOMBRE DE NOEUDS DE L'ELEMENT PARENT (POINTS)
!                D'INTERSECTION EXCLUS
! IN   NNO    : NOMBRE DE NOEUDS DU SOUS-ELEMENT DE REFERENCE
! IN   JCNSET : ADRESSE DANS ZI DE LA CONNECTIVITE DES SOUS-ELEMENTS
! IN   CNINV  : TABLEAU DE LA CONNECTIVITE INVERSE
!
!      SORTIE :
!-------------
! OUT  CVOISE : TABLEAU DES VOISINS PAR SOUS-ELEMENT
!
! ......................................................................
!
!
!
!
    integer(kind=8) :: nbmav1, nbmav2
    integer(kind=8) :: ise, in, ino, jno, imav1, numav1, indma1
    integer(kind=8) :: insui, inosui, jnosui, imav2, numav2, indma2
!
! ----------------------------------------------------------------------
! -----------------  BOUCLE SUR LES NSE SIMPLEXES  ------------------
! ----------------------------------------------------------------------
!
    do ise = 1, nse
!
! --------------------------------------------------------------------
! ------------  BOUCLE SUR LES SOMMETS DU SOUS-ELEMENTS  -------------
! --------------------------------------------------------------------
!
        do in = 1, nno
!
! ------- RECUPERATION DE LA NUMEROTATION XFEM
!
            ino = zi(jcnset-1+nno*(ise-1)+in)
!
! ------- NUMEROTATION PROPRE A LA CONNECTIVITE INVERSE
!
            if (ino .lt. 1000) then
                jno = ino
            else
                jno = ino-1000+nnop
            end if
!
            nbmav1 = cninv(jno, 1)
!
! --------------------------------------------------------------------
! -------  BOUCLE SUR LES VOISINS POTENTIELS CONTENANT "JNO"  --------
! --------------------------------------------------------------------
!
            do imav1 = 1, nbmav1
!
                indma1 = imav1+1
                numav1 = cninv(jno, indma1)
!
! --------- ON S'ASSURE QUE LE VOISIN POTENTIEL N'EST PAS LE
! --------- SOUS-ELEMENT COURANT
!
                if (numav1 .ne. ise) then
!
!
! ----------- RESPECT DE LA NUMEROTATION AU SEIN DU SOUS-ELEMENT
!
                    if (in .eq. nno) then
                        insui = 1
                    else
                        insui = in+1
                    end if
!
! ----------- RECUPERATION DE LA NUMEROTATION XFEM
!
                    inosui = zi(jcnset-1+nno*(ise-1)+insui)
!
! ----------- NUMEROTATION PROPRE A LA CONNECTIVITE INVERSE
!
                    if (inosui .lt. 1000) then
                        jnosui = inosui
                    else
                        jnosui = inosui-1000+nnop
                    end if
!
                    nbmav2 = cninv(jnosui, 1)
!
! --------------------------------------------------------------------
! ----  BOUCLE SUR LES VOISINS POTENTIELS CONTENANT "JNOSUI"  --------
! --------------------------------------------------------------------
!
                    do imav2 = 1, nbmav2
!
                        indma2 = imav2+1
                        numav2 = cninv(jnosui, indma2)
!
! ------------- ON LOCALISE LE VOISIN SITUE EN VIS-À-VIS DE L'ARRETE
! ------------- [JNO,JNOSUI] (S'IL EXISTE), PUIS ON L'ECRIT DANS ZI()
!
                        if (numav2 .eq. numav1) then
                            cvoise(in, ise) = numav1
                        end if
!
                    end do
!
                end if
!
            end do
!
        end do
!
    end do
!
end subroutine
