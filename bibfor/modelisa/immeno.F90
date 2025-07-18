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
subroutine immeno(ncncin, nmabet, mailla, x3dca, noebe, &
                  numail, nbcnx, cxma, xyzma, itetra, &
                  xbar, immer)
    implicit none
!  DESCRIPTION : TENTATIVE D'IMMERSION D'UN NOEUD CABLE DANS LES MAILLES
!  -----------   APPARTENANT A LA STRUCTURE BETON
!                APPELANT : IMMECA
!
!  IN     : NCNCIN : CHARACTER*24 ,
!                    OBJET CONTENANT LA CONNECTIVITE INVERSE POUR LE
!                    GROUP_MA BETON
!  IN     : NMABET : CHARACTER*24 ,
!                    OBJET CONTENANT LES MAILLES BETON
!  IN     : MAILLA : CHARACTER*8 , SCALAIRE
!                    NOM DU CONCEPT MAILLAGE ASSOCIE A L'ETUDE
!  IN     : X3DCA  : REAL*8 , VECTEUR DE DIMENSION 3
!                    COORDONNEES DU NOEUD CABLE CONSIDERE
!  IN     : NOEBE  : INTEGER , SCALAIRE
!                    NUMERO DU NOEUD BETON LE PLUS PROCHE DU NOEUD CABLE
!                    CONSIDERE
!  OUT    : NUMAIL : INTEGER , SCALAIRE
!                    SI IMMERSION REUSSIE : NUMERO DE LA MAILLE DANS
!                    LAQUELLE EST REALISEE L'IMMERSION
!  OUT    : NBCNX  : INTEGER , SCALAIRE
!                    SI IMMERSION REUSSIE : NOMBRE DE NOEUDS DE LA
!                    MAILLE DANS LAQUELLE EST REALISEE L'IMMERSION
!  OUT    : CXMA   : INTEGER , VECTEUR DE DIMENSION AU PLUS NNOMAX
!                    SI IMMERSION REUSSIE : NUMEROS DES NOEUDS DE LA
!                    MAILLE DANS LAQUELLE EST REALISEE L'IMMERSION
!                    (TABLE DE CONNECTIVITE)
!  OUT    : XYZMA  : REAL*8 , TABLEAU DE DIMENSIONS (3,NNOMAX)
!                    SI IMMERSION REUSSIE : TABLEAU DES COORDONNEES
!                    DES NOEUDS DE LA MAILLE DANS LAQUELLE EST REALISEE
!                    L'IMMERSION
!  OUT    : ITETRA : INTEGER , SCALAIRE
!                    SI IMMERSION REUSSIE : INDICATEUR DU SOUS-DOMAINE
!                    TETRAEDRE AUQUEL APPARTIENT LE NOEUD CABLE
!                    ITETRA = 1            SI IMMERSION DANS UNE
!                                          MAILLE TETRAEDRE
!                    ITETRA = 1 OU 2       SI IMMERSION DANS UNE
!                                          MAILLE PYRAMIDE
!                    ITETRA = 1 OU 2 OU 3  SI IMMERSION DANS UNE
!                                          MAILLE PENTAEDRE
!                    ITETRA = 1 OU 2 OU 3  SI IMMERSION DANS UNE
!                          OU 4 OU 5 OU 6  MAILLE HEXAEDRE
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
#include "jeveux.h"
#include "asterfort/immehx.h"
#include "asterfort/immepn.h"
#include "asterfort/immepy.h"
#include "asterfort/immett.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
!
!
! ARGUMENTS
! ---------
    character(len=8) :: mailla
    character(len=24) :: ncncin, nmabet
    integer(kind=8) :: noebe, numail, numai0, nbcnx, cxma(*), itetra, immer
    real(kind=8) :: x3dca(*), xyzma(3, *), xbar(*)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: imail, inoma, jcoor, jcxma, noe, jdrvlc, jcncin, iadr, nbm
    integer(kind=8) :: jlimab
    character(len=24) :: conxma, coorno
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 1   ACCES AUX OBJETS DU CONCEPT MAILLAGE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    conxma = mailla//'.CONNEX'
    coorno = mailla//'.COORDO    .VALE'
    call jeveuo(coorno, 'L', jcoor)
    call jeveuo(jexatr(ncncin, 'LONCUM'), 'L', jdrvlc)
    call jeveuo(jexnum(ncncin, 1), 'L', jcncin)
    call jeveuo(nmabet, 'L', jlimab)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 2   TENTATIVE D'IMMERSION DU NOEUD CABLE CONSIDERE DANS LES MAILLES
!     APPARTENANT A LA STRUCTURE BETON
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!.... BOUCLE SUR LES MAILLES APPARTENANT A LA STRUCTURE BETON, POUR
!.... RETROUVER LE NOEUD BETON LE PLUS PROCHE DANS LES CONNECTIVITES
!
    nbm = zi(jdrvlc+noebe+1-1)-zi(jdrvlc+noebe-1)
    iadr = zi(jdrvlc+noebe-1)
    do imail = 1, nbm
        numai0 = zi(jcncin+iadr-1+imail-1)
        numail = zi(jlimab+numai0-1)
        call jelira(jexnum(conxma, numail), 'LONMAX', nbcnx)
        call jeveuo(jexnum(conxma, numail), 'L', jcxma)
!
!........RECUPERATION DES NUMEROS ET DES COORDONNEES DES NOEUDS
!........DE LA MAILLE
!
        do inoma = 1, nbcnx
            noe = zi(jcxma+inoma-1)
            cxma(inoma) = noe
            xyzma(1, inoma) = zr(jcoor+3*(noe-1))
            xyzma(2, inoma) = zr(jcoor+3*(noe-1)+1)
            xyzma(3, inoma) = zr(jcoor+3*(noe-1)+2)
        end do
!
!........TEST D'APPARTENANCE DU NOEUD CABLE AU DOMAINE GEOMETRIQUE
!........DEFINI PAR LA MAILLE
!
        if ((nbcnx .eq. 4) .or. (nbcnx .eq. 10)) then
            call immett(nbcnx, xyzma(1, 1), x3dca(1), itetra, xbar(1), &
                        immer)
        else if ((nbcnx .eq. 5) .or. (nbcnx .eq. 13)) then
            call immepy(nbcnx, xyzma(1, 1), x3dca(1), itetra, xbar(1), &
                        immer)
        else if ((nbcnx .eq. 6) .or. (nbcnx .eq. 15)) then
            call immepn(nbcnx, xyzma(1, 1), x3dca(1), itetra, xbar(1), &
                        immer)
        else
            call immehx(nbcnx, xyzma(1, 1), x3dca(1), itetra, xbar(1), &
                        immer)
        end if
!
        if (immer .ge. 0) goto 999
!
    end do
!
999 continue
    call jedema()
!
! --- FIN DE IMMENO.
end subroutine
