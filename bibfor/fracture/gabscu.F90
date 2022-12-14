! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine gabscu(lobj2, coorn, nomno, fond, xl,&
                  absgam)
    implicit none
!
!     ----------------------------------------------------------------
! FONCTION REALISEE:
!
!     POUR CHAQUE NOEUD DU FOND DE FISSURE GAMM0 ON CALCULE
!     SON ABSCISSE CURVILIGNE
!     TRAITEMENT PARTICULIER SI LA COURBE EST FERMEE
!
!     ------------------------------------------------------------------
! ENTREE:
!        LOBJ2  : NOMBRE DE NOEUD DE GAMM0
!        COORN  : NOM DE L'OBJET CONTENANT LES COORDONNEES DES NOEUDS
!        NOMNO  : NOM DE L'OBJET CONTENANT LES NOMS DES NOEUDS
!        FOND   : NOMS DES NOEUDS DU FOND DE FISSURE
!
! SORTIE:
!        XL     : LONGUEUR DE LA FISSURE
!        ABSGAM : ABSCISSE CURVILIGNE DES NOEUDS DU FOND DE FISSURE
!     ------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: nomno, coorn, numgam, absgam, fond
!
    integer :: lobj2, iadrco, iadrno, iadnum, iadabs
!
    real(kind=8) :: xi1, yi1, zi1, xj1, yj1, zj1, xij, yij, zij, xl
!
!
!-----------------------------------------------------------------------
    integer :: i, iret, j
!-----------------------------------------------------------------------
    call jemarq()
    call jeveuo(coorn, 'L', iadrco)
    call jeveuo(fond, 'L', iadrno)
!
! CALCUL DE LA LONGUEUR DU FOND DE FISSURE
!
! RECUPERATION DES NUMEROS DE NOEUDS DE GAMM0
!
!
    numgam = '&&LEGEND.NUMGAMM0'
    call wkvect(numgam, 'V V I', lobj2, iadnum)
    do j = 1, lobj2
        call jenonu(jexnom(nomno, zk8(iadrno+j-1)), zi(iadnum+j-1))
    end do
!
    xl = 0.d0
    do j = 1, lobj2-1
        xi1 = zr(iadrco+(zi(iadnum+j-1) -1)*3+1-1)
        yi1 = zr(iadrco+(zi(iadnum+j-1) -1)*3+2-1)
        zi1 = zr(iadrco+(zi(iadnum+j-1) -1)*3+3-1)
        xj1 = zr(iadrco+(zi(iadnum+j+1-1)-1)*3+1-1)
        yj1 = zr(iadrco+(zi(iadnum+j+1-1)-1)*3+2-1)
        zj1 = zr(iadrco+(zi(iadnum+j+1-1)-1)*3+3-1)
        xij = xj1-xi1
        yij = yj1-yi1
        zij = zj1-zi1
        xl = xl + sqrt(xij*xij + yij *yij +zij*zij)
    end do
!
!  CALCUL DE L'ABSCISSE CURVILIGNE DE CHAQUE NOEUD DE GAMM0
!
    absgam = '&&LEGEND.ABSGAMM0'
    call jeexin(absgam, iret)
    if (iret .eq. 0) then
        call wkvect(absgam, 'V V R', lobj2, iadabs)
!
        zr(iadabs) = 0.d0
        do i = 1, lobj2-1
            xi1 = zr(iadrco+(zi(iadnum+i-1) -1)*3+1-1)
            yi1 = zr(iadrco+(zi(iadnum+i-1) -1)*3+2-1)
            zi1 = zr(iadrco+(zi(iadnum+i-1) -1)*3+3-1)
            xj1 = zr(iadrco+(zi(iadnum+i+1-1)-1)*3+1-1)
            yj1 = zr(iadrco+(zi(iadnum+i+1-1)-1)*3+2-1)
            zj1 = zr(iadrco+(zi(iadnum+i+1-1)-1)*3+3-1)
            xij = xj1-xi1
            yij = yj1-yi1
            zij = zj1-zi1
            zr(iadabs+i+1-1) = zr(iadabs+i-1)+ sqrt(xij*xij + yij * yij + zij*zij)
        end do
    endif
!
! DESTRUCTION DES OBJETS DE TRAVAIL
!
    call jedetr(numgam)
!
    call jedema()
end subroutine
