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
subroutine cupivo(xjvmax, indic, nbliac, ajliai, spliai, &
                  spavan, deficu, resocu)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/cudisi.h"
#include "asterfort/cuimp2.h"
#include "asterfort/cutabl.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=24) :: resocu
    character(len=24) :: deficu
    real(kind=8) :: xjvmax
    integer(kind=8) :: nbliac
    integer(kind=8) :: indic
    integer(kind=8) :: ajliai
    integer(kind=8) :: spliai
    integer(kind=8) :: spavan
!
! ----------------------------------------------------------------------
!
! ROUTINE LIAISON_UNILATERALE (RESOLUTION)
!
! ELIMINATION DES PIVOTS NULS DANS LA MATRICE DE CONTACT
!
!
! ----------------------------------------------------------------------
!
!
! IN  XJVMAX : VALEUR DU PIVOT MAX
! OUT INDIC  : +1 ON A RAJOUTE UNE LIAISON
!              -1 ON A ENLEVE UNE LIAISON
! I/O NBLIAC : NOMBRE DE LIAISONS ACTIVES
! I/O AJLIAI : INDICE DANS LA LISTE DES LIAISONS ACTIVES DE LA DERNIERE
!              LIAISON CORRECTE DU CALCUL
!              DE LA MATRICE DE CONTACT ACM1AT
! I/O SPLIAI : INDICE DANS LA LISTE DES LIAISONS ACTIVES DE LA DERNIERE
!              LIAISON AYANT ETE CALCULEE POUR LE VECTEUR CM1A
! IN  SPAVAN : INDICE DE DEBUT DE TRAITEMENT DES LIAISONS
! IN  DEFICU : SD DE DEFINITION (ISSUE D'AFFE_CHAR_MECA)
! IN  RESOCU : SD DE TRAITEMENT
!
!
!
!
!
    character(len=1) :: typesp
    character(len=19) :: liac, liot, matr, stoc, ouvert
    integer(kind=8) :: jliac, jliot, jvale, jva, jouv
    integer(kind=8) :: nbbloc
    real(kind=8) :: copmax
    integer(kind=8) :: kk1, kk2, kk1f, kk2f
    integer(kind=8) :: nbote, pivot, nbliai, lliac, ii
    integer(kind=8) :: niv, ifm
    integer(kind=8) :: bloc, dercol
    integer(kind=8) :: nnocu
    integer(kind=8), pointer :: scbl(:) => null()
    integer(kind=8), pointer :: scib(:) => null()
    integer(kind=8), pointer :: scde(:) => null()
!
! ----------------------------------------------------------------------
!
    call infniv(ifm, niv)
    call jemarq()
!
! --- LECTURE DES STRUCTURES DE DONNEES
!
    liac = resocu(1:14)//'.LIAC'
    liot = resocu(1:14)//'.LIOT'
    matr = resocu(1:14)//'.MATC'
    stoc = resocu(1:14)//'.SLCS'
!
    call jeveuo(liac, 'E', jliac)
    call jeveuo(liot, 'E', jliot)
    call jeveuo(stoc//'.SCIB', 'L', vi=scib)
    call jeveuo(stoc//'.SCBL', 'L', vi=scbl)
    call jeveuo(stoc//'.SCDE', 'L', vi=scde)
! ======================================================================
! --- INITIALISATION DES VARIABLES
! ======================================================================
    nnocu = cudisi(deficu, 'NNOCU')
    nbliai = nnocu
    typesp = 'S'
    copmax = xjvmax*1.0d-08
    pivot = 0
    nbbloc = scde(3)
    ouvert = '&&ELPIV2.TRAV'
    call wkvect(ouvert, 'V V L', nbbloc, jouv)
!
    do kk1 = spavan+1, nbliac
        do kk2 = 1, nbliac
            if (kk2 .gt. kk1) then
                kk1f = kk2
                kk2f = kk1
            else
                kk1f = kk1
                kk2f = kk2
            end if
            ii = scib(kk1f)
            dercol = scbl(ii)
            bloc = dercol*(dercol+1)/2
            if (.not. zl(jouv-1+ii)) then
                if ((ii .gt. 1) .and. (kk1f .ne. (spavan+1))) then
                    call jelibe(jexnum(matr//'.UALF', (ii-1)))
                    zl(jouv-2+ii) = .false.
                end if
                call jeveuo(jexnum(matr//'.UALF', ii), 'E', jvale)
                zl(jouv-1+ii) = .true.
            end if
!
            jva = jvale-1+(kk1f-1)*(kk1f)/2-bloc+kk2f
!
            if (abs(zr(jva)) .lt. copmax) then
                pivot = 1
            else
                pivot = 0
                goto 10
            end if
        end do
        if (pivot .eq. 1) then
!
            lliac = zi(jliac-1+kk1)
!
            zi(jliot+4*nbliai) = zi(jliot+4*nbliai)+1
            nbote = zi(jliot+4*nbliai)
            zi(jliot-1+nbote) = zi(jliac-1+kk1)
!
            call cutabl(indic, nbliac, ajliai, spliai, resocu, &
                        typesp, kk1, lliac)
            call cuimp2(ifm, lliac, typesp, 'PIV', resocu)
            goto 40
        end if
10      continue
    end do
!
40  continue
    call jedetr(ouvert)
    call jedema()
!
end subroutine
