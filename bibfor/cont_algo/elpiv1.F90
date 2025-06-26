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

subroutine elpiv1(xjvmax, indic, nbliac, ajliai, spliai, &
                  spavan, noma, sdcont_defi, sdcont_solv)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cfdisd.h"
#include "asterfort/cfimp2.h"
#include "asterfort/cftabl.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    character(len=8) :: noma
    character(len=24) :: sdcont_solv, sdcont_defi
    real(kind=8) :: xjvmax
    integer(kind=8) :: nbliac
    integer(kind=8) :: indic
    integer(kind=8) :: ajliai, spliai
    integer(kind=8) :: spavan
!
! ----------------------------------------------------------------------
!
! ROUTINE CONTACT (METHODES DISCRETES - UTILITAIRE)
!
! ELIMINATION DES PIVOTS NULS DANS LA MATRICE DE CONTACT
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
! IN  NOMA   : NOM DU MAILLAGE
! IN  DEFICO : SD DE DEFINITION DU CONTACT (ISSUE D'AFFE_CHAR_MECA)
! IN  RESOCO : SD DE TRAITEMENT NUMERIQUE DU CONTACT
!                'E': RESOCO(1:14)//'.LIAC'
!                'E': RESOCO(1:14)//'.LIOT'
!
!
!
!
!
    character(len=1) :: typesp
    character(len=19) :: liac, liot, macont, stoc, ouvert
    integer(kind=8) :: jliac, jliot, jvale, jva, jouv
    integer(kind=8) :: nbbloc, nbliai
    real(kind=8) :: copmax
    integer(kind=8) :: kk1, kk2, kk1f, kk2f
    integer(kind=8) :: nbote, lliac
    integer(kind=8) :: iblc
    integer(kind=8) :: niv, ifm
    integer(kind=8) :: bloc, dercol
    aster_logical :: pivnul
    integer(kind=8), pointer :: scbl(:) => null()
    integer(kind=8), pointer :: scde(:) => null()
    integer(kind=8), pointer :: scib(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('CONTACT', ifm, niv)
!
! --- LECTURE DES STRUCTURES DE DONNEES
!
    liac = sdcont_solv(1:14)//'.LIAC'
    liot = sdcont_solv(1:14)//'.LIOT'
    macont = sdcont_solv(1:14)//'.MATC'
    stoc = sdcont_solv(1:14)//'.SLCS'
    call jeveuo(liac, 'E', jliac)
    call jeveuo(liot, 'E', jliot)
    call jeveuo(stoc//'.SCIB', 'L', vi=scib)
    call jeveuo(stoc//'.SCBL', 'L', vi=scbl)
    call jeveuo(stoc//'.SCDE', 'L', vi=scde)
!
! --- INITIALISATIONS
!
    nbliai = cfdisd(sdcont_solv, 'NBLIAI')
    typesp = 'S'
    copmax = xjvmax*1.0d-08
    pivnul = .false.
    nbbloc = scde(3)
!
! --- BLOC EN LECTURE
!
    ouvert = '&&ELPIV2.TRAV'
    call wkvect(ouvert, 'V V L', nbbloc, jouv)
!
! --- DETECTION DES PIVOTS NULS
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
            iblc = scib(kk1f)
            dercol = scbl(iblc)
            bloc = dercol*(dercol+1)/2
!
! ------- ON ACCEDE AU BLOC
!
            if (.not. zl(jouv-1+iblc)) then
                if ((iblc .gt. 1) .and. (kk1f .ne. (spavan+1))) then
                    call jelibe(jexnum(macont//'.UALF', (iblc-1)))
                    zl(jouv-2+iblc) = .false.
                end if
                call jeveuo(jexnum(macont//'.UALF', iblc), 'E', jvale)
                zl(jouv-1+iblc) = .true.
            end if
!
! ------- ACCES A LA DIAGONALE
!
            jva = jvale-1+(kk1f-1)*(kk1f)/2-bloc+kk2f
!
! ------- PIVOT NUL ?
!
            if (abs(zr(jva)) .lt. copmax) then
                pivnul = .true.
            else
                pivnul = .false.
                goto 10
            end if
        end do
!
! ----- ON SUPPRIME LA LIAISON
!
        if (pivnul) then
            lliac = zi(jliac-1+kk1)
            zi(jliot+4*nbliai) = zi(jliot+4*nbliai)+1
            nbote = zi(jliot+4*nbliai)
            zi(jliot-1+nbote) = zi(jliac-1+kk1)
            call cftabl(indic, nbliac, ajliai, spliai, &
                        sdcont_solv, typesp, kk1, &
                        lliac)
            call cfimp2(sdcont_defi, sdcont_solv, noma, lliac, &
                        'PIV')
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
