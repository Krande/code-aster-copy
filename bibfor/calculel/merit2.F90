! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine merit2(modele, nchar, lchar, cara, time, &
                  matel, prefch, numero, base)
    implicit none
!
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
    character(len=8) :: modele, cara
    character(len=19) :: matel, prefch
    character(len=*) :: lchar(*)
    character(len=24) :: time
    character(len=1) :: base
    integer :: nchar, numero
! ----------------------------------------------------------------------
!
!     CALCUL DES MATRICES ELEMENTAIRES DE RIGIDITE THERMIQUE (2)
!        ( ISO_FACE, 'RIGI_THER_ECHA_R/F' , 'RIGI_THER_PARO_R/F' )
!
!     LES RESUELEM PRODUITS S'APPELLENT :
!           PREFCH(1:8).ME000I , I=NUMERO+1,NUMERO+N
!
!     ENTREES:
!
!     LES NOMS QUI SUIVENT SONT LES PREFIXES UTILISATEUR K8:
!        MODELE : NOM DU MODELE
!        NCHAR  : NOMBRE DE CHARGES
!        LCHAR  : LISTE DES CHARGES
!        CARA   : CHAMP DE CARAC_ELEM
!        MATEL  : NOM DU MATR_ELEM (N RESUELEM) PRODUIT
!        PREFCH : PREFIXE DES NOMS DES RESUELEM STOCKES DANS MATEL
!        NUMERO : NUMERO D'ORDRE A PARTIR DUQUEL ON NOMME LES RESUELEM
!        TIME   : CHAMPS DE TEMPSR
!
!     SORTIES:
!        MATEL  : EST REMPLI.
!
! ----------------------------------------------------------------------
!
!     FONCTIONS EXTERNES:
!     -------------------
!
!     VARIABLES LOCALES:
!     ------------------
!
!
    character(len=8) :: nomcha, lpain(4), lpaout(1), k8bid
    character(len=16) :: option
    character(len=24) :: ligrel(2), lchin(5), lchout(1), chgeom, chcara(18)
    integer :: iret, ilires, icha
! ----------------------------------------------------------------------
    integer :: nbchmx
!-----------------------------------------------------------------------
    integer :: iret3, k
!-----------------------------------------------------------------------
    parameter(nbchmx=2)
    integer :: nligr(nbchmx)
    character(len=6) :: nompar(nbchmx), nomchp(nbchmx), nomopt(nbchmx)
    data nomchp/'.COEFH', '.HECHP'/
    data nomopt/'_COEH_', '_PARO_'/
    data nompar/'PCOEFH', 'PHECHP'/
    data nligr/1, 2/
!
!
!     -- ON VERIFIE LA PRESENCE PARFOIS NECESSAIRE DE CARA_ELEM
    call jemarq()
    if (modele(1:1) .ne. ' ') then
    else
        call utmess('F', 'CALCULEL3_50')
    end if
!
    call megeom(modele, chgeom)
    call mecara(cara, chcara)
!
    call jeexin(matel//'.RERR', iret)
    if (iret .gt. 0) then
        call jedetr(matel//'.RERR')
        call jedetr(matel//'.RELR')
    end if
    call memare('V', matel, modele, ' ', cara, &
                'RIGI_THER')
!
    lpaout(1) = 'PMATTTR'
    lchout(1) = prefch(1:8)//'.ME000'
    ilires = 0
    if (lchar(1) (1:8) .ne. '        ') then
        ligrel(1) = modele(1:8)//'.MODELE'
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PTEMPSR'
        lchin(2) = time
!
        do icha = 1, nchar
            nomcha = lchar(icha)
            ligrel(2) = nomcha(1:8)//'.CHTH.LIGRE'
            call dismoi('TYPE_CHARGE', nomcha, 'CHARGE', repk=k8bid)
            if (k8bid(5:7) .eq. '_FO') then
                option = 'RIGI_THER_    _F'
                lpain(3) = '      F'
            else
                option = 'RIGI_THER_    _R'
                lpain(3) = '      R'
            end if
            do k = 1, nbchmx
                lchin(3) = nomcha//'.CHTH'//nomchp(k)//'.DESC'
                call jeexin(lchin(3), iret3)
                if (iret3 .gt. 0) then
                    option(10:15) = nomopt(k)
                    lpain(3) (1:6) = nompar(k)
                    ilires = ilires+1
                    call codent(ilires+numero, 'D0', lchout(1) (12:14))
                    call calcul('S', option, ligrel(nligr(k)), 3, lchin, &
                                lpain, 1, lchout, lpaout, base, &
                                'OUI')
                    call reajre(matel, lchout(1), base)
                end if
            end do
        end do
    end if
    call jedema()
end subroutine
