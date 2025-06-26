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
subroutine merit2(modelZ, caraElemZ, &
                  loadNameZ, &
                  timeMap, matrElem, resuElemPref, &
                  indxMatrElem, jvBase)
!
    implicit none
!
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
#include "jeveux.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ
    character(len=*), intent(in) :: loadNameZ
    character(len=24), intent(in) :: timeMap
    character(len=19), intent(in) :: matrElem, resuElemPref
    integer(kind=8), intent(in) :: indxMatrElem
    character(len=1), intent(in) :: jvBase
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
    character(len=8) :: model, caraElem, loadName
    character(len=8) :: lpain(4), lpaout(1), answer
    character(len=16) :: option
    character(len=24) :: ligrel(2), lchin(5), lchout(1), chgeom, chcara(18), ligrmo
    integer(kind=8) :: iret, ilires, iNeutType
    integer(kind=8), parameter :: nbNeutType = 2
    character(len=6), parameter :: loadField(nbNeutType) = (/'.COEFH', '.HECHP'/)
    character(len=6), parameter :: loadOption(nbNeutType) = (/'_COEH_', '_PARO_'/)
    character(len=6), parameter :: paraName(nbNeutType) = (/'PCOEFH', 'PHECHP'/)
    integer(kind=8), parameter :: ligrelType(nbNeutType) = (/1, 2/)
!
!
!     -- ON VERIFIE LA PRESENCE PARFOIS NECESSAIRE DE CARA_ELEM
    call jemarq()

    model = modelz
    caraElem = caraElemZ
    loadName = loadNameZ
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)

    if (model .eq. ' ') then
        call utmess('F', 'CALCULEL3_50')
    end if
!
    call megeom(modelZ, chgeom)
    call mecara(caraElemZ, chcara)
!
    call jeexin(matrElem//'.RERR', iret)
    if (iret .gt. 0) then
        call jedetr(matrElem//'.RERR')
        call jedetr(matrElem//'.RELR')
    end if
    call memare('V', matrElem, modelZ, 'RIGI_THER')
!
    lpaout(1) = 'PMATTTR'
    lchout(1) = resuElemPref(1:8)//'.ME000'
    ilires = 0
    if (loadName .ne. " ") then
        ligrel(1) = ligrmo
        ligrel(2) = loadName(1:8)//'.CHTH.LIGRE'
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PINSTR'
        lchin(2) = timeMap
        call dismoi('TYPE_CHARGE', loadName, 'CHARGE', repk=answer)
        if (answer(5:7) .eq. '_FO') then
            option = 'RIGI_THER_    _F'
            lpain(3) = '      F'
        else
            option = 'RIGI_THER_    _R'
            lpain(3) = '      R'
        end if
        do iNeutType = 1, nbNeutType
            lchin(3) = loadName//'.CHTH'//loadField(iNeutType)//'.DESC'
            call jeexin(lchin(3), iret)
            if (iret .gt. 0) then
                option(10:15) = loadOption(iNeutType)
                lpain(3) (1:6) = paraName(iNeutType)
                ilires = ilires+1
                call codent(ilires+indxMatrElem, 'D0', lchout(1) (12:14))
                call calcul('S', option, ligrel(ligrelType(iNeutType)), 3, lchin, &
                            lpain, 1, lchout, lpaout, jvBase, &
                            'OUI')
                call reajre(matrElem, lchout(1), jvBase)
            end if
        end do
    end if

    call jedema()
end subroutine
