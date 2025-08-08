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
subroutine merit1(modelZ, caraElemZ, matecoZ, &
                  loadNameZ, &
                  timeMap, matrElem, resuElemPref, &
                  indxMatrElem, jvBase)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ, matecoZ
    character(len=*), intent(in) :: loadNameZ
    character(len=24), intent(in) :: timeMap
    character(len=19), intent(in) :: matrElem, resuElemPref
    integer(kind=8), intent(in) :: indxMatrElem
    character(len=1), intent(in) :: jvBase
! ----------------------------------------------------------------------
!
!     CALCUL DES MATRICES ELEMENTAIRES DE RIGIDITE THERMIQUE (1)
!            ( ISO     , 'RIGI_THER'  )
!            ( CAL_TI  , 'THER_DDLM_R')
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
!        MATE   : CHAMP DE MATERIAUX
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
    character(len=8) :: lpain(6), lpaout(1)
    character(len=16) :: option
    character(len=24) :: chgeom, chharm, lchin(6), lchout(1)
    character(len=24) :: ligrmo, ligrch, chcara(18)
    integer(kind=8) :: ilires, iret
    integer(kind=8), parameter :: nbHarm = 0

! -----------------------------------------------------------------------
    call jemarq()

    model = modelz
    caraElem = caraElemZ
    loadName = loadNameZ
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)

! - Create input fields
    call megeom(modelZ, chgeom)
    call mecara(caraElemZ, chcara)
    call meharm(modelZ, nbHarm, chharm)
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
    if (model .ne. ' ') then
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = matecoZ(1:24)
        lpain(3) = 'PCACOQU'
        lchin(3) = chcara(7)
        lpain(4) = 'PINSTR'
        lchin(4) = timeMap
        lpain(5) = 'PHARMON'
        lchin(5) = chharm
        lpain(6) = 'PCAMASS'
        lchin(6) = chcara(12)
        option = 'RIGI_THER'
        ilires = ilires+1
        call codent(ilires+indxMatrElem, 'D0', lchout(1) (12:14))
        call calcul('S', option, ligrmo, 6, lchin, &
                    lpain, 1, lchout, lpaout, jvBase, &
                    'OUI')
        call reajre(matrElem, lchout(1), jvBase)
    end if
    if (loadName .ne. " ") then
        lpain(1) = 'PDDLMUR'
        call exisd('CHAMP_GD', loadName//'.CHTH.CMULT', iret)
        if (iret .ne. 0) then
            lchin(1) = loadName//'.CHTH.CMULT     '
            ilires = ilires+1
            call codent(ilires+indxMatrElem, 'D0', lchout(1) (12:14))
            ligrch = loadName//'.CHTH.LIGRE'
            option = 'THER_DDLM_R'
            call calcul('S', option, ligrch, 1, lchin, &
                        lpain, 1, lchout, lpaout, jvBase, &
                        'OUI')
            call reajre(matrElem, lchout(1), jvBase)
        end if
    end if
!
    call jedema()
end subroutine
