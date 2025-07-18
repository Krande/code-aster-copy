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

subroutine merige(model_, cara_elem_, sigg, strx, matel, &
                  base, nh, deplr, mateco)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/exixfe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/reajre.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nh
    character(len=1) :: base
    character(len=*) :: sigg, strx
    character(len=19) :: matel
    character(len=*), intent(in) :: model_
    character(len=*), intent(in) :: cara_elem_
    character(len=*), optional, intent(in) :: deplr
    character(len=*), optional, intent(in) :: mateco
!
!     CALCUL DES MATRICES ELEMENTAIRES DE RIGIDITE GEOMETRIQUE
!
!     ------------------------------------------------------------------
! IN  : MODELE : NOM DU MODELE
! IN  : CARA   : CHAMP DE CARAC_ELEM
! IN  : SIGG   : CHAMP DE CONTRAINTES AUX POINTS DE GAUSS
! IN  : NH     : NUMERO DE L'HARMONIQUE DE FOURIER
! VAR : MATEL  : NOM DU MATEL (N RESUELEM) PRODUIT
! IN  : BASE   : BASE POUR LA CREATION DE MATEL ('G'/'V')
! ----------------------------------------------------------------------
    character(len=8) ::  lpain(14), lpaout(1)
    character(len=24) :: lchin(14), lchout(1)
!
    character(len=16) :: option
    character(len=24) :: ligrmo, chgeom, chcara(18), chharm
    character(len=19) :: pintto, cnseto, heavto, loncha, basloc, lsn, lst, stano, pmilto, hea_no
    character(len=8) :: modele, cara
!
!-----------------------------------------------------------------------
    integer(kind=8) :: icode, ier, nbpara
!-----------------------------------------------------------------------
    call jemarq()
!
    modele = model_
    cara = cara_elem_
    option = 'RIGI_GEOM'
    if (modele(1:1) .eq. ' ') then
        call utmess('F', 'CALCULEL2_82')
    end if
    call detrsd('MATR_ELEM', matel)
    call mecham(option, modele, cara, nh, chgeom, &
                chcara, chharm, icode)
!
    call memare(base, matel, modele, option)
!
!  -----CAS DU MODELE X-FEM-----------------------
    call exixfe(modele, ier)
    if (ier .ne. 0) then
!
        pintto = modele(1:8)//'.TOPOSE.PIN'
        cnseto = modele(1:8)//'.TOPOSE.CNS'
        heavto = modele(1:8)//'.TOPOSE.HEA'
        loncha = modele(1:8)//'.TOPOSE.LON'
        pmilto = modele(1:8)//'.TOPOSE.PMI'
        hea_no = modele(1:8)//'.TOPONO.HNO'
        basloc = modele(1:8)//'.BASLOC'
        lsn = modele(1:8)//'.LNNO'
        lst = modele(1:8)//'.LTNO'
        stano = modele(1:8)//'.STNO'
!
        ligrmo = modele//'.MODELE'
!
! ----- REMPLISSAGE DES CHAMPS D'ENTREE
!
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PCONTRR'
        lchin(2) = sigg
        lpain(3) = 'PPINTTO'
        lchin(3) = pintto
        lpain(4) = 'PHEAVTO'
        lchin(4) = heavto
        lpain(5) = 'PLONCHA'
        lchin(5) = loncha
        lpain(6) = 'PCNSETO'
        lchin(6) = cnseto
        lpain(7) = 'PBASLOR'
        lchin(7) = basloc
        lpain(8) = 'PLSN'
        lchin(8) = lsn
        lpain(9) = 'PLST'
        lchin(9) = lst
        lpain(10) = 'PSTANO'
        lchin(10) = stano
        lpain(11) = 'PPMILTO'
        lchin(11) = pmilto
        lpain(12) = 'PSTRXRR'
        lchin(12) = strx
        lpain(13) = 'PHEA_NO'
        lchin(13) = hea_no
        nbpara = 13
!
! --- CHAMPS DE SORTIE
!
        lpaout(1) = 'PMATUUR'
        lchout(1) = matel(1:15)//'.ME001'
!
        option = 'RIGI_GEOM'
!
        call calcul('S', option, ligrmo, nbpara, lchin, &
                    lpain, 1, lchout, lpaout, base, &
                    'OUI')
        call reajre(matel, lchout(1), base)
!
    else if (ier .eq. 0) then
!
        lpaout(1) = 'PMATUUR'
        lchout(1) = matel(1:8)//'.ME001'
!
        ligrmo = modele//'.MODELE'
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PCONTRR'
        lchin(2) = sigg
        lpain(3) = 'PCAORIE'
        lchin(3) = chcara(1)
        lpain(4) = 'PCADISK'
        lchin(4) = chcara(2)
        lpain(5) = 'PCAGNPO'
        lchin(5) = chcara(6)
        lpain(6) = 'PCACOQU'
        lchin(6) = chcara(7)
        lpain(7) = 'PEFFORR'
        lchin(7) = sigg
        lpain(8) = 'PHARMON'
        lchin(8) = chharm
        lpain(9) = 'PNBSP_I'
        lchin(9) = chcara(16)
        lpain(10) = 'PSTRXRR'
        lchin(10) = strx
        lpain(11) = 'PFIBRES'
        lchin(11) = chcara(17)
        lpain(12) = 'PCACABL'
        lchin(12) = chcara(10)
        nbpara = 12
        if (present(deplr)) then
            if (deplr .ne. ' ') then
                nbpara = nbpara+1
                lpain(nbpara) = 'PDEPLPR'
                lchin(nbpara) = deplr
            end if
        end if
        if (present(mateco)) then
            if (mateco .ne. ' ') then
                nbpara = nbpara+1
                lpain(nbpara) = 'PMATERC'
                lchin(nbpara) = mateco
            end if
        end if

        option = 'RIGI_GEOM'
        call calcul('S', option, ligrmo, nbpara, lchin, &
                    lpain, 1, lchout, lpaout, base, &
                    'OUI')
        call reajre(matel, lchout(1), base)
!
    end if
!
    call jedema()
end subroutine
