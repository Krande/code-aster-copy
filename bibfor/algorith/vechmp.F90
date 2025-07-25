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

subroutine vechmp(nomo, mate, mateco, carele, varplu, lxfem, &
                  partps, nbin_maxi, lpain, lchin, lastin)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/mecoor.h"
    integer(kind=8) :: nbin_maxi, lastin
    character(len=8) :: lpain(nbin_maxi)
    character(len=19) :: lchin(nbin_maxi)
    character(len=8) :: nomo
    aster_logical :: lxfem
    real(kind=8) :: partps(3)
    character(len=19) :: varplu
    character(len=24) :: mate, carele, mateco
!
! ----------------------------------------------------------------------
!
! CALCUL DES VECTEURS ELEMENTAIRES DES CHARGEMENTS MECANIQUES
! DE NEUMANN
!
! PREPARATION DES CHAMPS D'ENTREE STANDARDS
!
! ----------------------------------------------------------------------
!
!
! IN  NOMO   : NOM DU MODELE
! IN  PARTPS : TABLEAU DONNANT T+, DELTAT ET THETA (POUR LE THM)
! IN  CARELE : CARACTERISTIQUES DES POUTRES ET COQUES
! IN  MATE   : MATERIAU CODE
! IN  VARPLU : VARIABLES DE COMMANDE A L'INSTANT T+
! IN  LXFEM  : .TRUE. SI XFEM
! IN  nbin_maxi   : NOMBRE MAXI DE CHAMPS D'ENTREE
! OUT LPAIN  : LISTE DES PARAMETRES IN
! OUT LCHIN  : LISTE DES CHAMPS IN
! OUT LASTIN : NOMBRE EFFECTIF DE CHAMPS IN
!
! ----------------------------------------------------------------------
!
    character(len=8) :: nomcmp(3)
    character(len=19) :: ligrmo
    character(len=19) :: chgeom, chcara(18), chtime
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ligrmo = nomo(1:8)//'.MODELE'
!
! --- CHAMP DE GEOMETRIE
!
    call mecoor(nomo, chgeom)
!
! --- CHAMP DE CARACTERISTIQUES ELEMENTAIRES
!
    call mecara(carele, chcara)
!
! --- CREATION DE LA CARTE DES INSTANTS
!
    chtime = '&&VECHMP.CH_INST_R'
    nomcmp(1) = 'INST'
    nomcmp(2) = 'DELTAT'
    nomcmp(3) = 'THETA'
    call mecact('V', chtime, 'LIGREL', ligrmo, 'INST_R', &
                ncmp=3, lnomcmp=nomcmp, vr=partps)
!
! --- CHAMPS D'ENTREES STANDARDS
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PINSTR'
    lchin(2) = chtime
    lpain(3) = 'PMATERC'
    lchin(3) = mateco(1:19)
    lpain(4) = 'PCACOQU'
    lchin(4) = chcara(7)
    lpain(5) = 'PCAGNPO'
    lchin(5) = chcara(6)
    lpain(6) = 'PCADISM'
    lchin(6) = chcara(3)
    lpain(7) = 'PCAORIE'
    lchin(7) = chcara(1)
    lpain(8) = 'PCACABL'
    lchin(8) = chcara(10)
    lpain(9) = 'PCAARPO'
    lchin(9) = chcara(9)
    lpain(10) = 'PCAGNBA'
    lchin(10) = chcara(11)
    lpain(11) = 'PVARCPR'
    lchin(11) = varplu
    lpain(12) = 'PCAMASS'
    lchin(12) = chcara(12)
    lpain(13) = 'PCAGEPO'
    lchin(13) = chcara(5)
    lpain(14) = 'PNBSP_I'
    lchin(14) = chcara(16)
    lpain(15) = 'PFIBRES'
    lchin(15) = chcara(17)
    lpain(16) = 'PCINFDI'
    lchin(16) = chcara(15)
    lastin = 16
!
! --- CHAMPS NECESSAIRES POUR XFEM
!
    if (lxfem) then
        lpain(17) = 'PPINTTO'
        lchin(17) = nomo(1:8)//'.TOPOSE.PIN'
        lpain(18) = 'PCNSETO'
        lchin(18) = nomo(1:8)//'.TOPOSE.CNS'
        lpain(19) = 'PHEAVTO'
        lchin(19) = nomo(1:8)//'.TOPOSE.HEA'
        lpain(20) = 'PLONCHA'
        lchin(20) = nomo(1:8)//'.TOPOSE.LON'
        lpain(21) = 'PLSN'
        lchin(21) = nomo(1:8)//'.LNNO'
        lpain(22) = 'PLST'
        lchin(22) = nomo(1:8)//'.LTNO'
        lpain(23) = 'PSTANO'
        lchin(23) = nomo(1:8)//'.STNO'
        lpain(24) = 'PPMILTO'
        lchin(24) = nomo(1:8)//'.TOPOSE.PMI'
        lpain(25) = 'PFISNO'
        lchin(25) = nomo(1:8)//'.FISSNO'
        lpain(26) = 'PPINTER'
        lchin(26) = nomo(1:8)//'.TOPOFAC.OE'
        lpain(27) = 'PAINTER'
        lchin(27) = nomo(1:8)//'.TOPOFAC.AI'
        lpain(28) = 'PCFACE'
        lchin(28) = nomo(1:8)//'.TOPOFAC.CF'
        lpain(29) = 'PLONGCO'
        lchin(29) = nomo(1:8)//'.TOPOFAC.LO'
        lpain(30) = 'PBASECO'
        lchin(30) = nomo(1:8)//'.TOPOFAC.BA'
        lastin = 30
    end if
!
! --- PCOMPOR UTILE POUR POUTRES MULTI-FIBRES
!
    lastin = lastin+1
    lpain(lastin) = 'PCOMPOR'
    lchin(lastin) = mate(1:8)//'.COMPOR'
!
    ASSERT(lastin .le. nbin_maxi)
!
    call jedema()
end subroutine
