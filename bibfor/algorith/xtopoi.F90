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

subroutine xtopoi(noma, modele)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/dbgcal.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/xchkgp.h"
    character(len=8) :: noma, modele
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM - PREPARATION)
!
! AJOUTER À LA SD FISS_XFEM LES DONNÉES TOPOLOGIQUES CONCERNANT
! LA DÉCOUPE DES ÉLÉMENTS POUR L'INTÉGRATION
!
! ----------------------------------------------------------------------
!
!
!  IN  MODELE : NOM DE L'OBJET MODELE
!  I/O FISS   : NOM DE LA SD FISS_XFEM
!
!
!
!
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=7, nbin=3)
    character(len=8) :: lpaout(nbout), lpain(nbin), licmp(2)
    character(len=16) :: option
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=19) :: ligrel, chgeom, painto
    character(len=19) :: pintto, cnseto, heavto, loncha, pmilto, joncno
    aster_logical :: debug
    integer(kind=8) :: ifm, niv, ifmdbg, nivdbg, ima, nbma
    integer(kind=8) :: jcesd, jcesl, iad
    integer(kind=8), pointer :: nbsp(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('XFEM', ifm, niv)
    call infdbg('PRE_CALCUL', ifmdbg, nivdbg)
!
! --- INITIALISATIONS
!
    ligrel = modele//'.MODELE'
    option = 'TOPOSE'
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    end if
    call jeveuo(modele//'.MODELE    .LGRF', 'L', vk8=lgrf)
    chgeom = lgrf(1)//'.COORDO'
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! --- RECUPERATION DES DONNEES XFEM (TOPOSE)
!
    pintto = modele(1:8)//'.TOPOSE.PIN'
    cnseto = modele(1:8)//'.TOPOSE.CNS'
    heavto = modele(1:8)//'.TOPOSE.HEA'
    loncha = modele(1:8)//'.TOPOSE.LON'
    pmilto = modele(1:8)//'.TOPOSE.PMI'
    painto = modele(1:8)//'.TOPOSE.PAI'
    joncno = modele(1:8)//'.TOPOSE.PJO'
!
    licmp(1) = 'NPG_DYN'
    licmp(2) = 'NCMP_DYN'
    call cescre('V', heavto, 'ELEM', noma, 'DCEL_I', &
                2, licmp, [0], [-1], [-2])
    call jeveuo(heavto//'.CESD', 'L', jcesd)
    call jeveuo(heavto//'.CESV', 'E', vi=cesv)
    call jeveuo(heavto//'.CESL', 'E', jcesl)
!
! --- RECUPERATION DU NOMBRE DE FISSURES VUES
!
    call jeveuo('&&XTYELE.NBSP', 'L', vi=nbsp)
!
! --- REMPLISSAGE DES SOUS POINT DE HEAVTO
!
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    do ima = 1, nbma
        call cesexi('S', jcesd, jcesl, ima, 1, &
                    1, 1, iad)
        zl(jcesl-1-iad) = .true.
        cesv(1-1-iad) = nbsp(ima)
!
    end do
!
! --- CREATION DES LISTES DES CHAMPS IN
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PLEVSET'
    lchin(2) = modele(1:8)//'.LNNO'
    lpain(3) = 'PFISCO'
    lchin(3) = modele(1:8)//'.FISSCO'
!
! --- CREATION DES LISTES DES CHAMPS OUT
!
    lpaout(1) = 'PPINTTO'
    lchout(1) = pintto
    lpaout(2) = 'PCNSETO'
    lchout(2) = cnseto
    lpaout(3) = 'PHEAVTO'
    lchout(3) = heavto
    lpaout(4) = 'PLONCHA'
    lchout(4) = loncha
    lpaout(5) = 'PPMILTO'
    lchout(5) = pmilto
    lpaout(6) = 'PAINTTO'
    lchout(6) = painto
    lpaout(7) = 'PJONCNO'
    lchout(7) = joncno
!
! --- APPEL A CALCUL
!
    call calcul('C', option, ligrel, nbin, lchin, &
                lpain, nbout, lchout, lpaout, 'G', &
                'OUI')
!
    if (debug) then
        call dbgcal(option, ifmdbg, nbin, lpain, lchin, &
                    nbout, lpaout, lchout)
    end if
!
!   vérification de la cohérence entre le nombre de points de Gauss générés
!   par la découpe en sous-éléments et le nombre maximal de points de Gauss
!   défini par la famille XFEM de l'élément parent
    call xchkgp(modele)
!
    call jedema()
end subroutine
