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

subroutine xtopoh(noma, modele)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/dbgcal.h"
#include "asterfort/infdbg.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/xcelno.h"
!
    character(len=8) :: modele, noma
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM - PREPARATION)
!
! AJOUTER À LA SD FISS_XFEM LES IDENTIFIANTS DES DOMAINES PAR NOEUD/FACETTE/SOUS-EL.
! C EST A DIRE LE CODAGE ASSOCIE
!   * A L INTERSECTION ENTRE LES DOMAINES DE DISCONTINUITE ET LE SUPPORT DU NOEUD X-FEM
!   * A CHAQUE DOMAINE MAITRE/ESCL.VU PAR CHAQUE FACETTE
!   * A CHAQUE  DOMAINE DE DISCONINUITE AUQUEL APPARTIENT LE SOUS-ELEMENT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   REMARQUE IMPORTANTE :
!     LA CONTIGUITE DU CODAGE PAR NOEUD EST CALCULEE DANS XCELNO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ----------------------------------------------------------------------
!
!
!  IN  MODELE : NOM DE L'OBJET MODELE
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=3, nbin=8)
    character(len=8) :: lpaout(nbout), lpain(nbin), licmp(2)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=19) :: cnseto, heavto, loncha, heavtn, ligrel
    character(len=19) ::  heavts, heavtf, heavfa, longco, lnno, fisno, fissco
    aster_logical :: debug
    integer(kind=8) :: ifm, niv, ifmdbg, nivdbg, ima, nbma
    integer(kind=8) :: jcesd, jcesl, iad
    character(len=16) :: option
    integer(kind=8), pointer :: nbsp(:) => null()
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
    option = 'TOPONO'
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    end if
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
! --- RECUPERATION DES DONNEES XFEM (TOPOSE)
!
    cnseto = modele(1:8)//'.TOPOSE.CNS'
    heavto = modele(1:8)//'.TOPOSE.HEA'
    loncha = modele(1:8)//'.TOPOSE.LON'
    heavfa = modele(1:8)//'.TOPOFAC.HE'
    longco = modele(1:8)//'.TOPOFAC.LO'
    heavtn = modele(1:8)//'.TOPONO.HNO'
    heavts = modele(1:8)//'.TOPONO.HSE'
    heavtf = modele(1:8)//'.TOPONO.HFA'
    fissco = modele(1:8)//'.FISSCO'
    lnno = modele(1:8)//'.LNNO'
    fisno = modele(1:8)//'.FISSNO'
!
! --- CREATION DES LISTES DES CHAMPS IN
!
    lpain(1) = 'PCNSETO'
    lchin(1) = cnseto
    lpain(2) = 'PHEAVTO'
    lchin(2) = heavto
    lpain(3) = 'PLONCHA'
    lchin(3) = loncha
    lpain(4) = 'PHEAVFA'
    lchin(4) = heavfa
    lpain(5) = 'PLONGCO'
    lchin(5) = longco
    lpain(6) = 'PFISNO'
    lchin(6) = fisno
    lpain(7) = 'PLEVSET'
    lchin(7) = lnno
    lpain(8) = 'PFISCO'
    lchin(8) = fissco
!
! --- CREATION DES LISTES DES CHAMPS OUT
!
    lpaout(1) = 'PHEA_NO'
    lchout(1) = heavtn
    lpaout(2) = 'PHEA_SE'
    lchout(2) = heavts
    lpaout(3) = 'PHEA_FA'
    lchout(3) = heavtf
!
! --- POUR LE MULTI-HEAVISIDE, LE CHAMP DE SORTIE PHEA_FA EST
! --- DUPLIQUÉ PAR LE NOMBRE DE FISSURES VUES
!
    call jeveuo('&&XTYELE.NBSP', 'L', vi=nbsp)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    licmp(1) = 'NPG_DYN'
    licmp(2) = 'NCMP_DYN'
    call cescre('V', heavtf, 'ELEM', noma, 'DCEL_I', &
                2, licmp, [0], [-1], [-2])
    call jeveuo(heavtf//'.CESD', 'L', jcesd)
    call jeveuo(heavtf//'.CESV', 'E', vi=cesv)
    call jeveuo(heavtf//'.CESL', 'E', jcesl)
!
! --- REMPLISSAGE DES SOUS-POINTS
    do ima = 1, nbma
        call cesexi('S', jcesd, jcesl, ima, 1, &
                    1, 1, iad)
        zl(jcesl-1-iad) = .true.
        cesv(1-1-iad) = nbsp(ima)
    end do
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
! --- CONCATENATION DU CHAMP ELNO
!        CONSTRUCTION DE LA LISTE DE DOMAINES VU PAR UN NOEU XFEM
!    call imprsd('CHAMP',modele(1:8)//'.TOPONO.HNO',29,'CHAMP HEAVISIDE PAR NOEUD')
!    call imprsd('CHAMP',modele(1:8)//'.FISSNO',28,'FISSNO')
!
    call xcelno(noma, modele, lchout(1), option, lpaout(1))
!
!    call imprsd('CHAMP',modele(1:8)//'.TOPONO.HNO',30,'CHAMP HEAVISIDE PAR NOEUD')
!    call imprsd('CHAMP',modele(1:8)//'.TOPONO.HSE',31,'CHAMP HEAVISIDE PAR SSE')
!    call imprsd('CHAMP',modele(1:8)//'.TOPONO.HFA',32,'CHAMP HEAVISIDE PAR FACETTE')
!
    call jedema()
end subroutine
