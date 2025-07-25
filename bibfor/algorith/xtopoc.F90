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

subroutine xtopoc(modele, decou)
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
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/xchdec.h"
    character(len=8) :: modele, decou
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (METHODE XFEM - PREPARATION)
!
! AJOUTER À LA SD FISS_XFEM LES DONNÉES TOPOLOGIQUES CONCERNANT
! LES FACETTES DE CONTACT
!
! ----------------------------------------------------------------------
!
!
!  IN  MODELE : NOM DE L'OBJET MODELE
!
!
!
!
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=7, nbin=15)
    character(len=8) :: lpaout(nbout), lpain(nbin), noma, licmp(2)
    character(len=8) :: nomfis, cpar
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=19) :: ligrel, chgeom
    character(len=19) :: lnno, grlnno, ltno, grltno, fissco, champ(7)
    character(len=19) :: pint, cnset, nit, phe, aint, milt, stano
    aster_logical :: debug
    character(len=16) :: option
    integer(kind=8) :: ifmdbg, nivdbg
    integer(kind=8) :: jcesd, jcesl, iad, i, nbma, ima
    integer(kind=8) :: jmofis, jnfiss, nfiss, ifiss, ityp
    character(len=16) :: typdis, memtyp
    character(len=19) :: typenr, chdec
    integer(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: nbsp(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('PRE_CALCUL', ifmdbg, nivdbg)
!
! --- INITIALISATIONS
!
    ligrel = modele//'.MODELE'
    call jeveuo(modele//'.MODELE    .LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)
    chgeom = noma//'.COORDO'
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    end if
    option = 'TOPOFA'
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout, &
                lchout)
!
    chdec = '&&XTOPOC.CHDEC'
    call xchdec(modele, decou, chdec)
!
! --- RECUPERATION DES DONNEES XFEM
!
    lnno = modele(1:8)//'.LNNO'
    ltno = modele(1:8)//'.LTNO'
    grlnno = modele(1:8)//'.GRLNNO'
    grltno = modele(1:8)//'.GRLTNO'
    fissco = modele(1:8)//'.FISSCO'
    pint = modele(1:8)//'.TOPOSE.PIN'
    cnset = modele(1:8)//'.TOPOSE.CNS'
    nit = modele(1:8)//'.TOPOSE.LON'
    phe = modele(1:8)//'.TOPOSE.HEA'
    aint = modele(1:8)//'.TOPOSE.PAI'
    milt = modele(1:8)//'.TOPOSE.PMI'
    stano = modele(1:8)//'.STNO'
    champ(1) = modele(1:8)//'.TOPOFAC.PI'
    champ(2) = modele(1:8)//'.TOPOFAC.AI'
    champ(3) = modele(1:8)//'.TOPOFAC.CF'
    champ(4) = modele(1:8)//'.TOPOFAC.LO'
    champ(5) = modele(1:8)//'.TOPOFAC.BA'
    champ(6) = modele(1:8)//'.TOPOFAC.OE'
    champ(7) = modele(1:8)//'.TOPOFAC.HE'
!
! --- CREATION D UNE CARTE POUR INFO TYPE DE DISCONTINUITE
!
    call jeveuo(modele//'.FISS', 'L', jmofis)
    call jeveuo(modele//'.NFIS', 'L', jnfiss)
    nfiss = zi(jnfiss)
    do ifiss = 1, nfiss
        nomfis = zk8(jmofis-1+ifiss)
        if (ifiss .gt. 1) memtyp = typdis
        call dismoi('TYPE_DISCONTINUITE', nomfis, 'FISS_XFEM', repk=typdis)
!
!       SI TYPE DE DISCONTINUITE DIFFERENTES, CELA DOIT ETRE DETECTE
!       EN AMONT DANS DEFI_FISS_XFEM
        if (ifiss .gt. 1) then
            if (memtyp .eq. 'FISSURE') typdis = 'FISSURE'
        end if
    end do
    typenr = '&&XTOPOC.TYPENR'
    if (typdis .eq. 'FISSURE') ityp = 1
    if (typdis .eq. 'INTERFACE') ityp = 2
    if (typdis .eq. 'COHESIF') ityp = 3
    cpar = 'X1'
    call mecact('V', typenr, 'MODELE', ligrel, 'NEUT_I', &
                ncmp=1, nomcmp=cpar, si=ityp)
!
! --- POUR LE MULTI-HEAVISIDE, TOUS LES CHAMPS DE SORTIE SONT
! --- DUPLIQUÉS PAR LE NOMBRE DE FISSURES VUES
!
    call jeveuo('&&XTYELE.NBSP', 'L', vi=nbsp)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
    licmp(1) = 'NPG_DYN'
    licmp(2) = 'NCMP_DYN'
!
    do i = 1, 7
        call cescre('V', champ(i), 'ELEM', noma, 'DCEL_I', &
                    2, licmp, [0], [-1], [-2])
        call jeveuo(champ(i)//'.CESD', 'L', jcesd)
        call jeveuo(champ(i)//'.CESV', 'E', vi=cesv)
        call jeveuo(champ(i)//'.CESL', 'E', jcesl)
!
! --- REMPLISSAGE DES SOUS-POINTS DE CHAMP(I)
!
        do ima = 1, nbma
            call cesexi('S', jcesd, jcesl, ima, 1, &
                        1, 1, iad)
            zl(jcesl-1-iad) = .true.
            cesv(1-1-iad) = nbsp(ima)
            if (i .eq. 7) cesv(1-1-iad) = nbsp(ima)**2
        end do
!
    end do
!
! --- CREATION DES LISTES DES CHAMPS IN
!
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PLSN'
    lchin(2) = lnno
    lpain(3) = 'PLST'
    lchin(3) = ltno
    lpain(4) = 'PGRADLN'
    lchin(4) = grlnno
    lpain(5) = 'PGRADLT'
    lchin(5) = grltno
    lpain(6) = 'PFISCO'
    lchin(6) = fissco
    lpain(7) = 'PPINTTO'
    lchin(7) = pint
    lpain(8) = 'PCNSETO'
    lchin(8) = cnset
    lpain(9) = 'PLONCHA'
    lchin(9) = nit
    lpain(10) = 'PHEAVTO'
    lchin(10) = phe
    lpain(11) = 'PAINTTO'
    lchin(11) = aint
    lpain(12) = 'PPMILTO'
    lchin(12) = milt
    lpain(13) = 'PTYPDIS'
    lchin(13) = typenr
    lpain(14) = 'PSTANO'
    lchin(14) = stano
    lpain(15) = 'PDECOU'
    lchin(15) = chdec
!
! --- CREATION DES LISTES DES CHAMPS OUT
!
    lpaout(1) = 'PPINTER'
    lchout(1) = champ(1)
    lpaout(2) = 'PAINTER'
    lchout(2) = champ(2)
    lpaout(3) = 'PCFACE'
    lchout(3) = champ(3)
    lpaout(4) = 'PLONGCO'
    lchout(4) = champ(4)
    lpaout(5) = 'PBASECO'
    lchout(5) = champ(5)
    lpaout(6) = 'PGESCLA'
    lchout(6) = champ(6)
    lpaout(7) = 'PHEAVFA'
    lchout(7) = champ(7)
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
    call detrsd('CHAM_ELEM', chdec)
!
    call jedema()
end subroutine
