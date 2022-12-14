! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine xgecfi(modele, depgeo)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/copisd.h"
#include "asterfort/dbgcal.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/inical.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/megeom.h"
#include "asterfort/getvid.h"
#include "asterfort/rcmfmc.h"
    character(len=8) :: modele
    character(len=19) :: depgeo
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM-GG
!
! CREATION DES SD CONTENANT LA GEOMETRIE DES
!                FACETTES DE CONTACT (ESCLAVES ET MAITRE)
!
! ----------------------------------------------------------------------
!
!
! ----------------------------------------------------------------------
! ROUTINE SPECIFIQUE A L'APPROCHE <<GRANDS GLISSEMENTS AVEC XFEM>>,
! TRAVAIL EFFECTUE EN COLLABORATION AVEC I.F.P.
! ----------------------------------------------------------------------
!
! IN  MODELE : NOM DU MODELE
! IN  DEPGEO : CHAMP DE DEPLACEMENTS
!
!
!
!
    integer :: nbout, nbin
    parameter    (nbout=2, nbin=13)
    character(len=8) :: lpaout(nbout), lpain(nbin), licmp(2), noma, materi
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    character(len=16) :: option
    character(len=19) :: ligrel, pinter, faclon, newges, newgem
    character(len=19) :: gesclo, ltno, fissno, hea_no, hea_fa
    character(len=19) :: chgeom, basloc, lnno, stano, mate
    character(len=1) :: base
    aster_logical :: debug
    integer :: ifmdbg, nivdbg, iret, nbma, ima
    integer :: jcesd, jcesl, iad, n1
    integer, pointer :: cesv(:) => null()
    character(len=8), pointer :: lgrf(:) => null()
    integer, pointer :: cesd2(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('PRE_CALCUL', ifmdbg, nivdbg)
!
! --- INITIALISATIONS
!
    if (nivdbg .ge. 2) then
        debug = .true.
    else
        debug = .false.
    endif
    base = 'V'
    option = 'GEOM_FAC'
!
! --- INITIALISATION DES CHAMPS POUR CALCUL
!
    call inical(nbin, lpain, lchin, nbout, lpaout,&
                lchout)
    ligrel = modele//'.MODELE'
    pinter = modele//'.TOPOFAC.PI'
    faclon = modele//'.TOPOFAC.LO'
    gesclo = modele//'.TOPOFAC.OE'
    ltno = modele//'.LTNO'
    fissno = modele(1:8)//'.FISSNO'
!    heavfa = modele(1:8)//'.TOPOFAC.HE'
    newges = modele//'.TOPOFAC.GE'
    newgem = modele//'.TOPOFAC.GM'
    hea_no = modele//'.TOPONO.HNO'
    hea_fa = modele//'.TOPONO.HFA'
!  AJOUT DES CHAMPS NECESSAIRES POUR L ASSEMBLAGE DU CHAMP SINGULIER
    call megeom(modele, chgeom)
    lnno = modele(1:8)//'.LNNO'
    stano = modele(1:8)//'.STNO'
    basloc = modele(1:8)//'.BASLOC'
    materi = ' '
    call getvid(' ', 'CHAM_MATER', scal=materi, nbret=n1)
    if (n1 .ne. 0) then
       call rcmfmc(materi, mate, l_ther_ = ASTER_FALSE)
    else
       mate = ' '
    endif

!
    call jeexin(newges//'.CESD', iret)
    if (iret .eq. 0) then
!
! --- RECOPIE DU NOMBRE DE SOUS POINTS DE TOPOFAC.OE DANS TOPOFAC.GE/M
!
        call jeveuo(modele//'.MODELE    .LGRF', 'L', vk8=lgrf)
        noma = lgrf(1)
        call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
        call celces(gesclo, 'V', '&&XGECFI.GESCLO')
        call jeveuo('&&XGECFI.GESCLO    .CESD', 'L', vi=cesd2)
        licmp(1) = 'NPG_DYN'
        licmp(2) = 'NCMP_DYN'
        call cescre('V', newges, 'ELEM', noma, 'DCEL_I',&
                    2, licmp, [0], [-1], [-2])
        call jeveuo(newges//'.CESD', 'L', jcesd)
        call jeveuo(newges//'.CESV', 'E', vi=cesv)
        call jeveuo(newges//'.CESL', 'E', jcesl)
        do ima = 1, nbma
            call cesexi('S', jcesd, jcesl, ima, 1,&
                        1, 1, iad)
            zl(jcesl-1-iad) = .true.
            cesv(1-1-iad) = cesd2(5+4*(ima-1)+2)
        end do
        call copisd('CHAM_ELEM_S', 'V', newges, newgem)
        call detrsd('CHAM_ELEM_S', '&&XGECFI.GESCLO')
    endif
!
! --- CREATION DES LISTES DES CHAMPS IN ET OUT
!
    lpain(1) = 'PDEPLA'
    lchin(1) = depgeo
    lpain(2) = 'PPINTER'
    lchin(2) = pinter
    lpain(3) = 'PLONGCO'
    lchin(3) = faclon
    lpain(4) = 'PGESCLO'
    lchin(4) = gesclo
    lpain(5) = 'PLST'
    lchin(5) = ltno
    lpain(6) = 'PFISNO'
    lchin(6) = fissno
    lpain(7) = 'PHEA_FA'
    lchin(7) = hea_fa
    lpain(8) = 'PHEA_NO'
    lchin(8) = hea_no
    lpain(9) = 'PGEOMER'
    lchin(9) = chgeom
    lpain(10) = 'PBASLOR'
    lchin(10) = basloc
    lpain(11) = 'PSTANO'
    lchin(11) = stano
    lpain(12) = 'PLSN'
    lchin(12) = lnno
    lpain(13) = 'PMATERC'
    lchin(13) = mate
!
    lpaout(1) = 'PNEWGES'
    lchout(1) = newges
    lpaout(2) = 'PNEWGEM'
    lchout(2) = newgem
!
! --- APPEL A CALCUL
!
    call calcul('C', option, ligrel, nbin, lchin,&
                lpain, nbout, lchout, lpaout, base,&
                'OUI')
    if (debug) then
        call dbgcal(option, ifmdbg, nbin, lpain, lchin,&
                    nbout, lpaout, lchout)
    endif
!
    call jedema()
end subroutine
