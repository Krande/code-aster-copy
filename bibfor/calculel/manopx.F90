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

subroutine manopx(model, ligrel, option, param, chsgeo, exixfm, &
                  kecono)
! person_in_charge: samuel.geniaut at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/indk32.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/modat2.h"
#include "asterfort/nucalc.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8), intent(in) :: model
    character(len=19) :: ligrel, chsgeo
    character(len=16) :: option
    character(len=8) :: param
    character(len=3) :: exixfm
    character(len=24) :: kecono
! ------------------------------------------------------------------
! BUT: REGARDER DANS LE LIGREL S'IL Y A DES ELEMENTS XFEM
!      ET SI LE CHAMP ELGA (OPTION/PARAM) UTILISE UNE FAMILLE XFEM..
!
!      * CALCULER UN OBJET (KECONO) DISANT SI CHAQUE GREL
!          PEUT ETRE STOCKE "ECONOMIQUE"
!      SI FAMILLE XFEM :
!        * CALCULER LE CHAMP SIMPLE CHSGEO CONTENANT LES COORDONNNEES
!          DES POINTS DE GAUSS DES FAMILLES XFEM...
! ------------------------------------------------------------------
!     ARGUMENTS:
!     ----------
! LIGREL  IN/JXIN  K19 : LIGREL
! OPTION,PARAM  IN  K* : OPTION ET PARAMETRE PERMETTANT DE DETERMINER
!                        LA FAMILLE DE PG UTILISEE.
! EXIXFM  OUT K3 : 'OUI' : IL EXISTE DES GRELS AVEC FAMILLE XFEM...
!                  'NON' SINON
! KECONO  IN/JXOUT K24 : VECTEUR D'ENTIERS LONG=NBGREL(LIGREL)
!      V(IGR) = 1 : LE GREL IGR UTILISE UNE FAMILLE XFEM...
!             = 0 SINON
! CHSGEO  IN/JXOUT K19 : CHAM_ELEM_S (GEOM_R) DE TYPE 'ELGA'
!
! REMARQUE :
!   L'OBJET CHSGEO N'EST CREE QUE SI EXIXFM='OUI'
! ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbout, nbin
    parameter(nbout=1, nbin=6)
    character(len=8) :: lpaout(nbout), lpain(nbin)
    character(len=19) :: lchout(nbout), lchin(nbin)
!
    integer(kind=8) :: iopt, iopt1, nute, numc, igr, nbgrel
    integer(kind=8) :: jecono, imolo, jmolo, nec, kfpg
    integer(kind=8) :: igd, nblfpg, nbfam, jfpgl
    integer(kind=8) :: k, nuflpg, nufgpg
    character(len=8) :: nomgd, elrefe, ma
    character(len=16) :: nofpg, nomte
    character(len=19) :: chgeom
    character(len=32) :: noflpg
    character(len=32), pointer :: pnlocfpg(:) => null()
    integer(kind=8), pointer :: nolocfpg(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=ma)
    call jelira(ligrel//'.LIEL', 'NMAXOC', nbgrel)
!
    call jeveuo('&CATA.TE.PNLOCFPG', 'L', vk32=pnlocfpg)
    call jelira('&CATA.TE.NOLOCFPG', 'LONMAX', nblfpg)
    call jeveuo('&CATA.TE.NOLOCFPG', 'L', vi=nolocfpg)
!
!
!     1. CALCUL DE KECONO ET EXIXFM :
!     ------------------------------------------------------------------
    exixfm = 'NON'
    call wkvect(kecono, 'V V I', nbgrel, jecono)
    do igr = 1, nbgrel
        zi(jecono-1+igr) = 1
    end do
    call jenonu(jexnom('&CATA.OP.NOMOPT', 'XFEM_XPG'), iopt1)
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), iopt)
    do igr = 1, nbgrel
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
!
!       L'ELEMENT SAIT-IL CALCULER XFEM_XPG ?
        numc = nucalc(iopt1, nute, 1)
        if (numc .lt. 0) goto 2
!
        imolo = modat2(iopt, nute, param)
        if (imolo .eq. 0) goto 2
!
        call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', jmolo)
        igd = zi(jmolo-1+2)
        call jenuno(jexnum('&CATA.GD.NOMGD', igd), nomgd)
        call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
        kfpg = zi(jmolo-1+4+nec+1)
!
!       -- FAMILLE "LISTE"
        if (kfpg .lt. 0) then
!          FAMILLE "LISTE" :
            call jelira(jexnum('&CATA.TE.FPG_LISTE', -kfpg), 'LONMAX', nbfam)
            nbfam = nbfam-1
            call jeveuo(jexnum('&CATA.TE.FPG_LISTE', -kfpg), 'L', jfpgl)
            elrefe = zk8(jfpgl-1+nbfam+1)
            do k = 1, nbfam
                noflpg = nomte//elrefe//zk8(jfpgl-1+k)
                nuflpg = indk32(pnlocfpg, noflpg, 1, nblfpg)
                nufgpg = nolocfpg(nuflpg)
                call jenuno(jexnum('&CATA.TM.NOFPG', nufgpg), nofpg)
                if (nofpg(9:12) .eq. 'XFEM') then
                    exixfm = 'OUI'
                    zi(jecono-1+igr) = 0
                end if
            end do
!
!       -- FAMILLE "ORDINAIRE"
        else
            call jenuno(jexnum('&CATA.TM.NOFPG', kfpg), nofpg)
            if (nofpg(9:12) .eq. 'XFEM') then
                exixfm = 'OUI'
                zi(jecono-1+igr) = 0
            end if
        end if
!
2       continue
    end do
    if (exixfm .eq. 'NON') goto 999
    if (model .eq. ' ') then
        call utmess('F', 'CALCULEL2_27')
    end if
!
!
!     2. CALCUL DE CHSGEO :
!     ------------------------------------------------------------------
    chgeom = '&&MANOPX.CHGEOM'
    lpain(1) = 'PGEOMER'
    lchin(1) = ma//'.COORDO'
    lpain(2) = 'PPINTTO'
    lchin(2) = model//'.TOPOSE.PIN'
    lpain(3) = 'PCNSETO'
    lchin(3) = model//'.TOPOSE.CNS'
    lpain(4) = 'PHEAVTO'
    lchin(4) = model//'.TOPOSE.HEA'
    lpain(5) = 'PLONCHA'
    lchin(5) = model//'.TOPOSE.LON'
    lpain(6) = 'PPMILTO'
    lchin(6) = model//'.TOPOSE.PMI'
    lpaout(1) = 'PXFGEOM'
    lchout(1) = chgeom
!
    call calcul('S', 'XFEM_XPG', ligrel, nbin, lchin, &
                lpain, nbout, lchout, lpaout, 'V', &
                'OUI')
    call celces(chgeom, 'V', chsgeo)
    call detrsd('CHAMP', chgeom)
!
999 continue
    call jedema()
end subroutine
