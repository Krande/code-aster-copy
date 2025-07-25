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

subroutine xnpgxx(model, ligrel, option, param, chsnpg, exixfm)

    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/indk32.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/modat2.h"
#include "asterfort/nbelem.h"
#include "asterfort/nucalc.h"
#include "asterfort/typele.h"
!
    character(len=8), intent(in) :: model
    character(len=19) :: ligrel, chsnpg
    character(len=16) :: option
    character(len=8) :: param
    character(len=3) :: exixfm
! ------------------------------------------------------------------
! BUT: REGARDER DANS LE LIGREL S'IL Y A DES ELEMENTS XFEM
!      ET SI LE CHAMP ELGA (OPTION/PARAM) UTILISE UNE FAMILLE XFEM..
!      SI FAMILLE XFEM :
!        * CALCULER LE CHAMP SIMPLE CHSNPG CONTENANT LE NOMBRE DE
!          POINTS DE GAUSS DES FAMILLES XFEM...
! ------------------------------------------------------------------
!     ARGUMENTS:
!     ----------
! LIGREL  IN/JXIN  K19 : LIGREL
! OPTION,PARAM  IN  K* : OPTION ET PARAMETRE PERMETTANT DE DETERMINER
!                        LA FAMILLE DE PG UTILISEE.
! EXIXFM  OUT K3 : 'OUI' : IL EXISTE DES GRELS AVEC FAMILLE XFEM...
!                  'NON' SINON
! CHSNPG  IN/JXOUT K19 : CHAM_ELEM_S (NEUT_R) DE TYPE 'ELEM', CONTENANT
!                        LE NOMBRE DE POINTS DE GAUSS DE LA FAMILLE 'XFEM'
!
! REMARQUE :
!   L'OBJET CHSGEO N'EST CREE QUE SI EXIXFM='OUI'
! ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbflmx
    parameter(nbflmx=20)
    character(len=8) :: lifapg(nbflmx)
!
    integer(kind=8) :: iopt, iopt1, nute, numc, igr, nbgrel
    integer(kind=8) :: imolo, jmolo, nec, kfpg, kfam
    integer(kind=8) :: igd, nblfpg, nbfam, nel, jliel, jfpgl, jcesdlon, jcesllon, jcesd, jcesl
    integer(kind=8) :: ndime, irese, nspg, nse, npg
    integer(kind=8) :: ima, iadlon, iad, iel
    integer(kind=8) :: k, nuflpg, nufgpg
    character(len=8) :: nomgd, elrese(6), elrefe, mesh, famil, noma
    character(len=8), pointer :: typma(:) => null()
    character(len=16) :: nofpg, nomte
    character(len=19) :: chslon
    character(len=24) :: chlong
    character(len=32) :: noflpg
    character(len=32), pointer :: pnlocfpg(:) => null()
    integer(kind=8), pointer :: nolocfpg(:) => null()
    integer(kind=8), pointer :: tmfpg(:) => null()
    integer(kind=8), pointer :: cesvlon(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
!     ------------------------------------------------------------------
    call jemarq()
!
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=mesh)

    exixfm = 'NON'

!   le modele comporte-t-il des elements X-FEM ?
    call dismoi('EXI_XFEM', ligrel, 'LIGREL', repk=exixfm)
    if (exixfm .eq. 'NON') then
        goto 999
    end if

    call jelira(ligrel//'.LIEL', 'NMAXOC', nbgrel)
!
    call jeveuo('&CATA.TE.PNLOCFPG', 'L', vk32=pnlocfpg)
    call jelira('&CATA.TE.NOLOCFPG', 'LONMAX', nblfpg)
    call jeveuo('&CATA.TE.NOLOCFPG', 'L', vi=nolocfpg)
!
    call jeveuo('&CATA.TM.TMFPG', 'L', vi=tmfpg)
    call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typma)
!
! -- construction d'un champ simple pour parcourir PLONCHA
    chlong = model//'.TOPOSE.LON'
    chslon = '&&XNPGSE.CHSLON'
    call celces(chlong, 'V', chslon)
!
    call jeveuo(chslon//'.CESD', 'L', jcesdlon)
    call jeveuo(chslon//'.CESL', 'L', jcesllon)
    call jeveuo(chslon//'.CESV', 'L', vi=cesvlon)
!
! -- allocation du CHAM_ELEM_S chsnpg
    call cescre('V', chsnpg, 'ELEM', mesh, 'NEUT_I', &
                0, ' ', [-1], [-1], [-1])
!
    call jeveuo(chsnpg//'.CESD', 'L', jcesd)
    call jeveuo(chsnpg//'.CESL', 'E', jcesl)
    call jeveuo(chsnpg//'.CESV', 'E', vi=cesv)
!
!     1. CALCUL DE EXIXFM :
!     ------------------------------------------------------------------
    call jenonu(jexnom('&CATA.OP.NOMOPT', 'XFEM_XPG'), iopt1)
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), iopt)
    do igr = 1, nbgrel
        nel = nbelem(ligrel, igr)
        call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', jliel)
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
!
!       L'ELEMENT SAIT-IL CALCULER XFEM_XPG ?
        numc = nucalc(iopt1, nute, 1)
        if (numc .lt. 0) cycle
!
        imolo = modat2(iopt, nute, param)
        if (imolo .eq. 0) cycle
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
                lifapg(k) = nofpg(9:16)
            end do
!
!       -- FAMILLE "ORDINAIRE"
        else
            nbfam = 1
            call jenuno(jexnum('&CATA.TM.NOFPG', kfpg), nofpg)
            lifapg(1) = nofpg(9:16)
        end if
!
!       --BOUCLE SUR LA/LES FAMILLE(S) :
        do kfam = 1, nbfam
            famil = lifapg(kfam)
!
            if (famil(1:4) .eq. 'XFEM') then
                exixfm = 'OUI'
!
!               RECHERCHE DE LA FAMILLE XINT
!
!               type de maille associe au type d'element
                noma = typma(nute)
!               recuperationn de la dimension topologique de la maille
                call dismoi('DIM_TOPO', noma, 'TYPE_MAILLE', repi=ndime)
!               calcul du decalage a appliquer si la maille est quadratique
                if (.not. ismali(noma)) then
                    irese = 3
                else
                    irese = 0
                end if
!
!               construction du nom de la famille de points de Gauss
!               pour le sous-element
                noflpg = nomte//elrese(ndime+irese)//'XINT'
!
!               recherche de cette famille dans la liste des familles
!               de points de Gauss
                nuflpg = indk32(pnlocfpg, noflpg, 1, nblfpg)
!
!               Assertion : on a trouve la famille
                ASSERT(nuflpg .ne. 0)

!               recuperation du nombre de point Gauss de la famille XINT,
!               i.e. du nombre de point de Gauss par sous-element
                nufgpg = nolocfpg(nuflpg)
                nspg = tmfpg(nufgpg)
!
!               calcul du nombre de points de Gauss pour chaque element
!               du groupe d'elements
                do iel = 1, nel
                    ima = zi(jliel-1+iel)
                    if (ima .lt. 0) cycle
!
!                   recuperation du nombre de sous-element de l'element
                    call cesexi('C', jcesdlon, jcesllon, ima, 1, 1, 1, iadlon)
                    ASSERT(iadlon .gt. 0)
                    nse = cesvlon(iadlon)

!                   calcul du nombre de points de Gauss de l'element
                    npg = nse*nspg

!                   stockage du nombre de points de Gauss de l'element
                    call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
                    iad = abs(iad)
                    zl(jcesl-1+iad) = .true.
                    cesv(iad) = npg
                end do
            end if
        end do
    end do
!
    call detrsd('CHAM_ELEM_S', chslon)
!
999 continue
    call jedema()
end subroutine
