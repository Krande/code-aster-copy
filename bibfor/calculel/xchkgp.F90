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

subroutine xchkgp(model)
! person_in_charge: patrick.massin at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/celces.h"
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
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8) :: model
! ------------------------------------------------------------------
! BUT: VERIFIER QUE LE CHAMP ELGA (OPTION/PARAM), UTILISE UNE FAMILLE
!      XFEM CONTENANT SUFFISAMENT DE POINTS POUR STOCKER LES DONNÉES
!      DES TOUS LES SOUS-ELEMENTS
! ------------------------------------------------------------------
!     ARGUMENTS:
!     ----------
! MODEL  IN/JXIN  K8   : MODELE
! OPTION,PARAM  IN  K* : OPTION ET PARAMETRE PERMETTANT DE DETERMINER
!                        LA FAMILLE DE PG UTILISEE.
! ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: nbflmx
    parameter(nbflmx=20)
    character(len=8) :: lifapg(nbflmx)
    integer(kind=8) :: linbpg(nbflmx)
!
    integer(kind=8) :: iopt, nute, igr, nbgrel
    integer(kind=8) :: imolo, jmolo, nec, kfpg, kfam
    integer(kind=8) :: igd, nblfpg, nbfam, nel, jliel, jfpgl, jcesdlon, jcesllon
    integer(kind=8) :: ndime, irese, nspg, nse, npg, npgfam
    integer(kind=8) :: ima, iadlon, iel
    integer(kind=8) :: k, nuflpg, nufgpg, vali(2)
    character(len=8) :: nomgd, elrese(6), elrefe, ma, famil, noma, nomail, param
    character(len=8), pointer :: typma(:) => null()
    character(len=16) :: nofpg, nomte, valk(2), option, pheno
    character(len=19) :: ligrel, chslon
    character(len=24) :: chlong
    character(len=32) :: noflpg
    character(len=32), pointer :: pnlocfpg(:) => null()
    integer(kind=8), pointer :: nolocfpg(:) => null()
    integer(kind=8), pointer :: tmfpg(:) => null()
    integer(kind=8), pointer :: cesvlon(:) => null()
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
!     ------------------------------------------------------------------
    call jemarq()
!
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=ma)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrel)
    call dismoi('PHENOMENE', model, 'MODELE', repk=pheno)
!
    call jelira(ligrel//'.LIEL', 'NMAXOC', nbgrel)
!
    call jeveuo('&CATA.TE.PNLOCFPG', 'L', vk32=pnlocfpg)
    call jelira('&CATA.TE.NOLOCFPG', 'LONMAX', nblfpg)
    call jeveuo('&CATA.TE.NOLOCFPG', 'L', vi=nolocfpg)
!
    call jeveuo('&CATA.TM.TMFPG', 'L', vi=tmfpg)
    call jeveuo('&CATA.TE.TYPEMA', 'L', vk8=typma)
!
!   définition de l'option et du paramètre permettant d'identifier
!   un champ porté par la famille XFEM de l'élément
    option = 'TOU_INI_ELGA'
    param = 'PSIEF_R'
    if (pheno .eq. 'THERMIQUE') then
        option = 'TEMP_ELGA'
        param = 'PTEMP_R'
    end if
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
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), iopt)
    do igr = 1, nbgrel
        nel = nbelem(ligrel, igr)
        call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', jliel)
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
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
!               stockage du nom de la famille courante
                call jenuno(jexnum('&CATA.TM.NOFPG', nufgpg), nofpg)
                lifapg(k) = nofpg(9:16)
!               stockage du nombre de points de la famille courante
                linbpg(k) = tmfpg(nufgpg)
            end do
!
!       -- FAMILLE "ORDINAIRE"
        else
            nbfam = 1
!           stockage du nom de la famille
            call jenuno(jexnum('&CATA.TM.NOFPG', kfpg), nofpg)
            lifapg(1) = nofpg(9:16)
!           stockage du nombre de points de la famille
            linbpg(1) = tmfpg(kfpg)
        end if
!
!       --BOUCLE SUR LA/LES FAMILLE(S) :
        do kfam = 1, nbfam
            famil = lifapg(kfam)
            npgfam = linbpg(kfam)
!
            if (famil(1:4) .eq. 'XFEM') then

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

!               recuperation du nombre de points de Gauss de la famille XINT,
!               i.e. du nombre de points de Gauss par sous-element
                nufgpg = nolocfpg(nuflpg)
                nspg = tmfpg(nufgpg)
!
!               calcul du nombre de points de Gauss pour chaque element
!               du groupe d'elements
                do iel = 1, nel
                    ima = zi(jliel-1+iel)
                    if (ima .le. 0) cycle
!
!                   recuperation du nombre de sous-element de l'element
                    call cesexi('C', jcesdlon, jcesllon, ima, 1, 1, 1, iadlon)
                    ASSERT(iadlon .gt. 0)
                    nse = cesvlon(iadlon)

!                   calcul du nombre de points de Gauss de l'element
                    npg = nse*nspg

!                   erreur fatale si le nombre de points de Gauss de l'élément
!                   dépasse le nombre de points de Gauss de la famille XFEM
                    if (npg .gt. npgfam) then
!                      récupération du nom de la maille
                        nomail = int_to_char8(ima)
!
                        vali(1) = npgfam
                        vali(2) = npg
                        valk(1) = nomte
                        valk(2) = nomail
                        call utmess('F', 'XFEM_56', nk=2, valk=valk, ni=2, vali=vali)
                    end if
!
                end do
            end if
        end do
    end do
!
    call detrsd('CHAM_ELEM_S', chslon)
!
    call jedema()
end subroutine
