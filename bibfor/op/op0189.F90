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

subroutine op0189()
! person_in_charge: sam.cuvilliez at edf.fr
!-----------------------------------------------------------------------
!
!     operateur : POST_VOISIN_CZM
!
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/cescar.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jeexin.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/knindi.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbgrel.h"
#include "asterfort/rsexch.h"
#include "asterfort/teattr.h"
#include "asterfort/typele.h"
#include "asterfort/trvocz.h"
#include "asterfort/utmess.h"
#include "asterfort/voiuti.h"
!
!-----------------------------------------------------------------------
!
    character(len=6) :: nompro
    parameter(nompro='OP0189')
!
    character(len=5) :: k5blan
    character(len=8) :: cartout, nomres, modele, mailla
    character(len=8) :: k8mod2, k8alia, k8prin, codvoi
    character(len=16) :: k16vge, k16tyel, k16bid
    character(len=16), pointer :: vvge(:) => null()
    character(len=19) :: ligrmo, cestmp
    character(len=24) :: k24ordr, chsief
!
    integer(kind=8) :: i, j, iret, iel, ima, nbordr, iorfin, jptvoi, jelvoi
    integer(kind=8) :: cpt, itypel, igrel, numa, numaco, ndim, dimmai, dimmod
    integer(kind=8) :: nbma, nbgrmo, nbel, nbelco, nbelma
    integer(kind=8) :: jrepmo, nbvois, nbvmas, numavo, jcesd, jcesl, iad
    integer(kind=8) :: jcelds, jcelvs
    integer(kind=8), pointer :: numeor(:) => null(), tab_typel(:) => null()
    integer(kind=8), pointer :: tab_tra(:) => null(), tab_numa(:) => null()
    integer(kind=8), pointer :: tab_maco(:) => null()
!
    real(kind=8) :: trxmoy
!
!   code des modelisations autorisees pour les elements "massifs"
    integer(kind=8) :: nmodmas
    parameter(nmodmas=4)
    character(len=3) :: lk3moma(nmodmas)
    character(len=8) :: lk8moma(nmodmas)
    data lk3moma/'CPL', 'DPL', 'AX_', '3D_'/
!
!   attributs 'TYPMOD2' des modelisations "cohesives"
    integer(kind=8) :: natczm
    parameter(natczm=2)
    character(len=8) :: lk8atczm(natczm)
    data lk8atczm/'ELEMJOIN', 'INTERFAC'/
!
!   pour la sd voisinage
    integer(kind=8) :: nvoima, nscoma
    parameter(nvoima=100, nscoma=4)
    integer(kind=8) :: livois(1:nvoima), tyvois(1:nvoima), nbnovo(1:nvoima)
    integer(kind=8) :: nbsoco(1:nvoima), lisoco(1:nvoima, 1:nscoma, 1:2)
    integer(kind=8) :: livmas(2)
    real(kind=8), pointer :: cesv(:) => null()
!
!-----------------------------------------------------------------------
!
    call jemarq()
!
! ---
!
    k5blan = '     '
    do i = 1, nmodmas
        lk8moma(i) = lk3moma(i)//k5blan
    end do
!
! --- recuperation du concept produit et des valeurs donnees aux mots-cles
!
    call getres(cartout, k16bid, k16bid)
    call getvid(' ', 'RESULTAT', scal=nomres)
!
! --- recuperation du dernier numero d'ordre "iorfin"
!
    k24ordr = nomres//'           .ORDR'
    call jeexin(k24ordr, iret)
    ASSERT(iret .ne. 0)
    call jelira(k24ordr, 'LONUTI', ival=nbordr)
    ASSERT(nbordr .gt. 0)
    call jeveuo(k24ordr, 'L', vi=numeor)
    iorfin = numeor(nbordr)
!
! --- recuperation du champ SIEF_ELGA au dernier numero d'ordre
!
    call rsexch('F', nomres, 'SIEF_ELGA', iorfin, chsief, iret)
    call jeveuo(chsief(1:19)//'.CELD', 'L', jcelds)
    call jeveuo(chsief(1:19)//'.CELV', 'L', jcelvs)
!
! --- recuperation des objets utiles (modele, maillage...)
!
    call dismoi('NOM_MODELE', nomres, 'RESULTAT', repk=modele)
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mailla)

    call jeveuo(modele//'.MAILLE', 'L', vi=tab_typel)
    call jelira(modele//'.MAILLE', 'LONUTI', ival=nbma)
    ligrmo = modele//'.MODELE'
    call jeveuo(ligrmo//'.REPE', 'L', jrepmo)
!
! --- recuperation de la dimension -> code des types de voisin
!
    call dismoi('DIM_GEOM', mailla, 'MAILLAGE', repi=dimmai)
    call dismoi('DIM_GEOM', modele, 'MODELE', repi=dimmod)
    ASSERT(dimmai .eq. dimmod)
    ndim = dimmai
!
    if (ndim .eq. 2) then
        codvoi = 'A2'
    else if (ndim .eq. 3) then
        codvoi = 'F3'
    else
        ASSERT(.false.)
    end if
!
! --- creation du cham_elem_s dans lequel on va ecrire par maille cohesive
! --- la triaxialite (moyennee) des elements massifs voisins
!
    cestmp = '&&'//nompro//'.CESTMP'
    call cescre('V', cestmp, 'ELEM', mailla, 'NEUT_R', &
                1, 'X1', [0], [-1], [-1])
    call jeveuo(cestmp//'.CESD', 'L', jcesd)
    call jeveuo(cestmp//'.CESV', 'E', vr=cesv)
    call jeveuo(cestmp//'.CESL', 'E', jcesl)
!
! --- allocation du tableau de travail "tab_tra"
!
    AS_ALLOCATE(vi=tab_tra, size=nbma)
    tab_tra(:) = 0
!
! -------------------------------------------------------------------
!     boucle sur les grels de ligrmo (1ere passe)
! -------------------------------------------------------------------
!
    nbelco = 0
    nbelma = 0

    nbgrmo = nbgrel(ligrmo)
    do igrel = 1, nbgrmo

        itypel = typele(ligrmo, igrel)
        call jenuno(jexnum('&CATA.TE.NOMTE', itypel), k16tyel)
!
!       ce sont des elements cohesifs
        call teattr('C', 'TYPMOD2', k8mod2, iret, typel=k16tyel)
        if (iret .eq. 0) then
            if (knindi(8, k8mod2, lk8atczm, natczm) .gt. 0) then

                nbelco = nbelco+nbelem(ligrmo, igrel)

            end if
        end if

!       ce sont des elements "massifs"
        call teattr('S', 'ALIAS8', k8alia, iret, typel=k16tyel)
        if (knindi(8, k8alia(3:5)//k5blan, lk8moma, nmodmas) .gt. 0) then

            call teattr('S', 'PRINCIPAL', k8prin, iret, typel=k16tyel)
            if (k8prin(1:3) .eq. 'OUI') then

                call jeveuo(jexnum(ligrmo//'.LIEL', igrel), 'L', vi=tab_numa)
                do iel = 1, nbelem(ligrmo, igrel)
                    numa = tab_numa(iel)
                    tab_tra(numa) = 1
                end do

                nbelma = nbelma+nbelem(ligrmo, igrel)

            end if

        end if

    end do
!
    if (nbelco .eq. 0) call utmess('F', 'PREPOST6_46', sk=modele)
    ASSERT(nbelma .gt. 0)
!
! --- allocation de la liste des mailles portant des elements cohesifs
!
    AS_ALLOCATE(vi=tab_maco, size=nbelco)
!
! -------------------------------------------------------------------
!     boucle sur les grels de ligrmo (2eme passe)
! -------------------------------------------------------------------
!
    cpt = 0
!
    do igrel = 1, nbgrmo

        itypel = typele(ligrmo, igrel)
        call jenuno(jexnum('&CATA.TE.NOMTE', itypel), k16tyel)
!
!       ce sont des elements cohesifs
        call teattr('C', 'TYPMOD2', k8mod2, iret, typel=k16tyel)
        if (iret .eq. 0) then
            if (knindi(8, k8mod2, lk8atczm, natczm) .gt. 0) then

                nbel = nbelem(ligrmo, igrel)
                call jeveuo(jexnum(ligrmo//'.LIEL', igrel), 'L', vi=tab_numa)
                tab_maco(cpt+1:cpt+nbel) = tab_numa(1:nbel)

                cpt = cpt+nbel

            end if
        end if

    end do
!
! -------------------------------------------------------------------
!                boucle sur les mailles cohesives
! -------------------------------------------------------------------
!
! --- recuperation de la sd_voisinage
!
    call jeveuo(ligrmo//'.NVGE', 'L', vk16=vvge)
    k16vge = vvge(1)
    call jeveuo(k16vge(1:12)//'.PTVOIS', 'L', jptvoi)
    call jeveuo(k16vge(1:12)//'.ELVOIS', 'L', jelvoi)
!
    do ima = 1, nbelco

        numaco = tab_maco(ima)

        call voiuti(numaco, codvoi, nvoima, nscoma, jrepmo, &
                    jptvoi, jelvoi, nbvois, livois, tyvois, &
                    nbnovo, nbsoco, lisoco)

        nbvmas = 0
        do j = 1, nbvois
            numavo = livois(j)
            if (tab_tra(numavo) .eq. 1) then
                nbvmas = nbvmas+1
                livmas(nbvmas) = numavo
            end if
        end do
!       1 ou 2 voisin(s) "massif" pour une maille cohesive
        ASSERT(nbvmas .eq. 1 .or. nbvmas .eq. 2)

!       calcul de la triaxialite moyennee sur les voisins "massif"
        call trvocz(ndim, nbvmas, livmas, jrepmo, &
                    jcelds, jcelvs, trxmoy)
!
        call cesexi('S', jcesd, jcesl, numaco, 1, 1, 1, iad)
        ASSERT(iad .lt. 0)
        iad = -iad
        zl(jcesl-1+iad) = .true.
        cesv(iad) = trxmoy

    end do
!
! --- creation de la carte en sortie
!
    call cescar(cestmp, cartout, 'G')
!
! --- menage
!
    AS_DEALLOCATE(vi=tab_tra)
    AS_DEALLOCATE(vi=tab_maco)
    call detrsd('CHAM_ELEM_S', cestmp)
!
    call jedema()
!
end subroutine
