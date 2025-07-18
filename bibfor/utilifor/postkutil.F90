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

subroutine postkutil(imater, nomres, nomfis, repmat, repmod)
    implicit none
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cncinv.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
#include "asterfort/xtmafi.h"
#include "asterfort/char8_to_int.h"
    integer(kind=8), intent(in) :: imater
    character(len=*), intent(in) :: nomres
    character(len=*), intent(in) :: nomfis
    character(len=*), intent(out) :: repmat
    character(len=*), intent(out) :: repmod
!
! person_in_charge: sam.cuvilliez at edf.fr
! ----------------------------------------------------------------------
!
! in:
!   resu   : nom d'une sd_resultat
!   nomfis : nom d'une sd_fond_fissure ou d'une sd_fiss_xfem
!   imater : 1 => on recherche une sd_mater
!            0 => on ne recherche pas
! out:
!   repmat : nom d'une sd_mater
!   repmod : "nom" d'une modelisation parmi '3D', 'AXIS', 'D_PLAN', 'C_PLAN'
!
! ------------------------------------------------------------------
!
! but : recuperer un nom de modélisation (parmi '3D', 'AXIS', 'D_PLAN',
! ---   'C_PLAN') et le nom d'une sd_mater dans la sd_resultat resu,
!       connaissant nomfis (sd_fond_fiss en fem ou sd_fiss_xfem en xfem)
!
!       On recupere la liste des mailles "voisines" du fond de fissure
!       (pour fem les mailles connectees aux noeuds du fond de fissure,
!       et pour xfem les mailles CTIP). On s'assure que la meme sd_mater
!       et la meme modelisation ont ete affectees sur ces mailles.
!
!       On renvoie le nom de cette sd_mater et le nom de cette modelisation.
!       Cette routine est utilisee depuis le corps de la macro POST_K1_K2_K3
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbma, ier, nfiss, iret, ima, nutyel, iad, i
    integer(kind=8) :: imodeli, ndim, jcesl, jcesd, nbcmp, icmp, vin_modeli(4)
    integer(kind=8) :: ier1, nbnof, ino, inoeu
    integer(kind=8) :: nmanoe, jmanoe, imanoe, itypma, ndime
    integer(kind=8) :: cpt_ma_fon, cpt_ma_noe, imaf, j, ima1, ima2, cpt_dbl
    integer(kind=8) :: nbma_tmp, nbma_fon
    character(len=8) :: nommod, nomchm, noma, vk8_typmod(4)
    character(len=8) :: vk8_mater(5)
    character(len=8) :: k8typmo, k8mater, k8noeu, k8typma
    character(len=16) :: ktyel
    character(len=19) :: chmat, cesmat, cnxinv
    character(len=24) :: mesmai, limafo
    aster_logical :: l_xfem, inco
    integer(kind=8), pointer :: vmatmp(:) => null()
    integer(kind=8), pointer :: vmafon(:) => null()
    integer(kind=8), pointer :: vtyele(:) => null()
    integer(kind=8), pointer :: vtypma(:) => null()
    character(len=8), pointer :: vcesv(:) => null()
    character(len=8), pointer :: vcesk(:) => null()
    character(len=8), pointer :: v8fiss(:) => null()
    character(len=8), pointer :: vnofon(:) => null()
!
    data vk8_typmod/'3D', 'AXIS', 'D_PLAN', 'C_PLAN'/
    data vin_modeli/3, 2, 2, 2/
!
!   --------------------------------------------------------------------
!   debut
!   --------------------------------------------------------------------
!
    call jemarq()
!
!   --------------------------------------------------------------------
!   prealables
!   --------------------------------------------------------------------
!
!   recup du modele et du vecteur .MAILLE
    call dismoi('NOM_MODELE', nomres, 'RESULTAT', repk=nommod)
    ASSERT((nommod .ne. '#AUCUN') .or. (nommod .ne. '#PLUSIEURS'))
    call jeveuo(nommod//'.MAILLE', 'E', vi=vtyele)
!
!   recup du maillage de definition du modele
    call dismoi('NOM_MAILLA', nommod, 'MODELE', repk=noma)
    call dismoi('NB_MA_MAILLA', noma, 'MAILLAGE', repi=nbma)
!
!   recup de la dimension du probleme
    call dismoi('DIM_GEOM', nommod, 'MODELE', repi=ndim)
    ASSERT((ndim .eq. 2) .or. (ndim .eq. 3))
!
!   recup de la carte cham_mater//'.CHAMP_MAT' sous forme de cham_elem_s
    call dismoi('CHAM_MATER', nomres, 'RESULTAT', repk=nomchm)
    ASSERT((nomchm .ne. '#AUCUN') .or. (nomchm .ne. '#PLUSIEURS'))
!
    chmat = nomchm//'.CHAMP_MAT'
    cesmat = '&&POSTKUTIL.CESMAT'
    call carces(chmat, 'ELEM', ' ', 'V', cesmat, 'A', iret)
    ASSERT(iret .eq. 0)
!
    call jeveuo(cesmat//'.CESV', 'L', vk8=vcesv)
    call jeveuo(cesmat//'.CESL', 'L', jcesl)
    call jeveuo(cesmat//'.CESD', 'L', jcesd)
    call jeveuo(cesmat//'.CESK', 'L', vk8=vcesk)
    nbcmp = zi(jcesd-1+2)
!
!   verifs sur le cham_mater
    ASSERT(vcesk(1) .eq. noma)
    ASSERT(vcesk(3) .eq. 'ELEM')
    ASSERT(zi(jcesd-1+1) .eq. nbma)
    ASSERT((nbcmp .eq. 1) .or. (nbcmp .eq. 5))
    ASSERT(zi(jcesd-1+3) .eq. 1)
    ASSERT(zi(jcesd-1+4) .eq. 1)
!
!   s'agit-il d'un modele fem ou x-fem
    l_xfem = .false.
    call jeexin(nommod//'.FISS', ier)
    if (ier .ne. 0) then
        l_xfem = .true.
    end if
!
!   --------------------------------------------------------------------
!   recup de la liste vmafon des mailles principales CTIP dans le cas xfem
!   --------------------------------------------------------------------
!
    if (l_xfem) then
!
!       verif sur la sd_fiss_xfem in
        call dismoi('NB_FISS_XFEM', nommod, 'MODELE', repi=nfiss)
        call jeveuo(nommod//'.FISS', 'L', vk8=v8fiss)
        ASSERT(indik8(v8fiss, nomfis, 1, nfiss) .gt. 0)
!
!       recup de la liste des mailles CTIP pour la fissure nomfis
        limafo = '&&XMOT2M.NUM_MAILLES'
        mesmai = '&&XMOT2M.MES_MAILLES'
        call xtmafi(ndim, [nomfis], 1, limafo, mesmai, nbma_fon, model=nommod, typ_enr='CTIP')
        call jeveuo(limafo, 'L', vi=vmafon)
!
!       menage
        call jedetr(mesmai)
!
!   --------------------------------------------------------------------
!   recup de la liste vmafon des mailles principales connectees aux noeuds
!   du fond dans le cas fem
!   --------------------------------------------------------------------
!
    else
!
!       creation de la connectivite inverse
        cnxinv = '&&POSTKUTIL.CNXINV'
        call cncinv(noma, [0], 0, 'V', cnxinv)
!
!       recup du vecteur .TYPMAIL
        call jeveuo(noma//'.TYPMAIL', 'L', vi=vtypma)
!
!       recup de la liste des noeuds du fond
        call jeexin(nomfis//'.FOND.NOEU', ier1)
        ASSERT(ier1 .gt. 0)
        if (ier1 .ne. 0) then
            call jeveuo(nomfis//'.FOND.NOEU', 'L', vk8=vnofon)
        end if
        nbnof = size(vnofon)
!
!       1ere boucle sur les noeuds du fond : dimensionnement de la liste
!       brute des mailles principales connectees aux noeud du fond
        cpt_ma_fon = 0
        do ino = 1, nbnof
!
!           recup des mailles connectees au noeud courant
            k8noeu = vnofon(ino)
            inoeu = char8_to_int(k8noeu)
            inoeu = char8_to_int(k8noeu)
            call jelira(jexnum(cnxinv, inoeu), 'LONMAX', nmanoe)
            call jeveuo(jexnum(cnxinv, inoeu), 'L', jmanoe)
!           on s'assure que le noeud n'est pas orphelin
            ASSERT(zi(jmanoe-1+1) .ne. 0)
!
!           boucle sur les mailles connectees au noeud courant
            cpt_ma_noe = 0
            do ima = 1, nmanoe
!
!               on ne garde que les mailles principales
                imanoe = zi(jmanoe-1+ima)
                itypma = vtypma(imanoe)
                call jenuno(jexnum('&CATA.TM.NOMTM', itypma), k8typma)
                call dismoi('DIM_TOPO', k8typma, 'TYPE_MAILLE', repi=ndime)
                if (ndime .eq. ndim) then
                    cpt_ma_noe = cpt_ma_noe+1
                end if
!
            end do
            ASSERT(cpt_ma_noe .gt. 0)
            cpt_ma_fon = cpt_ma_fon+cpt_ma_noe
!
        end do
        nbma_tmp = cpt_ma_fon
        AS_ALLOCATE(vi=vmatmp, size=nbma_tmp)
!
!       2eme boucle sur les noeuds du fond : remplissage de la liste
!       brute des mailles principales connectees aux noeud du fond
        cpt_ma_fon = 0
        do ino = 1, nbnof
!
!           recup des mailles connectees au noeud courant
            k8noeu = vnofon(ino)
            inoeu = char8_to_int(k8noeu)
            call jelira(jexnum(cnxinv, inoeu), 'LONMAX', nmanoe)
            call jeveuo(jexnum(cnxinv, inoeu), 'L', jmanoe)
!           on s'assure que le noeud n'est pas orphelin
            ASSERT(zi(jmanoe-1+1) .ne. 0)
!
!           boucle sur les mailles connectees au noeud courant
            do ima = 1, nmanoe
!
!               on ne garde que les mailles principales
                imanoe = zi(jmanoe-1+ima)
                ASSERT(imanoe .gt. 0)
                itypma = vtypma(imanoe)
                call jenuno(jexnum('&CATA.TM.NOMTM', itypma), k8typma)
                call dismoi('DIM_TOPO', k8typma, 'TYPE_MAILLE', repi=ndime)
                if (ndime .eq. ndim) then
                    cpt_ma_fon = cpt_ma_fon+1
                    vmatmp(cpt_ma_fon) = imanoe
                end if
!
            end do
!
        end do
!
!       suppression des doublons dans la liste brute et allocation
!       de liste sans doublons
        cpt_dbl = 0
        do i = 1, nbma_tmp
            ima1 = vmatmp(i)
            if (ima1 .ne. 0) then
                do j = i+1, nbma_tmp
                    ima2 = vmatmp(j)
                    if ((ima1 .eq. ima2) .and. (ima2 .ne. 0)) then
                        vmatmp(j) = 0
                        cpt_dbl = cpt_dbl+1
                    end if
                end do
            end if
        end do
        nbma_fon = nbma_tmp-cpt_dbl
        ASSERT(nbma_fon .gt. 0)
        AS_ALLOCATE(vi=vmafon, size=nbma_fon)
!
!       remplissage de la liste sans doublons
        cpt_ma_fon = 0
        do i = 1, nbma_tmp
            if (vmatmp(i) .ne. 0) then
                cpt_ma_fon = cpt_ma_fon+1
                vmafon(cpt_ma_fon) = vmatmp(i)
            end if
        end do
!
!       menage
        call jedetr(cnxinv)
        AS_DEALLOCATE(vi=vmatmp)
!
    end if
!
!   --------------------------------------------------------------------
!   recup du nom du materiau et de la modelisation sur les mailles
!   principales situees "au voisinage du fond de la fissure" (fem ou xfem)
!   --------------------------------------------------------------------
!
!   recup de la modelisation sur la 1ere maille de vmafon
    imodeli = 0
    imaf = vmafon(1)
    nutyel = vtyele(imaf)
    call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), ktyel)

    do i = 1, 4
        if (lteatt('TYPMOD', vk8_typmod(i), ktyel)) then
            imodeli = i
            exit
        end if
    end do
    inco = .False.
! .or. lteatt('INCO','C3B') .or. lteatt('INCO','C2 ') .or. lteatt('INCO','C2O')
    if (lteatt('INCO', 'C5GV', ktyel) .or. &
        lteatt('INCO', 'C3B', ktyel) .or. &
        lteatt('INCO', 'C2 ', ktyel) .or. &
        lteatt('INCO', 'C2O', ktyel)) then
        inco = .True.
    end if

    if (inco) then
        call utmess('F', 'RUPTURE1_74')
    end if
    ASSERT(imodeli .gt. 0)
    k8typmo = vk8_typmod(imodeli)

    ASSERT(ndim .eq. vin_modeli(imodeli))
!
!   recup du materiau sur la 1ere maille de vmafon
    vk8_mater(:) = ''
    imaf = vmafon(1)
    do icmp = 1, nbcmp
        call cesexi('S', jcesd, jcesl, imaf, 1, 1, icmp, iad)
        ASSERT(iad .gt. 0)
        vk8_mater(icmp) = vcesv(iad)
    end do
    if (nbcmp .gt. 1) then
        ASSERT(vk8_mater(2) .eq. 'TREF=>')
    end if
    k8mater = vk8_mater(1)
!
!   boucle sur les mailles de vmafon
    do ima = 1, nbma_fon
!
        imaf = vmafon(ima)
!
!       verif que toutes les mailles de vmafon sont affectees avec
!       la meme modelisation
        nutyel = vtyele(imaf)
        call jenuno(jexnum('&CATA.TE.NOMTE', nutyel), ktyel)
        ASSERT(lteatt('TYPMOD', k8typmo, ktyel))
!
!       verif que toutes les mailles de vmafon sont affectees avec
!       le meme materiau
        if (imater .eq. 1) then
            do icmp = 1, nbcmp
                call cesexi('S', jcesd, jcesl, imaf, 1, 1, icmp, iad)
                ASSERT(iad .gt. 0)
                if (vk8_mater(icmp) .ne. vcesv(iad)) then
                    call utmess('F', 'RUPTURE1_13')
                end if
            end do
        end if
!
    end do
!
!   --------------------------------------------------------------------
!   variables out
!   --------------------------------------------------------------------
!
    if (imater .eq. 1) then
        repmat = k8mater
    else
        repmat = ' '
    end if
    repmod = k8typmo
!
!   --------------------------------------------------------------------
!   menage
!   --------------------------------------------------------------------
!
    if (l_xfem) then
        call jedetr(limafo)
    else
        AS_DEALLOCATE(vi=vmafon)
    end if
!
    call detrsd('CHAM_ELEM_S', cesmat)
!
!   --------------------------------------------------------------------
!   fin
!   --------------------------------------------------------------------
!
    call jedema()
!
end subroutine
