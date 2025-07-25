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

subroutine irgmcg(chamsy, partie, ifi, nomcon, ordr, &
                  nbordr, coord, connx, point, nobj, &
                  nbel, nbcmpi, nomcmp, lresu, para, &
                  nomaou, versio)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/irgmg1.h"
#include "asterfort/irgmor.h"
#include "asterfort/irgmpv.h"
#include "asterfort/irgmtb.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=*) :: nomcon, chamsy, nomcmp(*), partie
    character(len=8) :: nomaou
    real(kind=8) :: coord(*), para(*)
    aster_logical :: lresu
    integer(kind=8) :: ifi, nbordr, nbcmpi
    integer(kind=8) :: versio
    integer(kind=8) :: ordr(*), connx(*), point(*)
!     NBRE, NOM D'OBJET POUR CHAQUE TYPE D'ELEMENT
    integer(kind=8) :: neletr
    parameter(neletr=8)
    integer(kind=8) :: ntyele, maxel, maxno
    parameter(ntyele=28)
    parameter(maxel=48)
    parameter(maxno=8)
    integer(kind=8) :: tdec(ntyele, maxel, maxno)
    integer(kind=8) :: typd(ntyele, 3)
    integer(kind=8) :: tord(neletr)
    integer(kind=8) :: nbel(ntyele), nbel2(ntyele), jel(ntyele)
    character(len=24) :: nobj(ntyele)
!
!        IMPRESSION D'UN CHAM_ELEM DE TYPE "ELGA","ELEM" AU FORMAT GMSH
!
!        CHAMSY : NOM DU CHAM_ELEM A ECRIRE
!        IFI    : NUMERO D'UNITE LOGIQUE DU FICHIER DE SORTIE GMSH
!        NOMCON : NOM DU CONCEPT A IMPRIMER
!        PARTIE : IMPRESSION DE LA PARTIE COMPLEXE OU REELLE DU CHAMP
!        NBORDR : NOMBRE DE NUMEROS D'ORDRE DANS LE TABLEAU ORDR
!        ORDR   : LISTE DES NUMEROS D'ORDRE A IMPRIMER
!        COORD  : VECTEUR COORDONNEES DES NOEUDS DU MAILLAGE
!        CONNX  : VECTEUR CONNECTIVITES DES NOEUDS DU MAILLAGE
!        POINT  : VECTEUR DU NOMBRE DE NOEUDS DES MAILLES DU MAILLAGE
!        NOBJ(i): NOM JEVEUX DEFINISSANT LES ELEMENTS DU MAILLAGE
!        NBEL(i): NOMBRE D'ELEMENTS DU MAILLAGE DE TYPE i
!        NBCMPI : NOMBRE DE COMPOSANTES DEMANDEES A IMPRIMER
!        NOMCMP : NOMS DES COMPOSANTES DEMANDEES A IMPRIMER
!        LRESU  : LOGIQUE INDIQUANT SI NOMCON EST UNE SD RESULTAT
!        PARA   : VALEURS DES VARIABLES D'ACCES (INST, FREQ)
!        NOMAOU : NOM DU MAILLAGE REDECOUPE
!        VERSIO : NUMERO DE LA VERSION GMSH UTILISEE (1 OU 2 )
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ior, i, j, k, ine, inoe, ima, listno(8), ix, nbno
    integer(kind=8) :: iq
    integer(kind=8) :: nbcmp, ipoin, iret, jcesc, jcesl
    integer(kind=8) :: jcesd, jtype
    integer(kind=8) :: icmp, ipt, isp, nbpt, nbsp, jnumol
    integer(kind=8) :: nbma, ncmpu, iad, nbcmpd, nbord2, iadmax, iadmm
    aster_logical :: iwri
    character(len=1) :: tsca
    character(len=8) :: k8b, nomgd, type, nocmp
    character(len=19) :: noch19, champs
    character(len=24) :: numold
    integer(kind=8), pointer :: cesc(:) => null()
    integer(kind=8), pointer :: cesd(:) => null()
    integer(kind=8), pointer :: cesl(:) => null()
    integer(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: vnocmp(:) => null()
    character(len=8), pointer :: cesk(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
! --- TABLEAU DES INFOS DE DECOUPAGE
    call irgmtb(tdec, typd, versio)
!
! --- ORDRE D'IMPRESSION DES VALEURS
    call irgmor(tord, versio)
!
    nbord2 = max(1, nbordr)
    numold = nomaou//'.NUMOLD         '
!
    AS_ALLOCATE(vi=cesd, size=nbord2)
    AS_ALLOCATE(vi=cesc, size=nbord2)
    AS_ALLOCATE(vi=cesv, size=nbord2)
    AS_ALLOCATE(vi=cesl, size=nbord2)
    call wkvect('&&IRGMCG.TYPE', 'V V K8', nbord2, jtype)
!
    nbcmp = 0
!
    do ior = 1, nbord2
        if (lresu) then
            call rsexch(' ', nomcon, chamsy, ordr(ior), noch19, &
                        iret)
            if (iret .ne. 0) goto 60
        else
            noch19 = nomcon
        end if
        call codent(ior, 'D0', k8b)
        champs = '&&IRGMCG.CH'//k8b
        call celces(noch19, 'V', champs)
        call jeveuo(champs//'.CESK', 'L', vk8=cesk)
        call jeveuo(champs//'.CESD', 'L', cesd(ior))
        call jeveuo(champs//'.CESC', 'L', cesc(ior))
        call jeveuo(champs//'.CESV', 'L', cesv(ior))
        call jeveuo(champs//'.CESL', 'L', cesl(ior))
        call jelira(champs//'.CESV', 'TYPE', cval=zk8(jtype+ior-1))
!
        nomgd = cesk(2)
        call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
        if (tsca .ne. 'R') then
            call utmess('F', 'ALGORITH2_63')
        end if
!
        type = cesk(3)
        if (type(1:4) .ne. 'ELGA' .and. type(1:4) .ne. 'ELEM') then
            call utmess('F', 'PREPOST2_56')
        end if
!
        if (ior .eq. 1) then
            jcesc = cesc(ior)
            jcesd = cesd(ior)
            jcesl = cesl(ior)
            nbma = zi(jcesd-1+1)
            nbcmp = zi(jcesd-1+2)
            ncmpu = 0
            AS_ALLOCATE(vk8=vnocmp, size=nbcmp)
            do icmp = 1, nbcmp
                do ima = 1, nbma
                    nbpt = zi(jcesd-1+5+4*(ima-1)+1)
                    nbsp = zi(jcesd-1+5+4*(ima-1)+2)
                    do ipt = 1, nbpt
                        do isp = 1, nbsp
                            call cesexi('C', jcesd, jcesl, ima, ipt, &
                                        isp, icmp, iad)
                            if (iad .gt. 0) goto 40
                        end do
                    end do
                end do
                goto 50
40              continue
                ncmpu = ncmpu+1
                vnocmp(ncmpu) = zk8(jcesc-1+icmp)
50              continue
            end do
        else
            if (zi(cesd(ior)-1+2) .ne. nbcmp) then
                call utmess('F', 'PREPOST2_53')
            end if
        end if
!
60      continue
    end do
!
! --- RECUPERATION DU TABLEAU DE CORRESPONDANCE ENTRE NUMERO DES
!     NOUVELLES MAILLES ET NUMERO DE LA MAILLE INITIALE
!     CREE PAR IRGMMA
!     ***************
    call jeveuo(numold, 'L', jnumol)
    do i = 1, ntyele
        nbel2(i) = 0
    end do
!
! --- BOUCLE SUR LE NOMBRE DE COMPOSANTES DU CHAM_ELEM
!     *************************************************
    if (nbcmpi .eq. 0) then
        nbcmpd = nbcmp
    else
        nbcmpd = nbcmpi
    end if
!
    do k = 1, nbcmpd
!
        if (nbcmpi .ne. 0) then
            do ix = 1, nbcmp
                if (vnocmp(ix) .eq. nomcmp(k)) then
                    icmp = ix
                    goto 80
                end if
            end do
            k8b = nomcmp(k)
            call utmess('F', 'PREPOST2_54', sk=k8b)
80          continue
        else
            icmp = k
        end if
        nocmp = vnocmp(icmp)
!
! ----- PREMIER PASSAGE POUR DETERMINER SI LE CHAMP A ECRIRE EXISTE
!       SUR LES POI1, SEG2, TRIA3, TETR4...
!       DONC ON  N'ECRIT RIEN
        iwri = .false.
!
! ----- BOUCLE SUR LES ELEMENTS DANS L'ORDRE DONNE PAR IRGMOR
!
        do ine = 1, neletr
!         I=NUM DE L'ELEMENT DANS LE CATALOGUE
            i = tord(ine)
            if (nbel(i) .ne. 0) then
                iadmm = 0
!           NBNO=NBRE DE NOEUDS DE CET ELEMENT
                nbno = typd(i, 3)
                call jeveuo(nobj(i), 'L', jel(i))
                do iq = 1, nbel(i)
                    ima = zi(jel(i)-1+iq)
                    call irgmg1(zi(jnumol), ima, nbord2, cesd, cesl, &
                                cesv, partie, jtype, nbno, icmp, &
                                ifi, iwri, iadmax)
                    iadmm = max(iadmax, iadmm)
                end do
                if (iadmm .gt. 0) nbel2(i) = nbel(i)
            end if
        end do
!
!
! ----- ECRITURE DE L'ENTETE DE View
!       ****************************
!
        call irgmpv(ifi, lresu, nomcon, chamsy, nbord2, &
                    para, nocmp, nbel2, .true._1, .false._1, &
                    .false._1, versio)
!
        iwri = .true.
!
! ----- BOUCLE SUR LES ELEMENTS DANS L'ORDRE DONNE PAR IRGMOR
!
        do ine = 1, neletr
!         I=NUM DE L'ELEMENT DANS LE CATALOGUE
            i = tord(ine)
            if (nbel2(i) .ne. 0) then
!           NBNO=NBRE DE NOEUDS DE CET ELEMENT
                nbno = typd(i, 3)
                call jeveuo(nobj(i), 'L', jel(i))
                do iq = 1, nbel(i)
                    ima = zi(jel(i)-1+iq)
                    ipoin = point(ima)
                    do inoe = 1, nbno
                        listno(inoe) = connx(ipoin-1+inoe)
                    end do
                    do j = 1, 3
                        write (ifi, 1000) (coord(3*(listno(inoe)-1)+j), &
                                           inoe=1, nbno)
                    end do
                    call irgmg1(zi(jnumol), ima, nbord2, cesd, cesl, &
                                cesv, partie, jtype, nbno, icmp, &
                                ifi, iwri, iadmax)
                end do
            end if
        end do
!
! ----- FIN D'ECRITURE DE View
!       **********************
!
        write (ifi, 1010) '$EndView'
!
    end do
!
    AS_DEALLOCATE(vi=cesc)
    AS_DEALLOCATE(vi=cesd)
    AS_DEALLOCATE(vi=cesv)
    AS_DEALLOCATE(vi=cesl)
    AS_DEALLOCATE(vk8=vnocmp)
    call jedetr('&&IRGMCG.TYPE')
    call jedema()
!
1000 format(1p, 4(e15.8, 1x))
1010 format(a8)
!
end subroutine
