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

subroutine verif_affe_carte(ligrmo, carte, comment, non_lin)
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jenonu.h"
#include "asterfort/jelira.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/jenuno.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/list_grma.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/etenca.h"
#include "asterfort/jeexin.h"
#include "asterfort/int_to_char8.h"
!
    character(len=19), intent(in) :: ligrmo
    character(len=19), intent(in) :: carte
    character(len=*), intent(in) :: comment
    aster_logical, intent(in), optional ::  non_lin
!
!-----------------------------------------------------------------------
!   But :
!     Emettre des alarmes concernant les affectations douteuses d'une carte
!
!   Entrees:
!     ligrmo     :  ligrel du modele
!     carte      :  sd_carte
!     comment    :  commentaire de description pour la sd_carte
!
!-----------------------------------------------------------------------
    character(len=3) :: tsca
    character(len=8) :: nomgd, nocmp, mailla, nommai
    character(len=24) :: lgrma(4)
    character(len=16) :: nomte
    character(len=80) :: valk(5)
    integer(kind=8) :: nbgrel, igrel, kcmp, nbcmp, nbop, nbte, k1, iexi, ient, iad1
    integer(kind=8) :: jnocmp, numgd, joptte, jligrmo, n1, kop, ioptte, joptmod, jvale
    integer(kind=8) :: jmodeloc, nbin, kin, moloc, nbma, ima, iret, te, nbmapb, nbgrma
    integer(kind=8) :: nucalc, k, kma, nec, nbma_verif, nbgdmx, code, decal, ico, kcmp2
    integer(kind=8), pointer :: a_un_sens(:) => null()
    integer(kind=8), pointer :: num_grel(:) => null()
    integer(kind=8), pointer :: numa_verif(:) => null()
    integer(kind=8), pointer :: desc(:) => null()
    integer(kind=8), pointer :: ptma(:) => null()
    integer(kind=8), pointer :: dg(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8)          :: list_ma_pb(5), typq4, typt3
    aster_logical    :: verif_coef_drz
    aster_logical    :: exiq4_drz_nook, exiq4_coef_drz
    aster_logical    :: exiq3_coef_drz

!-----------------------------------------------------------------------
!
    call jemarq()

    verif_coef_drz = ASTER_FALSE
    exiq4_drz_nook = ASTER_FALSE
    exiq4_coef_drz = ASTER_FALSE
    exiq3_coef_drz = ASTER_FALSE

    call dismoi('NOM_GD', carte, 'CARTE', repk=nomgd)
    call dismoi('NB_CMP_MAX', nomgd, 'GRANDEUR', repi=nbcmp)
    call dismoi('NUM_GD', nomgd, 'GRANDEUR', repi=numgd)
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', jnocmp)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
    call dismoi('NB_GREL', ligrmo, 'LIGREL', repi=nbgrel)
    call dismoi('NOM_MAILLA', ligrmo, 'LIGREL', repk=mailla)
    call dismoi('NB_MA_MAILLA', mailla, 'MAILLAGE', repi=nbma)

!   -- 1. recuperation des objets des catalogues EF:
!   -------------------------------------------------
    call jelira('&CATA.OP.NOMOPT', 'NOMMAX', nbop)
    call jelira('&CATA.TE.NOMTE', 'NOMMAX', nbte)
    call jeveuo('&CATA.TE.OPTTE', 'L', joptte)

!   -- 2. On calcule 2 tableaux :
!         A_UN_SENS(igrel,kcmp)  : si 0=non , si 1=oui
!         NUM_GREL(2*(ima-1)+1)  : igrel associe a la maille ima
!         NUM_GREL(2*(ima-1)+2)  : te    associe a la maille ima
!   --------------------------------------------------------------------
    AS_ALLOCATE(vi=a_un_sens, size=nbgrel*nbcmp)
    a_un_sens = 0
    AS_ALLOCATE(vi=num_grel, size=2*nbma)
    num_grel = 0

    do igrel = 1, nbgrel
        call jeveuo(jexnum(ligrmo//'.LIEL', igrel), 'L', jligrmo)
        call jelira(jexnum(ligrmo//'.LIEL', igrel), 'LONMAX', n1)
        te = zi(jligrmo-1+n1)
! ----- Loop on elements
        do k = 1, n1-1
            ima = zi(jligrmo-1+k)
            if (ima .gt. 0) then
                num_grel(2*(ima-1)+1) = igrel
                num_grel(2*(ima-1)+2) = te
            end if
        end do
        do kop = 1, nbop
            ioptte = zi(joptte-1+(te-1)*nbop+kop)
            if (ioptte .gt. 0) then
                call jeveuo(jexnum('&CATA.TE.OPTMOD', ioptte), 'L', joptmod)
                nucalc = zi(joptmod-1+1)
                if (nucalc .ge. 0) then
                    nbin = zi(joptmod-1+2)
                    do kin = 1, nbin
                        moloc = zi(joptmod-1+3+kin)
                        call jeveuo(jexnum('&CATA.TE.MODELOC', moloc), 'L', jmodeloc)
                        if (zi(jmodeloc-1+2) .eq. numgd) then
                            ASSERT(zi(jmodeloc-1+4) .gt. 0)
                            do kcmp = 1, nbcmp
                                if (exisdg(zi(jmodeloc-1+5), kcmp)) then
                                    a_un_sens((igrel-1)*nbcmp+kcmp) = 1
                                end if
                            end do
                        end if
                    end do
                end if
            end if
        end do
    end do
!
!   3. On parcourt les CMPS affectees volontairement dans la carte (sauf TOUT='OUI')
!
!   3.1 : on repère les mailles a vérifier (numa_verif(*)) :
    call jeveuo(carte//'.DESC', 'L', vi=desc)
    call jeveuo(carte//'.VALE', 'L', jvale)
    nbgdmx = desc(2)
!   si la carte est constante (TOUT='OUI'), on ne vérifie pas
    if (nbgdmx .eq. 1 .and. desc(3+1) .eq. 1) goto 999
!
    call etenca(carte, ligrmo, iret)
    if (iret .gt. 0) goto 999
    call jeexin(carte//'.PTMA', iexi)
    if (iexi .eq. 0) goto 999
!
    nbma_verif = 0
    AS_ALLOCATE(vi=numa_verif, size=nbma)
    call jeveuo(carte//'.PTMA', 'L', vi=ptma)
    do ima = 1, nbma
        ient = ptma(ima)
        ! si la maille n'est pas affectee :
        if (ient .eq. 0) cycle
        ! si la maille est affectee par TOUT='OUI' :
        code = desc(3+2*(ient-1)+1)
        if (code .eq. 1) cycle
        ! on ne vérifie pas les mailles qui ne sont pas affectées dans le modèle
        ! (on ne saurait pas remplir le champ nomte du message)
        igrel = num_grel(2*(ima-1)+1)
        if (igrel .eq. 0) cycle
        !
        nbma_verif = nbma_verif+1
        numa_verif(nbma_verif) = ima
    end do
    !
    call jeveuo(mailla//'.TYPMAIL', 'L', vi=typmail)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD4'), typq4)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA3'), typt3)
!
!   3.2 : on verifie les mailles a verifier (cmp par cmp) :
    do kcmp = 1, nbcmp
        nocmp = zk8(jnocmp-1+kcmp)
        verif_coef_drz = ASTER_FALSE
        ! Exceptions :
        ! E1) PESA_R / ROTA_R sont en général utilisés sans préciser les mailles
        if (nomgd .eq. 'PESA_R') cycle
        if (nomgd .eq. 'ROTA_R') cycle
        if (nomgd .eq. 'CAGNPO_R') cycle
        !
        ! E2) Valeurs fournies par le code d'AFFE_CHAR_MECA
        if (nomgd(1:5) .eq. 'FORC_' .and. nocmp .eq. 'REP') cycle
        if (nomgd(1:5) .eq. 'FORC_' .and. nocmp .eq. 'PLAN') cycle
        if (nomgd .eq. 'VENTCX_F' .and. nocmp .eq. 'FCXP') cycle
        !
        ! E3) Valeurs fournies en loucede par le code d'AFFE_CARA_ELEM
        if (nomgd .eq. 'CAMA_R' .and. nocmp .eq. 'C') cycle
        if (nomgd .eq. 'CACOQU_R' .and. nocmp .eq. 'KAPPA') cycle
        if (nomgd .eq. 'CACOQU_R' .and. nocmp .eq. 'CTOR') verif_coef_drz = ASTER_TRUE
        if (nomgd .eq. 'CAORIE_R' .and. nocmp .eq. 'ALPHA') cycle
        !
        if (nomgd .eq. 'CINFDI_R' .and. nocmp(1:3) .eq. 'REP') cycle
        if (nomgd .eq. 'CINFDI_R' .and. nocmp(1:3) .eq. 'SYM') cycle
        !
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'TSEC') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'HY1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'HZ1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'HY2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'HZ2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'EPY1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'EPY2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'EPZ1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'EPZ2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'EP1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'EP2') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'R1') cycle
        if (nomgd .eq. 'CAGEPO_R' .and. nocmp .eq. 'R2') cycle
        !
        nbmapb = 0
        do kma = 1, nbma_verif
            ima = numa_verif(kma)
            igrel = num_grel(2*(ima-1)+1)
            if ((a_un_sens((igrel-1)*nbcmp+kcmp) .eq. 1) .and. (.not. verif_coef_drz)) cycle
            !
            ient = ptma(ima)
            decal = 3+2*nbgdmx+nec*(ient-1)
            dg => desc(decal+1:decal+nec)
            if (.not. exisdg(dg, kcmp)) cycle
            ! si la cmp est nulle, on n'alarme pas :
            !   comptage des cmps presentes pour pouvoir acceder a la valeur
            ico = 0
            do kcmp2 = 1, kcmp
                if (.not. (exisdg(dg, kcmp2))) cycle
                ico = ico+1
            end do
            iad1 = (ient-1)*nbcmp+ico
            !
            if (tsca .eq. 'R') then
                if (verif_coef_drz) then
                    if (zr(jvale-1+iad1) .lt. 0.d0) then
                        if (typmail(ima) .eq. typt3) then
                            exiq3_coef_drz = exiq3_coef_drz .or. ASTER_TRUE
                        else if (typmail(ima) .eq. typq4) then
                            exiq4_coef_drz = (exiq4_coef_drz .or. ASTER_TRUE)
                            exiq4_drz_nook = (exiq4_coef_drz) .and. &
                                             (zr(jvale-1+iad1) .gt. -1.0d12) .and. &
                                             (zr(jvale-1+iad1) .lt. -1.0d2)
                        end if
                    else
                        cycle
                    end if
                else
                    if (zr(jvale-1+iad1) .eq. 0.d0) cycle
                end if
            else if (tsca .eq. 'C') then
                if (abs(zc(jvale-1+iad1)) .eq. 0.d0) cycle
            else if (tsca .eq. 'I') then
                if (zi(jvale-1+iad1) .eq. 0) cycle
            else if (tsca .eq. 'L') then
                if (.not. zl(jvale-1+iad1)) cycle
            else if (tsca(1:2) .eq. 'K8') then
                if (zk8(jvale-1+iad1) .eq. ' ') cycle
                if (zk8(jvale-1+iad1) .eq. '&FOZERO') cycle
                if (zk8(jvale-1+iad1) .eq. 'GLOBAL') cycle
                if (zk8(jvale-1+iad1) .eq. 'LOCAL') cycle
                if (zk8(jvale-1+iad1) .eq. 'VENT') cycle
                if (zk8(jvale-1+iad1) .eq. 'LOCAL_PR') cycle
            else if (tsca(1:3) .eq. 'K16') then
                if (zk16(jvale-1+iad1) .eq. ' ') cycle
                if (zk16(jvale-1+iad1) .eq. '&FOZERO') cycle
            else if (tsca(1:3) .eq. 'K24') then
                if (zk24(jvale-1+iad1) .eq. ' ') cycle
                if (zk24(jvale-1+iad1) .eq. '&FOZERO') cycle
            else
                ASSERT(ASTER_FALSE)
            end if
            !
            ! s'il y a un problème
            nbmapb = nbmapb+1
            te = num_grel(2*(ima-1)+2)
            if (nbmapb .eq. 1) call jenuno(jexnum('&CATA.TE.NOMTE', te), nomte)
            if (nbmapb .le. 5) list_ma_pb(nbmapb) = ima
        end do
        !
        ! Message d'alarme en cas de probleme :
        if (nbmapb .gt. 0) then
            valk(1) = carte
            valk(2) = comment
            valk(3) = nomgd
            valk(4) = nocmp
            valk(5) = nomte
            if (present(non_lin)) then
                if (exiq3_coef_drz .or. exiq4_coef_drz) then
                    call utmess('F', 'CALCULEL_45', nk=5, valk=valk, si=nbmapb)
                end if
            else if (exiq3_coef_drz .and. exiq4_coef_drz) then
                call utmess('A', 'CALCULEL_42', nk=5, valk=valk, si=nbmapb)
                cycle
            else if (exiq3_coef_drz .and. .not. exiq4_coef_drz) then
                call utmess('F', 'CALCULEL_43', nk=5, valk=valk, si=nbmapb)
            else if (.not. exiq3_coef_drz .and. exiq4_drz_nook) then
                call utmess('A', 'CALCULEL_44', nk=5, valk=valk, si=nbmapb)
                cycle
            else
                if (exiq3_coef_drz .or. exiq4_coef_drz) cycle
                if (a_un_sens((igrel-1)*nbcmp+kcmp) .eq. 1) cycle

                if (nocmp .eq. 'C_METR') then
                    call utmess('F', 'MODELISA10_1')
                else
                    call utmess('A', 'CALCULEL_40', nk=5, valk=valk, si=nbmapb)
                end if
            end if
            do k = 1, min(5, nbmapb)
                valk = ' '
                ima = list_ma_pb(k)
                nommai = int_to_char8(ima)
                valk(1) = nommai

                call list_grma(mailla, ima, 4, lgrma, nbgrma)
                do k1 = 1, min(nbgrma, 3)
                    valk(1+k1) = lgrma(k1)
                end do
                if (nbgrma .gt. 3) valk(5) = '...'
                call utmess('I', 'CALCULEL_41', nk=5, valk=valk)
            end do
        end if
    end do
    AS_DEALLOCATE(vi=numa_verif)
    !
999 continue
    AS_DEALLOCATE(vi=a_un_sens)
    AS_DEALLOCATE(vi=num_grel)
    !
    call jedema()
end subroutine
