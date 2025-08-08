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

subroutine rc3600_chtotab(nomtb, conceptin, nsymb, champ)
!
! --------------------------------------------------------------------------------------------------
!
!                       OPERATEUR POST_RCCM , MOMENT EQUIVALENT
!
!     CONSTRUCTION DE LA TABLE A PARTIR D'UN CHAMP
!
!        IN     : conceptin (K8)  : objet en entrée de POST_RCCM
!                 nsymb     (K16) : si RESULTAT en entrée : nom du symbolique du champ
!                                 : sinon ' '
!                 MESMAI (K24) : OBJET DES NOMS DE MAILLE
!                 NBMA  (I)    : NOMBRE DE MAILLES UTILISATEUR
!                 CHAMP : nom du champ à traiter
!
!        IN/OUT : NOMTB (K24)  : OBJET TABLE
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
!
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/indiis.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/reliem1.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbcrsv.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: conceptin, nomtb
    character(len=16) :: nsymb
    character(len=19) :: champ
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbcmp = 4
    integer(kind=8) :: nbma, nbpara, jlma, iret, jcesl, jcesl2, jcesd, jcesd2, nbpara_tmp
    integer(kind=8) :: jconx2, nbcmpx, nr, nk, kk, ima, nbmax, indma, nbpt, nbcmpt
    integer(kind=8) :: nbabsc, icmp, iad, ipt, inot, ispt, kcp, iexi
    real(kind=8) :: table_valr(5), val_cmp(nbcmp), val_absc(4)
    complex(kind=8) :: cbid
    character(len=8) :: motcle(3), typmcl(3), noma, kma, kno
    character(len=16) :: table_valk(4)
    character(len=19) :: chames, chamescurv, ligrelCham
    character(len=24) :: mesmai
    aster_logical :: l_abscurv, l_alarm

    real(kind=8), pointer :: cesv(:) => null()
    character(len=8), pointer :: cesc(:) => null()
    real(kind=8), pointer :: cesv2(:) => null()
    integer(kind=8), pointer :: connex(:) => null()

    character(len=9) :: parata1(9), parata2(8), parata(9), parata_tmp(9)
    character(len=3) :: typarata1(9), typarata2(8), typarata(9), nomcmp(nbcmp)
    data nomcmp/'MT ', 'MFY', 'MFZ', 'MEQ'/
    data parata1/'RESULTAT ', 'NOM_CHAM ', 'MAILLE   ', 'NOEUD    ', 'ABSC_CURV', &
        'MT       ', 'MFY      ', 'MFZ      ', 'MEQ      '/
    data parata2/'CHAM_GD  ', 'MAILLE   ', 'NOEUD    ', 'ABSC_CURV', &
        'MT       ', 'MFY      ', 'MFZ      ', 'MEQ      '/
    data typarata1/'K8 ', 'K16', 'K8 ', 'K8 ', 'R  ', 'R  ', 'R  ', 'R  ', 'R  '/
    data typarata2/'K8 ', 'K8 ', 'K8 ', 'R  ', 'R  ', 'R  ', 'R  ', 'R  '/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

!   0 - Initialisations
    cbid = (0.d0, 0.d0)

!   1 - Création de la table

    if (nsymb .eq. ' ') then
        nbpara = 8
        parata(1:8) = parata2(1:8)
        typarata(1:8) = typarata2(1:8)
        nr = 5
        nk = 3
        table_valk(1) = conceptin
        kk = 1
    else
        nbpara = 9
        parata(1:9) = parata1(1:9)
        typarata(1:9) = typarata1(1:9)
        nr = 5
        nk = 4
        table_valk(1) = conceptin
        table_valk(2) = nsymb
        kk = 2
    end if
    call tbcrsv(nomtb, 'G', nbpara, parata, typarata, 0)

!   2 - Liste des mailles

    mesmai = '&&RC3600.MES_MAILLES'
    motcle(1) = 'TOUT'
    motcle(2) = 'MAILLE'
    motcle(3) = 'GROUP_MA'
    typmcl(1) = 'TOUT'
    typmcl(2) = 'MAILLE'
    typmcl(3) = 'GROUP_MA'
    call dismoi('NOM_LIGREL', champ, 'CHAMP', repk=ligrelCham)
    call dismoi('NOM_MAILLA', ligrelCham, 'LIGREL', repk=noma)
    call reliem1(ligrelCham, noma, 'NU_MAILLE', 'ZONE_ANALYSE', 1, &
                 3, motcle, typmcl, mesmai, nbma)
    call jeveuo(mesmai, 'L', jlma)

!   3 - Abscisses curvilignes

    call jeexin(noma//'.ABSC_CURV .VALE', iexi)
    if (iexi .eq. 0) then
        l_abscurv = .false.
        call utmess('A', 'POSTRELE_16')
    else
        l_abscurv = .true.
        chamescurv = '&&RC3600.ACURV.CES'
        call carces(noma//'.ABSC_CURV', 'ELEM', ' ', 'V', chamescurv, ' ', iret)
        ASSERT(iret .eq. 0)
        call jeveuo(chamescurv//'.CESV', 'L', vr=cesv2)
        call jeveuo(chamescurv//'.CESL', 'L', jcesl2)
        call jeveuo(chamescurv//'.CESD', 'L', jcesd2)
    end if

!   4 - Champ à traiter

    call jeveuo(noma//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)

    chames = '&&RC3600.CES       '
    call celces(champ, 'V', chames)

    call jeveuo(chames//'.CESV', 'L', vr=cesv)
    call jeveuo(chames//'.CESL', 'L', jcesl)
    call jeveuo(chames//'.CESD', 'L', jcesd)
    call jeveuo(chames//'.CESC', 'L', vk8=cesc)
!
!   NOMBRE DE MAILLES MAX DU CHAMP : NBMAX
    nbmax = zi(jcesd)
!   NOMBRE DE COMPOSANTES MAX DU CHAMP : NBCMPX
    nbcmpx = zi(jcesd+1)
    ASSERT(nbcmpx .eq. nbcmp)

!   5 - parcours de mailles et remplissage de la table
    do ima = 1, nbmax
!       SI LA MAILLE FAIT PARTIE DES MAILLES DESIREES,
!       ON POURSUIT, SINON ON VA A LA MAILLE SUIVANTE:
        indma = indiis(zi(jlma), ima, 1, nbma)
        if (indma .eq. 0) cycle

        kma = int_to_char8(ima)
        table_valk(kk+1) = kma

!       NOMBRE DE POINTS DE LA MAILLE IMA : NBPT
        nbpt = zi(jcesd+5+4*(ima-1))
!       NOMBRE DE COMPOSANTES PORTEES PAR LES POINTS DE LA MAILLE IMA
        nbcmpt = zi(jcesd+5+4*(ima-1)+2)
        ASSERT(nbcmpt .eq. nbcmp)

!       ON PARCOURT LES COMPOSANTES DE LA CARTE DES ABSCISSES CURVILIGNES
        if (l_abscurv) then
            nbabsc = 0
            val_absc(:) = -1.d0
            l_alarm = .false.
            do icmp = 1, nbcmp
                !           Valeur de icmp au point ipt de la maille ima: zr(jcesv2+iad-1)
                call cesexi('C', jcesd2, jcesl2, ima, 1, 1, icmp, iad)
                if (iad .gt. 0) then
                    val_absc(icmp) = cesv2(iad)
                    nbabsc = nbabsc+1
                else
                    if ((.not. l_alarm) .and. icmp .le. nbpt) then
                        call utmess('A', 'POSTRELE_4', sk=kma)
                        l_alarm = .true.
                    end if
                end if
            end do
            ASSERT(nbabsc .eq. 0 .or. nbabsc .eq. nbpt)
        end if
!       la composante i de val_absc correspond au point ipt
!
!       ON PARCOURT LES POINTS DE LA MAILLE IMA
        do ipt = 1, nbpt
!           NUMERO DU POINT (DU MAILLAGE GLOBAL): INOT
            inot = connex(zi(jconx2-1+ima)+ipt-1)
            kno = int_to_char8(inot)
            table_valk(kk+2) = kno
            ispt = 1
            kcp = 0
!           ON PARCOURT LES COMPOSANTES PORTEES PAR LE POINT IPT
            do icmp = 1, nbcmpt
!               Valeur de icmp au point ipt de la maille ima: zr(jcesv+iad-1)
                call cesexi('C', jcesd, jcesl, ima, ipt, ispt, icmp, iad)
                if (iad .gt. 0) then
                    kcp = kcp+1
                    val_cmp(kcp) = cesv(iad)
                    ASSERT(cesc(icmp) .eq. nomcmp(icmp))
                end if
            end do

            if (l_abscurv .and. val_absc(ipt) .gt. -1.d0) then
                table_valr(1) = val_absc(ipt)
                table_valr(2:nbcmp+1) = val_cmp(1:nbcmp)
                parata_tmp(:) = parata(:)
                nbpara_tmp = nbpara
            else
                table_valr(1:nbcmp) = val_cmp(1:nbcmp)
                parata_tmp(1:nbpara-nbcmp-1) = parata(1:nbpara-nbcmp-1)
                parata_tmp(nbpara-nbcmp:nbpara-1) = parata(nbpara-nbcmp+1:nbpara)
                nbpara_tmp = nbpara-1
            end if

            call tbajli(nomtb, nbpara_tmp, parata_tmp, [0], table_valr, &
                        [cbid], table_valk, 0)
        end do
    end do
!
    call jedema()
!
end subroutine
