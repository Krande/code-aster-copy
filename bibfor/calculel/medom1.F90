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
!
subroutine medom1(modele, mater, mateco, cara, kcha, nbLoad, &
                  result, nuord)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/rslesd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbLoad, nuord
    character(len=8) :: modele, cara, result
    character(len=24) :: mater, mateco
    character(len=19) :: kcha
!     SAISIE ET VERIFICATION DE LA COHERENCE DES DONNEES MECANIQUES
!     DU PROBLEME
!
! ----------------------------------------------------------------------
! OUT : MODELE : NOM DU MODELE
! OUT : MATER  : CHAMP MATERIAU
! OUT : MATECO : MATERIAU CODE
! OUT : CARA   : NOM DU CHAMP DE CARACTERISTIQUES
! IN  : KCHA   : NOM JEVEUX POUR STOCKER LES CHARGES
! OUT : NCHA   : NOMBRE DE CHARGES
! IN  : RESULT : NOM DE LA SD RESULTAT
! IN  : NUORD  : NUMERO D'ORDRE
! ----------------------------------------------------------------------
    integer(kind=8) :: iret
    integer(kind=8) :: iexcit, i, icha, ie, ikf, in
    integer(kind=8) :: n, n1, n2, n3, n5, nbLoadEff
!-----------------------------------------------------------------------
    character(len=8) :: k8b, nomo, materi, loadType
    character(len=8) :: blan8
    character(len=16) :: concep, nomcmd, phenom
    character(len=19) :: excit
    integer(kind=8), pointer :: infc(:) => null()
    character(len=24), pointer :: fcha(:) => null()
    character(len=24), pointer :: lcha(:) => null()
    aster_logical :: l_ther
!
    call jemarq()
!
    blan8 = ' '
    nbLoad = 0
    modele = ' '
    cara = ' '
    materi = ' '
    nomo = blan8
    iexcit = 1
    n1 = 0
!
    call getres(k8b, concep, nomcmd)
!
    if ((nomcmd .eq. 'CALC_CHAMP') .or. (nomcmd .eq. 'CALC_ERREUR') .or. &
        (nomcmd .eq. 'CALC_META') .or. (nomcmd .eq. 'CALC_G_XFEM') .or. &
        (nomcmd .eq. 'CALC_G')) then
!
!        RECUPERATION DU MODELE, MATERIAU, CARA_ELEM et EXCIT
!        POUR LE NUMERO d'ORDRE NUORD
!
        call rslesd(result, nuord, modele, materi, cara, &
                    excit, iexcit)
        call dismoi('PHENOMENE', modele, 'MODELE', repk=phenom)
        l_ther = ASTER_FALSE
        if (phenom(1:5) .eq. 'THERM') then
            l_ther = ASTER_TRUE
        end if

        if (materi .ne. blan8) then
            call rcmfmc(materi, mateco, l_ther_=l_ther)
        else
            mateco = ' '
        end if
    else
!
        call getvid(' ', 'MODELE', scal=modele, nbret=n1)
        call getvid(' ', 'CARA_ELEM', scal=cara, nbret=n2)
        call dismoi('EXI_RDM', modele, 'MODELE', repk=k8b)
        if ((n2 .eq. 0) .and. (k8b(1:3) .eq. 'OUI')) then
            call utmess('A', 'CALCULEL3_39')
        end if
        call dismoi('PHENOMENE', modele, 'MODELE', repk=phenom)
        l_ther = ASTER_FALSE
        if (phenom(1:5) .eq. 'THERM') then
            l_ther = ASTER_TRUE
        end if
!
        call getvid(' ', 'CHAM_MATER', scal=materi, nbret=n3)
        call dismoi('BESOIN_MATER', modele, 'MODELE', repk=k8b)
        if ((n3 .eq. 0) .and. (k8b(1:3) .eq. 'OUI')) then
            call utmess('A', 'CALCULEL3_40')
        end if
!
        if (n3 .ne. 0) then
            call rcmfmc(materi, mateco, l_ther_=l_ther)
        else
            mateco = ' '
        end if
    end if
!
    mater = materi
!
!     TRAITEMENT DU CHARGEMENT
!
!     SI IEXCIT=1 ON PREND LE CHARGEMENT DONNE PAR L'UTILISATEUR
    if (iexcit .eq. 1) then
        if (getexm('EXCIT', ' ') .eq. 0) then
            n5 = 0
        else
            call getfac('EXCIT', n5)
        end if
!
        if (n5 .ne. 0) then
            nbLoad = n5
            call jeexin(kcha//'.LCHA', iret)
            if (iret .ne. 0) then
                call jedetr(kcha//'.LCHA')
                call jedetr(kcha//'.FCHA')
            end if
            call wkvect(kcha//'.LCHA', 'V V K8', n5, icha)
            call wkvect(kcha//'.FCHA', 'V V K8', n5, ikf)
            do iexcit = 1, n5
                call getvid('EXCIT', 'CHARGE', iocc=iexcit, scal=zk8(icha+iexcit-1), nbret=n)
                call getvid('EXCIT', 'FONC_MULT', iocc=iexcit, scal=k8b, nbret=n)
                if (n .ne. 0) then
                    zk8(ikf+iexcit-1) = k8b
                end if
            end do
        else
            call jeexin(kcha//'.LCHA', iret)
            if (iret .ne. 0) then
                call jedetr(kcha//'.LCHA')
                call jedetr(kcha//'.FCHA')
            end if
            call wkvect(kcha//'.LCHA', 'V V K8', 1, icha)
            call wkvect(kcha//'.FCHA', 'V V K8', 1, ikf)
        end if
!
        if (nbLoad .gt. 0) then
!           VERIFICATION QUE LES CHARGES PORTENT SUR LE MEME MODELE.
            call dismoi('NOM_MODELE', zk8(icha), 'CHARGE', repk=nomo)
            do i = 1, nbLoad
                call dismoi('NOM_MODELE', zk8(icha-1+i), 'CHARGE', repk=k8b)
                if (k8b .ne. nomo) then
                    call utmess('F', 'CALCULEL3_41')
                end if
            end do
!           VERIFICATION QUE LES CHARGES PORTENT SUR LE MODELE
            if (n1 .ne. 0 .and. modele .ne. nomo) then
                call utmess('F', 'CALCULEL3_42')
            end if
        end if
!
!     SI IEXCIT=0 ON PREND LE CHARGEMENT PRESENT DANS LA SD
!
    else
!
        call jeveuo(excit//'.INFC', 'L', vi=infc)
        call jeveuo(excit//'.LCHA', 'L', vk24=lcha)
        call jeveuo(excit//'.FCHA', 'L', vk24=fcha)
        nbLoad = infc(1)

        if (nbLoad .eq. 0) then
            nbLoadEff = 1
        else
            nbLoadEff = nbLoad
        end if
!
        call jeexin(kcha//'.LCHA', iret)
        if (iret .ne. 0) then
            call jedetr(kcha//'.LCHA')
            call jedetr(kcha//'.FCHA')
        end if
        call wkvect(kcha//'.LCHA', 'V V K8', nbLoadEff, icha)
        call wkvect(kcha//'.FCHA', 'V V K8', nbLoadEff, ikf)
        call dismoi('PHENOMENE', modele, 'MODELE', repk=phenom)
        in = 0
        do i = 1, nbLoad
            call dismoi('TYPE_CHARGE', lcha(i), 'CHARGE', repk=loadType, arret='C', ier=ie)
            if ((ie .eq. 0) .and. (phenom(1:4) .eq. loadType(1:4))) then
                zk8(icha+in) = lcha(i) (1:8)
                zk8(ikf+in) = fcha(i) (1:8)
                in = in+1
            end if
        end do
        nbLoad = in
!
    end if
!
    call jedema()
end subroutine
