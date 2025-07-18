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
subroutine ccvrpu(resuin, lisord, nbordr)
    implicit none
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbordr
    character(len=8) :: resuin
    character(len=19) :: lisord
! person_in_charge: nicolas.sellenet at edf.fr
! ----------------------------------------------------------------------
!  CALC_CHAMP - VERIFICATION DES PARAMETRES UTILISATEURS :
!  -    -       - -              -          -
!    MODELE, CARA_ELEM, ...
! ----------------------------------------------------------------------
!
!  LE BUT DE CETTE ROUTINE EST DE VERIFIER QUE L'UTILISATEUR NE
!    SURCHARGE PAS LE MODELE, CARA_ELEM, CHAM_MATER OU LE CHARGEMENT
!  SI UN DE CES ELEMENTS EST PRESENT DANS LA SD ET QUE L'UTILISATEUR
!    LE DONNE EN ENTREE DE CALC_CHAMP, ON VERIFIE QUE C'EST LE MEME
!    SINON ON INTERDIT LA REENTRANCE
!
! IN  :
!   RESUIN K8   NOM DE LA SD IN
!   LISORD K19  NOM DE LA LISTE DES NUMEROS D'ORDRE
!   NBORDR I    NOMBRE DE NUMEROS D'ORDRE
! ----------------------------------------------------------------------
!
    integer(kind=8) :: jordr, iordr, numord, jpara, n1, n2, n3, nchalu, icharg
    integer(kind=8) :: lchalu, fchalu, nchasd, jfcha, ilu, isd
!
    character(len=8) :: k8b, modelu, carelu, chmatu, modelr, carelr, chmatr
    character(len=8) :: fonclu
    character(len=16) :: valk(3)
    character(len=19) :: kcha, kfon, excit
    character(len=24) :: excisd
    integer(kind=8), pointer :: infc(:) => null()
    character(len=24), pointer :: lcha(:) => null()
!
    kcha = '&&CCVRPU.CHARGE    '
    kfon = '&&CCVRPU.FONC_MULT '
!
    call jeveuo(lisord, 'L', jordr)
!
    modelu = ' '
    carelu = ' '
    chmatu = ' '
    call getvid(' ', 'MODELE', scal=modelu, nbret=n1)
    call getvid(' ', 'CARA_ELEM', scal=carelu, nbret=n2)
    call getvid(' ', 'CHAM_MATER', scal=chmatu, nbret=n3)
!
    nchalu = 0
    if (getexm('EXCIT', 'CHARGE') .eq. 1) then
        call getfac('EXCIT', nchalu)
!
        if (nchalu .ne. 0) then
            call wkvect(kcha, 'V V K8', nchalu, lchalu)
            call wkvect(kfon, 'V V K8', nchalu, fchalu)
!
            do icharg = 1, nchalu
                call getvid('EXCIT', 'CHARGE', iocc=icharg, scal=zk8(lchalu+icharg-1), nbret=n1)
!
                call getvid('EXCIT', 'FONC_MULT', iocc=icharg, scal=fonclu, nbret=n2)
!
                if (n2 .ne. 0) then
                    zk8(fchalu+icharg-1) = fonclu
                end if
            end do
        end if
    end if
!
    if (modelu .ne. ' ' .or. carelu .ne. ' ' .or. chmatu .ne. ' ' .or. nchalu .ne. 0) then
        do iordr = 1, nbordr
            numord = zi(jordr-1+iordr)
!
!         VERIFICATION DU MODELE
            if (modelu .ne. ' ') then
                call rsadpa(resuin, 'L', 1, 'MODELE', numord, &
                            0, sjv=jpara, styp=k8b)
                modelr = zk8(jpara)
                if (modelr .ne. ' ' .and. modelr .ne. modelu) then
                    valk(1) = 'MODELE'
                    valk(2) = modelr
                    valk(3) = modelu
                    call utmess('F', 'CALCULEL_33', nk=3, valk=valk)
                    ASSERT(.false.)
                end if
            end if
!
!         VERIFICATION DU CARAELEM
            if (carelu .ne. ' ') then
                call rsadpa(resuin, 'L', 1, 'CARAELEM', numord, &
                            0, sjv=jpara, styp=k8b)
                carelr = zk8(jpara)
                if (carelr .ne. ' ' .and. carelr .ne. carelu) then
                    valk(1) = 'CARA_ELEM'
                    valk(2) = carelr
                    valk(3) = carelu
                    call utmess('F', 'CALCULEL_33', nk=3, valk=valk)
                    ASSERT(.false.)
                end if
            end if
!
!         VERIFICATION DU CHAMATER
            if (chmatu .ne. ' ') then
                call rsadpa(resuin, 'L', 1, 'CHAMPMAT', numord, &
                            0, sjv=jpara, styp=k8b)
                chmatr = zk8(jpara)
                if (chmatr .ne. ' ' .and. chmatr .ne. chmatu) then
                    valk(1) = 'CHAM_MATER'
                    valk(2) = chmatr
                    valk(3) = chmatu
                    call utmess('F', 'CALCULEL_33', nk=3, valk=valk)
                    ASSERT(.false.)
                end if
            end if
!
!         VERIFICATION DU CHARGEMENT
            if (nchalu .ne. 0) then
                call rsadpa(resuin, 'L', 1, 'EXCIT', numord, &
                            0, sjv=jpara, styp=k8b)
                excisd = zk24(jpara)
                nchasd = 0
                if (excisd .ne. ' ') then
                    excit = excisd(1:19)
                    call jeveuo(excit//'.LCHA', 'L', vk24=lcha)
                    call jeveuo(excit//'.INFC', 'L', vi=infc)
                    call jeveuo(excit//'.FCHA', 'L', jfcha)
                    nchasd = infc(1)
                    if (nchasd .ne. nchalu) then
                        call utmess('F', 'CALCULEL_39')
                        ASSERT(.false.)
                    end if
                    do ilu = 1, nchalu
                        do isd = 1, nchasd
                            if (zk8(lchalu-1+ilu) .eq. lcha(isd) (1:8)) goto 30
                        end do
                        call utmess('F', 'CALCULEL_39')
30                      continue
                    end do
                end if
            end if
        end do
    end if
!
    call jedetr(kcha)
    call jedetr(kfon)
!
end subroutine
