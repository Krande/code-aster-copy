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
subroutine cbchei(load, mesh, model, valeType)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getfac.h"
#include "asterfort/cachei.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/alcart.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of load PRE_EPSI
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : mesh
! In  load             : load
! In  model            : model
! In  geomDime         : space dimension
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: keywordFact = 'PRE_EPSI'
    integer(kind=8) :: nbfac, iepsi, ncmp, cret, jcesd, jcesc, nbcmpch, i, j
    character(len=4) :: typch
    character(len=5) :: para
    character(len=19) :: carte, chames
    character(len=24) :: chepsi
    aster_logical :: compok
    character(len=8), pointer :: valv(:) => null()
    character(len=8), pointer :: vncmp(:) => null()
    integer(kind=8), parameter :: nbcmpdisp = 17
    character(len=8), parameter :: nomcmpdisp(nbcmpdisp) = (/ &
                                   'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ', &
                                   'EPX ', 'KY  ', 'KZ  ', 'EXX ', 'EYY ', 'EXY ', &
                                   'KXX ', 'KYY ', 'KXY ', &
                                   'GAX ', 'GAY '/)
!
! --------------------------------------------------------------------------------------------------
!
    call getfac(keywordFact, nbfac)
!
    if (nbfac .ne. 0) then
        para = 'EPSIN'

        iepsi = 0
        if (valeType .eq. 'REEL') then
            call getvid(keywordFact, 'EPSI', iocc=1, scal=chepsi, nbret=iepsi)
        end if

        if (iepsi .eq. 0) then
            call cachei(load, model, mesh, valeType, para, keywordFact)
        else
            if (nbfac .gt. 1) call utmess('F', 'CHARGES_5')

!
! ---       verification des composantes
!
            chames = '&&CHCHEI.CES'

            call dismoi('TYPE_CHAMP', chepsi, 'CHAMP', repk=typch)

            if (typch .eq. 'CART') then
                call carces(chepsi, 'ELEM', ' ', 'V', chames, &
                            ' ', cret)
            elseif (typch .eq. 'ELGA') then
                call celces(chepsi, 'V', chames)
            else
                call utmess('F', 'CHARGES_8', sk=typch)
            end if

            call jeveuo(chames//'.CESD', 'L', jcesd)
            call jeveuo(chames//'.CESC', 'L', jcesc)

            nbcmpch = zi(jcesd+1)
            do i = 1, nbcmpch
                compok = ASTER_FALSE
                do j = 1, nbcmpdisp
                    if (zk8(jcesc-1+i) .eq. nomcmpdisp(j)) then
                        compok = ASTER_TRUE
                        exit
                    end if
                end do
                if (.not. compok) then
                    call utmess('F', 'CHARGES_6', sk=zk8(jcesc-1+i))
                end if
            end do

            call jedetr(chames)
!
            carte = load//'.CHME.'//para
            call alcart('G', carte, mesh, 'NEUT_K8')
            call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
            call jeveuo(carte//'.VALV', 'E', vk8=valv)
!
            ncmp = 1
            vncmp(1) = 'Z1'
            valv(1) = chepsi(1:8)
            call nocart(carte, 1, ncmp)
        end if
    end if
!
end subroutine
