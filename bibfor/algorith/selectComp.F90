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

subroutine selectComp(chams0, nom_cham, type_cham, nbcmp, nom_cmp, ndim_type)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    character(len=19), intent(in) :: chams0
    character(len=*), intent(in) :: nom_cham, type_cham
    integer(kind=8), intent(out) :: nbcmp, ndim_type
    character(len=8), intent(out) :: nom_cmp(*)
!
! --------------------------------------------------------------------------------------------------
!
!                       MODI_REPERE : SELECTION DES COMPOSANTES A TRAITER
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jcesd, nbcmp0, icmp, icmp0
    character(len=8)    :: nomgd
    character(len=8), pointer  :: cesk(:) => null()
    character(len=8), pointer  :: cesc(:) => null()
    aster_logical :: l_cmp_ok

    character(len=8)    :: sief_cmp(6), epsi_cmp(6), coque_efge(8), coque_dege(8)
    character(len=8)    :: pout_effo(6)
    data sief_cmp/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
    data epsi_cmp/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
    data coque_efge/'NXX', 'NYY', 'NXY', 'MXX', 'MYY', 'MXY', 'QX', 'QY'/
    data coque_dege/'EXX', 'EYY', 'EXY', 'KXX', 'KYY', 'KXY', 'GAX', 'GAY'/
    data pout_effo/'N', 'VY', 'VZ', 'MT', 'MFY', 'MFZ'/
! --------------------------------------------------------------------------------------------------
    call jeveuo(chams0//'.CESK', 'L', vk8=cesk)
    call jeveuo(chams0//'.CESD', 'L', jcesd)
    call jeveuo(chams0//'.CESC', 'L', vk8=cesc)
    nomgd = cesk(2)

!   nombre de composantes présentes
    nbcmp0 = zi(jcesd-1+2)
!   noms des composantes : cesc(1:nbcmp)

    ndim_type = 3
    nbcmp = 6

    if (type_cham(1:5) .eq. 'TENS_') then
        if (type_cham(6:7) .eq. '2D') then
            ndim_type = 2
            nbcmp = 4
        end if
        if (nomgd(1:4) .eq. 'SIEF') then
            nom_cmp(1:nbcmp) = sief_cmp(1:nbcmp)
        elseif (nomgd(1:4) .eq. 'EPSI') then
            nom_cmp(1:nbcmp) = epsi_cmp(1:nbcmp)
        else
            call utmess('F', 'ALGORITH2_5', nk=2, valk=[nomgd, type_cham])
        end if
    elseif (type_cham(1:10) .eq. 'COQUE_GENE') then
        nbcmp = 8
        if (nomgd(1:4) .eq. 'SIEF') then
            nom_cmp(1:nbcmp) = coque_efge(1:nbcmp)
        elseif (nomgd(1:4) .eq. 'EPSI') then
            nom_cmp(1:nbcmp) = coque_dege(1:nbcmp)
        else
            call utmess('F', 'ALGORITH2_5', nk=2, valk=[nomgd, type_cham])
        end if
    elseif (type_cham(1:7) .eq. '1D_GENE') then
        if (nomgd(1:4) .eq. 'SIEF') then
            nom_cmp(1:nbcmp) = pout_effo(1:nbcmp)
        else
            call utmess('F', 'ALGORITH2_5', nk=2, valk=[nomgd, type_cham])
        end if
    else
        call utmess('F', 'ALGORITH2_14', sk=type_cham)
    end if

    do icmp = 1, nbcmp
        l_cmp_ok = ASTER_FALSE
        do icmp0 = 1, nbcmp0
            if (cesc(icmp0) .eq. nom_cmp(icmp)) then
                l_cmp_ok = ASTER_TRUE
                exit
            end if
        end do
        if (.not. l_cmp_ok) call utmess('F', 'ALGORITH2_12', nk=3, &
                                        valk=[nom_cham, type_cham, nom_cmp(icmp)])
    end do

end subroutine selectComp
