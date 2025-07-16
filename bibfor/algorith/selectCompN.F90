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

subroutine selectCompN(chams0, nom_cham, type_cham, nbcmp, nom_cmp, ndim_type)
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
    integer(kind=8) :: jcnsd, nbcmp0, icmp, icmp0, nbcmpD
    character(len=8)    :: nomgd
    character(len=8), pointer  :: cnsk(:) => null()
    character(len=8), pointer  :: cnsc(:) => null()
    aster_logical :: l_cmp_ok

    character(len=8)    :: sief_cmp(6), epsi_cmp(6), depl_cmp(6), forc_cmp(6), flux_cmp(3)
    data sief_cmp/'SIXX', 'SIYY', 'SIZZ', 'SIXY', 'SIXZ', 'SIYZ'/
    data epsi_cmp/'EPXX', 'EPYY', 'EPZZ', 'EPXY', 'EPXZ', 'EPYZ'/
    data depl_cmp/'DX', 'DY', 'DZ', 'DRX', 'DRY', 'DRZ'/
    data forc_cmp/'FX', 'FY', 'FZ', 'MX', 'MY', 'MZ'/
    data flux_cmp/'FLUX', 'FLUY', 'FLUZ'/
! --------------------------------------------------------------------------------------------------
    call jeveuo(chams0//'.CNSK', 'L', vk8=cnsk)
    call jeveuo(chams0//'.CNSD', 'L', jcnsd)
    call jeveuo(chams0//'.CNSC', 'L', vk8=cnsc)
    nomgd = cnsk(2)

!   nombre de composantes présentes
    nbcmp0 = zi(jcnsd-1+2)
!   noms des composantes : cnsc(1:nbcmp)

    ndim_type = 3
    nbcmp = 6
    if (type_cham(6:7) .eq. '2D') then
        ndim_type = 2
        nbcmp = 4
    end if

    if (type_cham(1:5) .eq. 'TENS_') then
        nbcmp = 6
        if (ndim_type .eq. 2) nbcmp = 4
        if (nomgd(1:4) .eq. 'SIEF') then
            nom_cmp(1:nbcmp) = sief_cmp(1:nbcmp)
        elseif (nomgd(1:4) .eq. 'EPSI') then
            nom_cmp(1:nbcmp) = epsi_cmp(1:nbcmp)
        else
            call utmess('F', 'ALGORITH2_5', nk=2, valk=[nomgd, type_cham])
        end if
        nbcmpD = nbcmp
    elseif (type_cham(1:5) .eq. 'VECT_') then
        nbcmp = ndim_type
        nbcmpD = nbcmp
        if (ndim_type .eq. 3) nbcmp = 6
        if (nomgd(1:4) .eq. 'DEPL') then
            nom_cmp(1:nbcmp) = depl_cmp(1:nbcmp)
        elseif (nomgd(1:4) .eq. 'FORC') then
            nom_cmp(1:nbcmp) = forc_cmp(1:nbcmp)
        elseif (nomgd(1:4) .eq. 'FLUX') then
            nbcmp = 3
            nom_cmp(1:nbcmp) = flux_cmp(1:nbcmp)
        else
            call utmess('F', 'ALGORITH2_5', nk=2, valk=[nomgd, type_cham])
        end if
    else
        call utmess('F', 'ALGORITH2_6', sk=type_cham)
    end if

    do icmp = 1, nbcmp
        l_cmp_ok = ASTER_FALSE
        do icmp0 = 1, nbcmp0
            if (cnsc(icmp0) .eq. nom_cmp(icmp)) then
                l_cmp_ok = ASTER_TRUE
                exit
            end if
        end do
        if (.not. l_cmp_ok .and. icmp .le. nbcmpD) then
            call utmess('F', 'ALGORITH2_12', nk=3, &
                        valk=[nom_cham, type_cham, nom_cmp(icmp)])
        elseif (.not. l_cmp_ok) then
            exit
        elseif (l_cmp_ok .and. icmp .eq. (nbcmpD+1)) then
            nbcmpD = nbcmp
        end if
    end do
    nbcmp = nbcmpD

end subroutine selectCompN
