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

subroutine ef0347(nomte)
!
! --------------------------------------------------------------------------------------------------
!
!     OPTION EFGE_ELNO
!
!     ELEMENTS :
!        POU_D_TG
!        POU_D_T
!        POU_D_E
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/moytem.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rcvalb.h"
!
    character(len=16) :: nomte
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nc, ii, npg, ipoids, imate, k, kk
    integer(kind=8) :: icgp, icontn, kp
!
    real(kind=8) :: d1b(7, 14), co(3), fs(14), xl
    real(kind=8) :: aa, alfay, alfaz, young, nu, gcis, ey, ez, xiy, xiz, temp, phiy, phiz
!
    integer(kind=8)             :: iret(2)
    character(len=2)    :: nomres(2)
    real(kind=8)        :: valres(2)
!
    integer(kind=8), parameter  :: nb_cara = 7
    real(kind=8)        :: vale_cara(nb_cara)
    character(len=8)    :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'EY1', 'EZ1'/
!
    aster_logical :: okelem
!
! --------------------------------------------------------------------------------------------------
!
    okelem = (nomte .eq. 'MECA_POU_D_TG') .or. &
             (nomte .eq. 'MECA_POU_D_T') .or. &
             (nomte .eq. 'MECA_POU_D_E')
    ASSERT(okelem)
!
    call elrefe_info(fami='RIGI', npg=npg, jpoids=ipoids)
    ASSERT((npg .eq. 2) .or. (npg .eq. 3))
!
    call jevech('PCONTRR', 'L', icgp)
    call jevech('PEFFORR', 'E', icontn)

! --------------------------------------------------------------------------------------------------
    if (nomte .eq. 'MECA_POU_D_E' .or. nomte .eq. 'MECA_POU_D_T') then
        ! Pour les POU_D_E et POU_D_T : Même chose que dans te0347
        nc = 6
        if (npg .eq. 2) then
            do ii = 1, nc
                zr(icontn-1+ii) = zr(icgp-1+ii)
                zr(icontn-1+ii+nc) = zr(icgp-1+ii+nc)
            end do
        else
            do ii = 1, nc
                zr(icontn-1+ii) = zr(icgp-1+ii)
                zr(icontn-1+ii+nc) = zr(icgp-1+ii+nc+nc)
            end do
        end if
!
! --------------------------------------------------------------------------------------------------
!&<
    else if (nomte .eq. 'MECA_POU_D_TG') then
        ! Pour les POU_D_TG : Même chose que dans te0347
        call jevech('PMATERC', 'L', imate)
        !
        nc = 7
        co(1:npg) = zr(ipoids:ipoids+npg-1)
        ! Longueur de l'élément
        xl = lonele()
        !  coefficient dependant de la temperature moyenne
        call moytem('RIGI', npg, 1, '+', temp, iret(1))
        nomres(1) = 'E'
        nomres(2) = 'NU'
        call rcvalb('RIGI', 1, 1, '+', zi(imate), ' ', 'ELAS', 1, 'TEMP', [temp], &
                    2, nomres, valres, iret, 1)
        young = valres(1)
        nu    = valres(2)
        gcis  = young/(2.d0*(1.d0+nu))
        !
        call poutre_modloc('CAGNPO', noms_cara, nb_cara, lvaleur=vale_cara)
        aa      = vale_cara(1)
        xiy     = vale_cara(2)
        xiz     = vale_cara(3)
        alfay   = vale_cara(4)
        alfaz   = vale_cara(5)
        ey      = vale_cara(6)
        ez      = vale_cara(7)
        !
        phiy = young*xiz*12.d0*alfay/(xl*xl*gcis*aa)
        phiz = young*xiy*12.d0*alfaz/(xl*xl*gcis*aa)
        !
        fs(:) = 0.0
        do kp = 1, npg
            call jsd1ff(kp, xl, phiy, phiz, d1b)
            do k = 1, 2*nc
                do kk = 1, nc
                    fs(k) = fs(k)+xl*zr(icgp+nc*(kp-1)+kk-1)*d1b(kk,k)*co(kp)*0.50d0
                end do
            end do
        end do
        ! Prise en compte du centre de torsion
        fs(4)  = fs(4)-ez*fs(2)+ey*fs(3)
        fs(11) = fs(11)-ez*fs(9)+ey*fs(10)
        ! Comme d'hab on change le signe sur le noeud 1
        do ii = 1, nc
            zr(icontn+ii-1) = -fs(ii)
        end do
        do ii = nc+1, 2*nc
            zr(icontn+ii-1) = fs(ii)
        end do
!&>
! --------------------------------------------------------------------------------------------------
    end if
!
end subroutine
