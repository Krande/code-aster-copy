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

subroutine ef0517(nomte)
!
!
! --------------------------------------------------------------------------------------------------
!
!                   EFGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=16) :: nomte
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jsd1ff.h"
#include "asterfort/lonele.h"
#include "asterfort/matela.h"
#include "asterfort/moytem.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/pmfitg.h"
#include "asterfort/pmfmats.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/jevech.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kp, ncomp, ii, jacf
    integer(kind=8) :: icgp, icontn, npg, istrxr, nc
    real(kind=8) :: flel(14), co(3), xl, d1b(7, 14), carsec(6)
!
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10)
!
    integer(kind=8) :: npge, kk, iret, ipoids, imate, ifgp
    real(kind=8) :: aa, alfay, alfaz, young, nu, gcis, ey, ez, xiy, xiz, temp, phiy, phiz
    character(len=8) :: mator
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: nb_cara = 4
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8), parameter :: noms_cara(nb_cara) = (/'AY1  ', 'AZ1  ', 'EY1  ', 'EZ1  '/)
!
! --------------------------------------------------------------------------------------------------
!&<
    if (nomte .eq. 'MECA_POU_D_TGM') then
        call jevech('PMATERC', 'L', imate)
        call jevech('PCONTRR', 'L', icgp)
        call jevech('PSTRXRR', 'L', istrxr)
        call jevech('PEFFORR', 'E', icontn)
        ! Récupération des caractéristiques des fibres
        call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, jacf=jacf)
        !
        nc = 7; npg = 3; ncomp = 18
        !
        call elrefe_info(fami='RIGI', npg=npge, jpoids=ipoids)
        ASSERT(npg .ge. npge)
        co(1:npg) = zr(ipoids:ipoids+npg-1)
        ! Longueur de l'élément
        xl = lonele()
        !  coefficient dependant de la temperature moyenne
        call moytem('RIGI', npg, 1, '+', temp, iret)
        call pmfmats(mator)
        ASSERT(mator .ne. ' ')
        call matela(zi(imate), mator, 1, temp, young, nu)
        gcis = young/(2.d0*(1.d0+nu))
        !
        call poutre_modloc('CAGNP2', noms_cara, nb_cara, lvaleur=vale_cara)
        alfay = vale_cara(1)
        alfaz = vale_cara(2)
        ey    = vale_cara(3)
        ez    = vale_cara(4)
        ! Récupération des caracteristiques de la section
        call pmfitg(tygrfi, nbfibr, nbcarm, zr(jacf), carsec)
        aa  = carsec(1)
        xiz = carsec(4)
        xiy = carsec(5)
        !
        phiy = young*xiz*12.d0*alfay/(xl*xl*gcis*aa)
        phiz = young*xiy*12.d0*alfaz/(xl*xl*gcis*aa)
        !
        flel(:) = 0.0d+0
        ! boucle sur les points de gauss
        do kp = 1, npg
            call jsd1ff(kp, xl, phiy, phiz, d1b)
            ifgp = ncomp*(kp-1)-1
            do kk = 1, 2*nc
                do ii = 1, nc
                    flel(kk) = flel(kk)+xl*0.5*zr(istrxr+ifgp+ii)*d1b(ii,kk)*co(kp)
                end do
            end do
        end do
        ! Prise en compte du centre de torsion
        flel(4)  = flel(4)-ez*flel(2)+ey*flel(3)
        flel(11) = flel(11)-ez*flel(9)+ey*flel(10)
        ! Comme d'hab on change le signe sur le noeud 1
        do ii = 1, nc
            zr(icontn+ii-1) = -flel(ii)
        end do
        do ii = nc+1, 2*nc
            zr(icontn+ii-1) = flel(ii)
        end do
!&>
! --------------------------------------------------------------------------------------------------
!
    else if (nomte .eq. 'MECA_POU_D_EM') then
        nc = 6
        npg = 2
        ncomp = 18
        call jevech('PSTRXRR', 'L', istrxr)
        call jevech('PEFFORR', 'E', icontn)
        kp = 1
        do ii = 1, nc
            zr(icontn-1+nc*(kp-1)+ii) = -zr(istrxr-1+ncomp*(kp-1)+ii)
        end do
        kp = 2
        do ii = 1, nc
            zr(icontn-1+nc*(kp-1)+ii) = zr(istrxr-1+ncomp*(kp-1)+ii)
        end do
!
! --------------------------------------------------------------------------------------------------
!
    else if (nomte .eq. 'MECA_POU_D_SQUE') then
        nc = 9
        npg = 2
        ncomp = 21
        call jevech('PSTRXRR', 'L', istrxr)
        call jevech('PEFFORR', 'E', icontn)
        kp = 1
        do ii = 1, 6
            zr(icontn-1+nc*(kp-1)+ii) = -zr(istrxr-1+ncomp*(kp-1)+ii)
        end do
        do ii = 7, nc
            zr(icontn-1+nc*(kp-1)+ii) = -zr(istrxr-1+ncomp*(kp-1)+12+ii)
        end do
        kp = 2
        do ii = 1, 6
            zr(icontn-1+nc*(kp-1)+ii) = zr(istrxr-1+ncomp*(kp-1)+ii)
        end do
        do ii = 7, nc
            zr(icontn-1+nc*(kp-1)+ii) = zr(istrxr-1+ncomp*(kp-1)+12+ii)
        end do

    else
        ASSERT(.false.)
    end if
!
end subroutine
