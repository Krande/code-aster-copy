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

subroutine rcdiff(imate, comp, temp, c, diff, difl_, difv_)
    implicit none
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rftDiffusion.h"
#include "asterfort/utmess.h"
    integer(kind=8), intent(in) :: imate
    real(kind=8), intent(in) :: temp, c
    character(len=16), intent(in) :: comp
    real(kind=8), intent(out) :: diff
    real(kind=8), intent(out), optional :: difl_, difv_
! ----------------------------------------------------------------------
!     CALCUL DU COEFFICIENT DE DIFFUSION POUR LES LOI DE TYPE SECHAGE
!
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  COMP    : COMPORTEMENT
! IN  TEMP    : TEMPERATURE
! IN  C       : CONCENTRATION EN EAU
! OUT DIFF    : VALEUR DU COEFFICIENT DE DIFFUSION
! OUT DIFL_   : VALEUR DU COEFFICIENT DE DIFFUSION LIQUIDE (OPTIONNEL)
! OUT DIFV_   : VALEUR DU COEFFICIENT DE DIFFUSION VAPEUR (OPTIONNEL)
! ----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nbres
    real(kind=8) :: rap
!-----------------------------------------------------------------------
    parameter(nbres=10)
    integer(kind=8) :: nbpar, kpg, spt
    real(kind=8) :: valres(nbres), valpar(2), tz0
    integer(kind=8) :: icodre(nbres)
    character(len=8) :: nompar(2), fami, poum
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
    real(kind=8) :: val_non_physique
    real(kind=8) :: difl, difv
!
!
    call rccoma(imate, comp(1:6), 1, phenom, icodre(1))
!
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    tz0 = r8t0()
!
    difl = 0.d0
    difv = 0.d0
    if (phenom .eq. 'SECH_GRANGER') then
        nbpar = 0

        nomres(1) = 'A'
        nomres(2) = 'B'
        nomres(3) = 'QSR_K'
        nomres(4) = 'TEMP_0_C'
        call rcvalb(fami, kpg, spt, poum, imate, &
                    ' ', phenom, nbpar, nompar, valpar, &
                    4, nomres, valres, icodre, 1)

        val_non_physique = max(valres(2)*c, -valres(3)*(1.d&
               &0/(temp+tz0)-1.d0/(valres(4)+tz0)))
        if (val_non_physique .gt. 1.d10) then
            call utmess('F', 'ALGORITH10_91', sk=phenom, sr=val_non_physique)
        end if

        difl = valres(1)*exp(valres(2)*c)*((temp+tz0)/(valres(4)+tz0))*exp(-valres(3)*(1.d&
               &0/(temp+tz0)-1.d0/(valres(4)+tz0)))
!
    else if (phenom .eq. 'SECH_MENSI') then
        nbpar = 0
        nomres(1) = 'A'
        nomres(2) = 'B'
        call rcvalb(fami, kpg, spt, poum, imate, &
                    ' ', phenom, nbpar, nompar, valpar, &
                    2, nomres, valres, icodre, 1)
        difl = valres(1)*exp(valres(2)*c)
!
    else if (phenom .eq. 'SECH_BAZANT') then
        nbpar = 1
        nompar(1) = 'SECH'
        valpar(1) = c
        nomres(1) = 'D1'
        nomres(2) = 'ALPHA_BAZANT'
        nomres(3) = 'N'
        nomres(4) = 'FONC_DESORP'
        call rcvalb(fami, kpg, spt, poum, imate, &
                    ' ', phenom, nbpar, nompar, valpar, &
                    4, nomres, valres, icodre, 1)
        rap = ((1.d0-valres(4))/0.25d0)**valres(3)
        difl = valres(1)*(valres(2)+(1.d0-valres(2))/(1.d0+rap))
!
    else if (phenom .eq. 'SECH_RFT') then
        call rftDiffusion(fami, kpg, spt, poum, imate, &
                          c, temp, diff, difl, difv)
!
    else if (phenom .eq. 'SECH_NAPPE') then
        nbpar = 2
        nompar(1) = 'SECH'
        valpar(1) = c
!       pour SECH_NAPPE il faut donner TEMP en paramètre pour avoir
!       la valeur du bon instant selon l'option
!       sinon il faudrait voir globalement s'il ne serait pas plus
!       correct de mettre poum = '-' pour l'option CHAR_THER_EVOLNI
        nompar(2) = 'TEMP'
        valpar(2) = temp
        nomres(1) = 'FONCTION'
        call rcvalb(fami, kpg, spt, poum, imate, &
                    ' ', phenom, nbpar, nompar, valpar, &
                    1, nomres, valres, icodre, 1)
        difl = valres(1)
!
    else
        call utmess('F', 'ALGORITH10_20', sk=comp)
    end if

    diff = difl+difv

    if (present(difl_)) then
        difl_ = difl
    end if

    if (present(difv_)) then
        difv_ = difv
    end if

!
!
end subroutine
