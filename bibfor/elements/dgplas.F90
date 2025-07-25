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

subroutine dgplas(ea, sya, eb, nub, ftj, fcj, &
                  num, nuf, a, b1, b, &
                  syt, syf, ef, dxd, drd, h, &
                  ipentetrac, ipenteflex, icisai, emaxm, emaxf, nnap, &
                  omx, rx, ry, np, dxp, pendt, &
                  drp, mp, pendf)
!
! aslint: disable=W1504
    implicit none
!
! PARAMETRES ENTRANTS
#include "asterfort/dgmmax.h"
#include "asterfort/dgmpla.h"
#include "asterfort/calc_myf_gf.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nnap, ilit, icisai, ipentetrac, ipenteflex
!
    real(kind=8) :: ea(*), sya(*), eb, nub, num, nuf, w, emaxm, emaxf, ef
    real(kind=8) :: a, b1, b, syt, syf, dxd, drd, h, c, omx, rx(*), ry(*)
    real(kind=8) :: rmesg(2), ftj, fcj
!
! PARAMETRES SORTANTS
    real(kind=8) :: pendt, pendf
!
!
! person_in_charge: sebastien.fayolle at edf.fr
! ----------------------------------------------------------------------
!
! BUT : DETERMINATION DES PENTES POST ENDOMMAGEMENT
!
! IN:
!       EA       : MODULES D YOUNG DES ACIERS
!       SYA      : LIMITES ELASTIQUES DES ACIERS
!       EB       : MODULE D YOUNG DU BETON
!       NUB      : COEFF DE POISSON DU BETON
!       FTJ      : LIMITE A LA TRACTION DU BETON
!       NUM      : COEFF DE POISSON EN MEMBRANE
!       NUF      : COEFF DE POISSON EN FLEXION
!       A
!       B1
!       B        : SECTIONS DES ACIERS
!       SYT      : SEUIL D'ENDOMMAGEMENT EN TRACTION
!       DXD      : DEPLACEMENT A L'APPARITION DE L'ENDOMMAGEMENT
!       DRD      : ROTATION A L'APPARITION DE L'ENDOMMAGEMENT
!       H        : EPAISSEUR DE LA PLAQUE
!       IPENTETRAC   : OPTION DE CALCUL DES PENTES POST ENDOMMAGEMENT TR
!                  1 : RIGI_ACIER
!                  2 : PLAS_ACIER
!                  3 : UTIL
!       IPENTEFLEX   : OPTION DE CALCUL DES PENTES POST ENDOMMAGEMENT FL
!                  1 : RIGI_INIT
!                  2 : UTIL
!                  3 : RIGI_ACIER
!       ICISAI   : INDICATEUR DE CISAILLEMENT
!       EMAXM    : DEFO GENE MAX EN MEMBRANE
!       EMAXF    : DEFO GENE MAX EN FLEXION
!       NNAP     : NOMBRE DE NAPPE
!       RX       : POSITION ADIMENSIONNEE DU LIT DE CABLES SUIVANT X
!       RY       : POSITION ADIMENSIONNEE DU LIT DE CABLES SUIVANT Y
!       NP       : EFFORT A PLASTICITE
!       DXP      : DEPLACEMENT A PLASTICITE
!       DRP      : ROTATION A PLASTICITE
!       MP       : MOMENT A PLASTICITE
!
! OUT:
!       PENDT    : PENTE POST ENDOMMAGEMENT EN MEMBRANNE
!       PENDF    : PENTE POST ENDOMMAGEMENT EN FLEXION

! IN/OUT :
!       SYF      : SEUIL D'ENDOMMAGEMENT EN FLEXION
! ----------------------------------------------------------------------
!
    real(kind=8) :: np, dxp, mp, drp, ya, efm
!
! - DETERMINATION DE LA PENTE POST ENDOMMAGEMENT EN TRACTION
    if (ipentetrac .eq. 3) then
        if (emaxm .lt. dxd) then
            rmesg(1) = emaxm
            rmesg(2) = dxd
            call utmess('F', 'ALGORITH6_5', nr=2, valr=rmesg)
        end if
        dxp = emaxm
        np = b*dxp
        pendt = (np-syt)/(dxp-dxd)
    else if (ipentetrac .eq. 1) then
        pendt = b
    else if (ipentetrac .eq. 2) then
        dxp = sya(1)/ea(1)
        do ilit = 2, nnap
            if (sya(ilit)/ea(ilit) .lt. dxp) then
                dxp = sya(ilit)/ea(ilit)
            end if
        end do
        np = b*dxp
        pendt = (np-syt)/(dxp-dxd)
    end if
!
! - ESSAI DE CISAILLEMENT PUR DANS LE PLAN
    if (icisai .eq. 1) then
        if (ipentetrac .eq. 1) then
            pendt = b
        else if (ipentetrac .eq. 3) then
            if (emaxm .lt. dxd) then
                rmesg(1) = emaxm
                rmesg(2) = dxd
                call utmess('F', 'ALGORITH6_5', nr=2, valr=rmesg)
            end if
            dxp = emaxm
            np = b*dxd+ftj/3.d0
        else if (ipentetrac .eq. 2) then
            dxp = sqrt(2.d0)*dxp+2.d0*dxd
            np = b*dxp+ftj/3.d0
        end if
    end if
!
! - DETERMINATION DE LA PENTE POST ENDOMMAGEMENT EN FLEXION

    if (ipenteflex .eq. 3) then
        drp = 0.d0
        call dgmmax(eb, nub, num, nuf, h, &
                    a, b1, b, mp, drp, &
                    w, c)
        pendf = c
    else if (ipenteflex .eq. 4) then
! PLAS_ACIER
        call dgmpla(eb, nub, ea, sya, num, &
                    nuf, h, a, b1, b, &
                    nnap, rx, ry, mp, drp, &
                    w)
        pendf = (mp-syf)/(drp-drd)
    else
        ya = rx(1)*h
        efm = 1./12.*ef*h**3-2*ea(1)*omx*(ya**2)
        efm = efm*12/(h**3)
        call calc_myf_gf(efm, ftj, fcj, h, ea(1), omx, &
                         ya, sya(1), ipenteflex, emaxf, &
                         syf, pendf)

    end if
!
end subroutine
