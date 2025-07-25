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

subroutine dgseui(em, num, ef, nuf, eb, &
                  nub, ftj, h, icisai, syt, &
                  dxd, syf, drd, pelast, &
                  pelasf)
!
    implicit none
!
! PARAMETRES ENTRANTS
    integer(kind=8) :: icisai
    real(kind=8) :: em, num, ef, nuf, h
    real(kind=8) :: eb, nub, ftj
!
! PARAMETRES SORTANTS
    real(kind=8) :: syt, dxd, syf, drd, pelast, pelasf
!
! person_in_charge: sebastien.fayolle at edf.fr
! ----------------------------------------------------------------------
!
! BUT : DETERMINATION DES SEUIL D ENDOMMAGEMENT
!
! IN:
!       ICOMPR : METHODE DE COMPRESSION
!       EM     : MODULE D YOUNG EN MEMBRANE
!       EF     : MODULE D YOUNG EN FLEXION
!       NUM    : COEFF DE POISSON EN MEMBRANE
!       NUF    : COEFF DE POISSON EN FLEXION
!       H      : EPAISSEUR DE LA PLAQUE
!       EB     : MODULE D YOUNG DU BETON
!       NUB    : COEFF DE POISSON DU BETON
!       FTJ    : LIMITE A LA TRACTION DU BETON
! OUT:
!       SYT    : SEUIL D'ENDOMMAGEMENT EN TRACTION
!       DXD    : DEPLACEMENT A L'APPARITION DE L'ENDOMMAGEMENT
!       SYF    : SEUIL D'ENDOMMAGEMENT EN FLEXION
!       DRD    : ROTATION A L'APPARITION DE L'ENDOMMAGEMENT
!       PELAST : PENTE ELASTIQUE EN TRACTION
!       PELASF : PENTE ELASTIQUE EN FLEXION
! ----------------------------------------------------------------------
!
! - DETERMINATION DES PARAMETRES EN MEMBRANE
! - SEUILS D'ENDOMMAGEMENT EN TRACTION PURE
    syt = ftj*em*h*(1.d0-nub**2)/(eb*(1.d0-nub*num))
!
! - DEPLACEMENT A L'APPARITION DE L'ENDOMMAGEMENT
    dxd = syt/(h*em)
!
    if (icisai .eq. 1) then
! - CALCUL DE LA PENTE ELASTIQUE EN CISAILLEMENT PUR DANS LE PLAN
        pelast = em/(1.d0+num)/2.d0*h
    else
! - PENTE ELASTIQUE EN TRACTION
        pelast = em*h
    end if
!
! - DETERMINATION DES PARAMETRES EN FLEXION
! - SEUIL D'ENDOMMAGEMENT EN FLEXION PURE
    syf = ftj*ef*h**2*(1.d0-nub**2)/(6.d0*eb*(1.d0-nub*nuf))
! - ROTATION A L'APPARITION DE L'ENDOMMAGEMENT
    drd = 12.d0*syf/(h**3*ef)
    pelasf = ef*h**3/12.d0
!
end subroutine
