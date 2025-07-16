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
subroutine burftm(cmp, ndim, vim, epsfm)
! person_in_charge: alexandre.foucault at edf.fr
! ----------------------------------------------------------------------
! ROUTINE DE COMPORTEMENT DE FLUAGE BETON_BURGER
! PASSAGE DES VARIABLES INTERNES AU VECTEUR DE DEFORMATIONS DE FLUAGE
!
! IN  CMP : COMPOSANTES A CALCULER :
!             'FP' : FLUAGE PROPRE UNIQUEMENT
!             'FD' : FLUAGE DE DESSICCATION
!     NDIM  : DIMENSION DU PROBLEME (2D OU 3D)
!     VIM   : VARIABLES INTERNES
! OUT EPSFM : VECTEUR DE DEFORMATIONS DE FLUAGE
!
!  STRUCTURE DES VARIABLES INTERNES MFRONT : VIM,VIP ( X = I ou F )
!
!     DANS LE CAS 2D :
!     VIX(1)     = ElasticStrainXX : DEFORMATION ELASTIQUE 11
!     VIX(2)     = ElasticStrainYY : DEFORMATION ELASTIQUE 22
!     VIX(3)     = ElasticStrainZZ : DEFORMATION ELASTIQUE 33
!     VIX(4)     = ElasticStrainXY : DEFORMATION ELASTIQUE 12
!     VIX(5)     = ESPHI   : DEFORMATION DE FLUAGE IRR SPHERIQUE
!     VIX(6)     = ELIM    : DEF. EQUIVALENTE DU FP IRREVERSIBLE MAX
!     VIX(7)    = EDEVIXX : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 11
!     VIX(8)    = EDEVIYY : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 22
!     VIX(9)    = EDEVIZZ : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 33
!     VIX(10)    = EDEVIXY : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 12
!     VIX(11)     = ESPHR   : DEFORMATION DE FLUAGE REV SPHERIQUE
!     VIX(12)    = EDEVRXX : DEFORMATION DE FLUAGE REV DEVIATORIQUE 11
!     VIX(13)    = EDEVRYY : DEFORMATION DE FLUAGE REV DEVIATORIQUE 22
!     VIX(14)    = EDEVRZZ : DEFORMATION DE FLUAGE REV DEVIATORIQUE 33
!     VIX(15)    = EDEVRXY : DEFORMATION DE FLUAGE REV DEVIATORIQUE 12
!     VIX(16)    = EdessXX : DEFORMATION DE FLUAGE DE DESSICCATION  11
!     VIX(17)    = EdessYY : DEFORMATION DE FLUAGE DE DESSICCATION  22
!     VIX(18)    = EdessZZ : DEFORMATION DE FLUAGE DE DESSICCATION  33
!     VIX(19)    = EdessXY : DEFORMATION DE FLUAGE DE DESSICCATION  12
!     VIX(20)    = EFXX    : DEFORMATION TOTALE DE FLUAGE PROPRE    11
!     VIX(21)    = EFYY    : DEFORMATION TOTALE DE FLUAGE PROPRE    22
!     VIX(22)    = EFZZ    : DEFORMATION TOTALE DE FLUAGE PROPRE    33
!     VIX(23)    = EFXY    : DEFORMATION TOTALE DE FLUAGE PROPRE    12
!     VIX(24)    = ShiftedHistoricalMinimumRelativeHumidity (rHmin)
!
!     DANS LE CAS 3D :
!     VIX(1)     = ElasticStrainXX : DEFORMATION ELASTIQUE 11
!     VIX(2)     = ElasticStrainYY : DEFORMATION ELASTIQUE 22
!     VIX(3)     = ElasticStrainZZ : DEFORMATION ELASTIQUE 33
!     VIX(4)     = ElasticStrainXY : DEFORMATION ELASTIQUE 12
!     VIX(5)     = ElasticStrainXZ : DEFORMATION ELASTIQUE 13
!     VIX(6)     = ElasticStrainYZ : DEFORMATION ELASTIQUE 23
!     VIX(7)     = ESPHI   : DEFORMATION DE FLUAGE IRR SPHERIQUE
!     VIX(8)     = ELIM    : DEF. EQUIVALENTE DU FP IRREVERSIBLE MAX
!     VIX(9)    = EDEVIXX : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 11
!     VIX(10)    = EDEVIYY : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 22
!     VIX(11)    = EDEVIZZ : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 33
!     VIX(12)    = EDEVIXY : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 12
!     VIX(13)    = EDEVIXZ : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 13
!     VIX(14)    = EDEVIYZ : DEFORMATION DE FLUAGE IRE DEVIATORIQUE 23
!     VIX(15)     = ESPHR   : DEFORMATION DE FLUAGE REV SPHERIQUE
!     VIX(16)    = EDEVRXX : DEFORMATION DE FLUAGE REV DEVIATORIQUE 11
!     VIX(17)    = EDEVRYY : DEFORMATION DE FLUAGE REV DEVIATORIQUE 22
!     VIX(18)    = EDEVRZZ : DEFORMATION DE FLUAGE REV DEVIATORIQUE 33
!     VIX(19)    = EDEVRXY : DEFORMATION DE FLUAGE REV DEVIATORIQUE 12
!     VIX(20)    = EDEVRXZ : DEFORMATION DE FLUAGE REV DEVIATORIQUE 13
!     VIX(21)    = EDEVRYZ : DEFORMATION DE FLUAGE REV DEVIATORIQUE 23
!     VIX(22)    = EdessXX : DEFORMATION DE FLUAGE DE DESSICCATION  11
!     VIX(23)    = EdessYY : DEFORMATION DE FLUAGE DE DESSICCATION  22
!     VIX(24)    = EdessZZ : DEFORMATION DE FLUAGE DE DESSICCATION  33
!     VIX(25)    = EdessXY : DEFORMATION DE FLUAGE DE DESSICCATION  12
!     VIX(26)    = EdessXZ : DEFORMATION DE FLUAGE DE DESSICCATION  13
!     VIX(27)    = EdessYZ : DEFORMATION DE FLUAGE DE DESSICCATION  23
!     VIX(28)    = EFXX    : DEFORMATION TOTALE DE FLUAGE PROPRE    11
!     VIX(29)    = EFYY    : DEFORMATION TOTALE DE FLUAGE PROPRE    22
!     VIX(30)    = EFZZ    : DEFORMATION TOTALE DE FLUAGE PROPRE    33
!     VIX(31)    = EFXY    : DEFORMATION TOTALE DE FLUAGE PROPRE    12
!     VIX(32)    = EFXZ    : DEFORMATION TOTALE DE FLUAGE PROPRE    13
!     VIX(33)    = EFYZ    : DEFORMATION TOTALE DE FLUAGE PROPRE    23
!     VIX(34)    = ShiftedHistoricalMinimumRelativeHumidity (rHmin)
! ----------------------------------------------------------------------
    implicit none
    character(len=*) :: cmp
    real(kind=8) :: vim(33), epsfm(6)
    integer(kind=8) :: ndim, i
!
    if (cmp(1:2) .eq. 'FP') then
        if (ndim .eq. 3) then
            do i = 1, 3
                epsfm(i) = (vim(15)+vim(7))
            end do
            epsfm(1) = epsfm(1)+vim(16)+vim(9)
            epsfm(2) = epsfm(2)+vim(17)+vim(10)
            epsfm(3) = epsfm(3)+vim(18)+vim(11)
            epsfm(4) = vim(19)+vim(12)
            epsfm(5) = vim(20)+vim(13)
            epsfm(6) = vim(21)+vim(14)
        else if (ndim .eq. 2) then
            do i = 1, 3
                epsfm(i) = (vim(11)+vim(5))
            end do
            epsfm(1) = epsfm(1)+vim(12)+vim(7)
            epsfm(2) = epsfm(2)+vim(13)+vim(8)
            epsfm(3) = epsfm(3)+vim(14)+vim(9)
            epsfm(4) = vim(15)+vim(10)
            epsfm(5) = 0.
            epsfm(6) = 0.
        end if
    else if (cmp(1:2) .eq. 'FD') then
        if (ndim .eq. 3) then
            epsfm(1) = vim(22)
            epsfm(2) = vim(23)
            epsfm(3) = vim(24)
            epsfm(4) = vim(25)
            epsfm(5) = vim(27)
            epsfm(6) = vim(26)
        else if (ndim .eq. 2) then
            epsfm(1) = vim(16)
            epsfm(2) = vim(17)
            epsfm(3) = vim(18)
            epsfm(4) = vim(19)
            epsfm(5) = 0.
            epsfm(6) = 0.
        end if
    end if
!
end subroutine
