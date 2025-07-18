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

subroutine lcmafl(fami, kpg, ksp, poum, nmater, &
                  imat, necoul, nbval, valres, nmat, &
                  itbint, nfs, nsg, hsri, nbsys)
    implicit none
!     MONOCRISTAL : RECUPERATION DU MATERIAU A T(TEMPD) ET T+DT(TEMPF)
!                  MATER(*,2) = COEF ECOULEMENT VISCOPLASTIQUE
!     ----------------------------------------------------------------
!     IN  IMAT   :  ADRESSE DU MATERIAU CODE
!         NMATER :  NOM DU MATERIAU
!         NMAT   :  DIMENSION  DE MATER
!         NECOUL :  NOM DE LA LOI D'ECOULEMENT
!         VALPAR :  VALEUR DES PARAMETRES
!         NOMPAR :  NOM DES PARAMETRES
!     OUT VALRES :  COEFFICIENTS MATERIAU A T
!         NBVAL  :  NOMBRE DE COEF MATERIAU LUS
!     ----------------------------------------------------------------
#include "asterfort/lcmhsr.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: kpg, ksp, nmat, imat, nbval, nbcoef, itbint, nfs, nsg
    integer(kind=8) :: iret2, nbsys
    real(kind=8) :: valres(nmat), hsri(nsg, nsg), h, e, nu, mu
    real(kind=8) :: tempf, valh(6), vallue(nmat), val(1)
    character(len=16) :: nomres(nmat)
    integer(kind=8) :: icodre(nmat)
    character(len=*) :: fami, poum
    character(len=16) :: nmater, necoul, necris, phenom
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
!
    if (necoul .eq. 'MONO_VISC1') then
        nbval = 3
        nomres(1) = 'N'
        nomres(2) = 'K'
        nomres(3) = 'C'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECOU_VISC1 A LE NUMERO 1
        valres(1) = 1
!
    end if
    if (necoul .eq. 'MONO_VISC2') then
        nbval = 5
        nomres(1) = 'N'
        nomres(2) = 'K'
        nomres(3) = 'C'
        nomres(4) = 'A'
        nomres(5) = 'D'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECOU_VISC2 A LE NUMERO 2
        valres(1) = 2
!
    end if
    if (necoul .eq. 'MONO_DD_CFC') then
        nbval = 6
        nomres(1) = 'TAU_F'
        nomres(2) = 'GAMMA0'
        nomres(3) = 'A'
        nomres(4) = 'B'
        nomres(5) = 'N'
        nomres(6) = 'Y'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECOU_DD_CFC A LE NUMERO 5
        valres(1) = 5
        nbval = nbval+1
        valres(nbval) = 0.d0
    end if
    if (necoul .eq. 'MONO_DD_CFC_IRRA') then
        nbval = 6
        nomres(1) = 'TAU_F'
        nomres(2) = 'GAMMA0'
        nomres(3) = 'A'
        nomres(4) = 'B'
        nomres(5) = 'N'
        nomres(6) = 'Y'
!
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECOU_DD_CFC_IRRA A LE NUMERO 8
        valres(1) = 8
        nbval = nbval+1
        valres(nbval) = 0.d0
    end if
    if (necoul .eq. 'MONO_DD_FAT') then
        nbval = 7
        nomres(1) = 'TAU_F'
        nomres(2) = 'GAMMA0'
        nomres(3) = 'BETA'
        nomres(4) = 'UN_SUR_D'
        nomres(5) = 'N'
        nomres(6) = 'GC0'
        nomres(7) = 'K'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECOU_ECP_CFC A LE NUMERO 6
        valres(1) = 6
!
        nbval = nbval+1
        valres(nbval) = 0.d0
!
    end if
    if (necoul(1:10) .eq. 'MONO_DD_CC') then
        nbval = 15
        nomres(1) = 'B'
        nomres(2) = 'GH'
        nomres(3) = 'DELTAG0'
        nomres(4) = 'TAU_0'
        nomres(5) = 'D'
        nomres(6) = 'GAMMA0'
        nomres(7) = 'N'
        nomres(8) = 'DEPDT'
        nomres(9) = 'Y_AT'
        nomres(10) = 'D_LAT'
        nomres(11) = 'K_F'
        nomres(12) = 'K_SELF'
        nomres(13) = 'TAU_F'
        nomres(14) = 'RHO_MOB'
        nomres(15) = 'K_BOLTZ'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
!
!       LECTURE OU CALCUL ET STOCKAGE DE MU
        nomres(16) = 'MU_MOY'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    1, nomres(16), vallue(16), icodre(16), 0)
        if (icodre(16) .eq. 0) then
            mu = vallue(16)
        else
            call rccoma(imat, 'ELAS', 1, phenom, icodre(1))
            if (phenom .eq. 'ELAS') then
                call rcvalb(fami, kpg, ksp, poum, imat, &
                            ' ', 'ELAS', 0, ' ', [0.d0], &
                            1, 'E', val, icodre, 1)
                e = val(1)
                call rcvalb(fami, kpg, ksp, poum, imat, &
                            ' ', 'ELAS', 0, ' ', [0.d0], &
                            1, 'NU', val, icodre, 1)
                nu = val(1)
                mu = e/(2.0d0+2.0d0*nu)
            else
                call utmess('F', 'COMPOR1_88', sk=necoul)
            end if
        end if
        call rcvarc('F', 'TEMP', poum, fami, kpg, &
                    ksp, tempf, iret2)
        if (iret2 .ne. 0) then
            call utmess('F', 'COMPOR1_82', sk=necoul)
        end if
!
        nbval = nbval+1
        vallue(nbval) = tempf
        nbval = nbval+1
        vallue(nbval) = mu
        nbval = nbval+1
        vallue(nbval) = 0.d0
        irr = 0
        if (necoul .eq. 'MONO_DD_CC_IRRA') then
            irr = 1
            vallue(nbval) = 1.d0
            nomres(1) = 'A_IRRA'
            nomres(2) = 'XI_IRRA'
            call rcvalb(fami, kpg, ksp, poum, imat, &
                        nmater, necoul, 0, ' ', [0.d0], &
                        2, nomres, vallue(nbval+1), icodre, 1)
            nbval = nbval+2
        end if
        valres(2:2-1+nbval) = vallue(1:nbval)
!         PAR CONVENTION ECOU_DD_CC A LE NUMERO 7
        nbval = nbval+1
        valres(1) = 7
        nbval = nbval+1
        valres(nbval) = 0.d0
    end if
!
    if (necoul .eq. 'MONO_DD_KR') then
        nbval = 10
        nomres(1) = 'K'
        nomres(2) = 'TAUR'
        nomres(3) = 'TAU0'
        nomres(4) = 'GAMMA0'
        nomres(5) = 'DELTAG0'
        nomres(6) = 'BSD'
        nomres(7) = 'GCB'
        nomres(8) = 'KDCS'
        nomres(9) = 'P'
        nomres(10) = 'Q'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION KOCKS_RAUCH A LE NUMERO 4
        valres(1) = 4
!
        call rcvarc('F', 'TEMP', poum, fami, kpg, &
                    ksp, tempf, iret2)
        if (iret2 .ne. 0) then
            call utmess('F', 'COMPOR1_82', sk=necoul)
        end if
        nbval = nbval+1
        valres(nbval) = tempf
!
!
!         DEFINITION DE LA MATRICE D'INTERACTION POUR KOCKS-RAUCH
        nomres(1) = 'H'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necoul, 0, ' ', [0.d0], &
                    1, nomres, val, icodre, 0)
        h = val(1)
        if (icodre(1) .eq. 0) then
            nbcoef = 1
            valh(1) = h
        else
            nomres(1) = 'H1'
            nomres(2) = 'H2'
            nomres(3) = 'H3'
            nomres(4) = 'H4'
            nomres(5) = 'H5'
            nomres(6) = 'H6'
            call rcvalb(fami, kpg, ksp, poum, imat, &
                        nmater, necoul, 0, ' ', [0.d0], &
                        6, nomres, valh, icodre, 0)
            if (icodre(5) .eq. 0) then
                nbcoef = 6
            else
                nbcoef = 4
            end if
!
        end if
        if (itbint .eq. 0) then
            necris = necoul
            call lcmhsr(necoul, necris, nbsys, nbcoef, valh, &
                        nsg, hsri)
        end if
!
    end if
end subroutine
