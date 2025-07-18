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
subroutine lcmaei(fami, kpg, ksp, poum, nmater, &
                  imat, necris, necoul, nbval, valres, &
                  nmat, itbint, nfs, nsg, hsri, &
                  ifa, nomfam, nbsys)
    implicit none
!     MONOCRISTAL : RECUPERATION DU MATERIAU A T(TEMPD) ET T+DT(TEMPF)
!                  MATER(*,2) = COEF ECRO ISOT ET CALCUL DE LA
!                  MATRICE D'INTERACTION HSR
!     ----------------------------------------------------------------
!     IN  IMAT   :  ADRESSE DU MATERIAU CODE
!         NMATER :  NOM DU MATERIAU
!         NMAT   :  DIMENSION  DE MATER
!         NECRIS :  NOM DE LA LOI D'ECOULEMENT
!         IFA    :  NUMERO DE LA FAMILLE DE GLISSEMENT
!         NBCOMM :  NOMBRE DE COEF MATERIAU PAR FAMILLE
!     OUT VALRES :  COEFFICIENTS MATERIAU
!     OUT NBVAL  :  NB DE COEFFICIENTS MATERIAU
!     OUT HSR    :  MATRICE D'INTERACTION
!     ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcmhsr.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
    integer(kind=8) :: kpg, ksp, itbint, nfs, nsg
    integer(kind=8) :: nmat, nbval, imat, i, nbsys, ifa, nbcoef
    real(kind=8) :: valh(6)
    real(kind=8) :: valres(nmat), hsri(nsg, nsg), h, e, nu, mu
    real(kind=8) :: vallue(nmat), val(1)
    character(len=*) :: fami, poum
    character(len=16) :: nomres(nmat)
    integer(kind=8) :: icodre(nmat)
    character(len=16) :: nmater, necris, nomfam, necoul, phenom
    aster_logical :: zecris
!     ----------------------------------------------------------------
!
    nbval = 0
!
    if (necoul .eq. 'MONO_DD_KR') goto 999
!
    if (necris .eq. 'MONO_ISOT1') then
        nbval = 3
        nomres(1) = 'R_0'
        nomres(2) = 'Q'
        nomres(3) = 'B'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necris, 0, ' ', [0.d0], &
                    3, nomres, vallue, icodre, 1)
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECRO_ISOT1 A LE NUMERO 1
        valres(1) = 1
!
    else if (necris .eq. 'MONO_ISOT2') then
        nbval = 5
        nomres(1) = 'R_0'
        nomres(2) = 'Q1'
        nomres(3) = 'B1'
        nomres(4) = 'Q2'
        nomres(5) = 'B2'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necris, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
!
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECRO_ISOT2 A LE NUMERO 2
        valres(1) = 2
!
    else if (necris .eq. 'MONO_DD_CFC') then
        nbval = 3
        nomres(1) = 'ALPHA'
        nomres(2) = 'BETA'
        nomres(3) = 'RHO_REF'
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necris, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
!         CALCUL ET STOCKAGE DE MU
        call rccoma(imat, 'ELAS', 1, phenom, icodre(1))
!
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
            call rcvalb(fami, kpg, ksp, poum, imat, &
                        ' ', phenom, 0, ' ', [0.d0], &
                        1, 'G_LN', val, icodre, 1)
            mu = val(1)
        end if
        vallue(4) = mu
        nbval = 4
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECRO_DD_CFC A LE NUMERO 3
        valres(1) = 3
!
    else if (necris .eq. 'MONO_DD_CFC_IRRA') then
        nbval = 3+8
        nomres(1) = 'ALPHA'
        nomres(2) = 'BETA'
        nomres(3) = 'RHO_REF'
        nomres(4) = 'RHO_VOID'
        nomres(5) = 'PHI_LOOP'
        nomres(6) = 'ALP_VOID'
        nomres(7) = 'ALP_LOOP'
        nomres(8) = 'RHO_SAT'
        nomres(9) = 'PHI_SAT'
        nomres(10) = 'XI_IRRA'
        nomres(11) = 'DZ_IRRA'
!
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necris, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
!
!         CALCUL ET STOCKAGE DE MU
        call rccoma(imat, 'ELAS', 1, phenom, icodre(1))
!
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
            call rcvalb(fami, kpg, ksp, poum, imat, &
                        ' ', phenom, 0, ' ', [0.d0], &
                        1, 'G_LN', val, icodre, 1)
            mu = val(1)
        end if
        nbval = nbval+1
        vallue(nbval) = mu
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECRO_DD_CFC A LE NUMERO 3
        valres(1) = 8
!
    else if (necris(1:10) .eq. 'MONO_DD_CC') then
        nbval = 1
        nomres(1) = 'TAU_F'
!         on limite au strict minimum. tout est dans l'écoulement
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necris, 0, ' ', [0.d0], &
                    nbval, nomres, vallue, icodre, 1)
        nbval = 1
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION ECRO_DD_CC A LE NUMERO 7
        valres(1) = 7
!
!
    else if (necris(1:11) .eq. 'MONO_DD_FAT') then
!         CALCUL ET STOCKAGE DE MU
        call rccoma(imat, 'ELAS', 1, phenom, icodre(1))
!
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
            call rcvalb(fami, kpg, ksp, poum, imat, &
                        ' ', phenom, 0, ' ', [0.d0], &
                        1, 'G_LN', val, icodre, 1)
            mu = val(1)
        end if
        vallue(1) = mu
        nbval = 1
        valres(2:2-1+nbval) = vallue(1:nbval)
        nbval = nbval+1
!         PAR CONVENTION MONO_DD_FAT A LE NUMERO 4
        valres(1) = 4
!
    end if
!
    if (itbint .eq. 0) then
!
!        DEFINITION DE LA MATRICE D'INTERACTION
!        SOIT UN SEUL COEF H, SOIT H1,...H4,H5,H6
        nomres(1) = 'H'
!
        call rcvalb(fami, kpg, ksp, poum, imat, &
                    nmater, necris, 0, ' ', [0.d0], &
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
                        nmater, necris, 0, ' ', [0.d0], &
                        6, nomres, valh, icodre, 0)
!           IL FAUT AU MOINS H1 A H4
            do i = 1, 4
                if (icodre(i) .ne. 0) then
                    ASSERT(.false.)
                end if
            end do
!
            if (icodre(5) .eq. 0) then
                if (icodre(6) .eq. 0) then
                    nbcoef = 6
                else
                    nbcoef = 5
                end if
            else
                nbcoef = 4
            end if
        end if
!
!        CETTE MATRICE EST LA MEME POUR DD_CFC ET ECP_CFC
!        AUSSI POUR NE PAS MODIFIER LA ROUTINE LCMHSR (ET EN DESSOUS)
!        DANS LE CAS ECP_CFC, ON PASSE EN ARGUMENT NECRIS=MONO_DD_CFC
!        ET ON RECTIFIE APRES
!
!        LA VARIABLE LOGIQUE ZECRIS MEMORISE LA SUBSTITUTION
!
        zecris = .false.
        if (necris .eq. 'MONO_DD_FAT') then
            zecris = .true.
            necris = 'MONO_DD_CFC'
        end if
!
        call lcmhsr(necoul, necris, nbsys, nbcoef, valh, &
                    nsg, hsri)
!
        if (zecris) then
            necris = 'MONO_DD_FAT'
        end if
!
    end if
999 continue
end subroutine
