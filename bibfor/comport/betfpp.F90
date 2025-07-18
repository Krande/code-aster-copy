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

subroutine betfpp(BEHinteg, &
                  materf, nmat, pc, pt, &
                  nseuil, fc, ft, dfcdlc, dftdlt, &
                  kuc, kut, ke)
!
    use Behaviour_type
!
    implicit none
!
!       BETON_DOUBLE_DP: CONVEXE ELASTO PLASTIQUE POUR (MATER,SIG,P1,P2)
!                   AVEC UN SEUIL EN COMPRESSION ET UN SEUIL EN TRACTION
!       CALCUL DES VALEURS DES COURBES D'ADOUCISSEMENT ET DES DERIVES
!       PAR RAPPORT AUX INCREMENTS DE MULTIPLICATEURS PLASTIQUES
!       IN  MATERF :  COEFFICIENTS MATERIAU
!           NMAT   :  DIMENSION MATERF
!           PC     :  MULTIPLICATEUR PLASTIQUE EN COMPRESSION
!           PT     :  MULTIPLICATEUR PLASTIQUE EN TRACTION
!           NSEUIL :  SEUIL D'ELASTICITE ACTIVE
!                     NSEUIL = 1  -->  SEUIL COMPRESSION ACTIF
!                     NSEUIL = 2  -->  SEUIL TRACTION ACTIF
!                     NSEUIL = 3  -->  SEUIL COMPRESSION ET TRACTION
!                                                             ACTIFS
!       OUT FC     :  ECCROUISSAGE EN COMPRESSION
!           FT     :  ECCROUISSAGE EN TRACTION
!           DFCDLC :  DERIVE DE LA COURBE D'ADOUCISSEMENT EN COMPRESSION
!           DFTDLT :  DERIVE DE LA COURBE D'ADOUCISSEMENT EN TRACTION
!           KUC    :  ECROUISSAGE ULTIME EN COMPRESSION
!           KUT    :  ECROUISSAGE ULTIME EN TRACTION
!       ----------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
    real(kind=8) :: un, zero, d13, deux
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(zero=0.d0)
    parameter(d13=.33333333333333d0)
!
    type(Behaviour_Integ), intent(in) :: BEHinteg
    integer(kind=8) :: nmat, nseuil
    real(kind=8) :: materf(nmat, 2)
    real(kind=8) :: pc, pt, dfcdlc, dftdlt, kuc, kut
!
    real(kind=8) :: lc, fcp, ftp, fc, ft
    real(kind=8) :: gc, gt, celas
    real(kind=8) :: e, ku, ke
    real(kind=8) :: lc0, epsi
    integer(kind=8) :: typcom, typtra, iadzi, iazk24
    character(len=8) :: nomail
!     ------------------------------------------------------------------
    integer(kind=8) :: n, nd
    common/tdim/n, nd
    data epsi/1.d-6/
!     ------------------------------------------------------------------
!
! --- INITIALISATION
!
    e = materf(1, 1)
    fcp = materf(1, 2)
    ftp = materf(2, 2)
    gc = materf(4, 2)
    gt = materf(5, 2)
    celas = materf(6, 2)
    typcom = int(materf(7, 2)+0.5d0)
    typtra = int(materf(8, 2)+0.5d0)
!
    kuc = zero
    ke = zero
    kut = zero
    fc = zero
    ft = zero
    dfcdlc = zero
    dftdlt = zero
!
! --- LONGUEUR CARACTERISTIQUE POUR LOI BETON LC
!
    if (materf(9, 2) .lt. zero) then
        if (BEHinteg%behavESVA%behavESVAGeom%lElemSize1) then
            lc = BEHinteg%behavESVA%behavESVAGeom%elemSize1
        else
            call utmess('F', 'COMPOR2_12')
        end if
    else
        lc = materf(9, 2)
    end if
!
!
! --- VALEUR ET DERIVEE DE LA COURBE D'ADOUCISSEMENT EN COMPRESSION
!
    if (nseuil .eq. 1 .or. nseuil .eq. 3 .or. nseuil .eq. 11 .or. nseuil .eq. 33) then
!
! -      COURBE POST PIC EN COMPRESSION LINEAIRE
!
        if (typcom .eq. 0) then
            ke = deux*(un-celas)*fcp/e
            ku = deux*gc/(lc*fcp)-(un+deux*celas)*ke*d13
            lc0 = (6.d0*e*gc)/(fcp*fcp*(11.d0-4.d0*celas*(un+celas)))
            if (lc .gt. lc0) then
                if (materf(9, 2) .lt. zero) then
                    call utmess('A', 'ALGORITH_44')
                else
                    call tecael(iadzi, iazk24)
                    nomail = zk24(iazk24-1+3) (1:8)
                    call utmess('A', 'ALGORITH_45', sk=nomail)
                end if
            end if
            if (pc .lt. ke) then
                fc = fcp*(celas+deux*(un-celas)*pc/ke+(celas-un)*pc*pc/(ke*ke))
                dfcdlc = fcp*deux*(un-celas)*(un-pc/ke)/ke
            else
                if (pc .lt. ku) then
                    fc = fcp*(pc-ku)/(ke-ku)
                    dfcdlc = fcp/(ke-ku)
                else
                    fc = fcp*epsi
!               FC = ZERO
                    dfcdlc = zero
                end if
            end if
!
! -      COURBE POST PIC EN COMPRESSION NON LINEAIRE
!
        else
            ke = deux*(un-celas)*fcp/e
            ku = 1.5d0*gc/(lc*fcp)-0.5d0*celas*ke
            lc0 = (1.5d0*e*gc)/(fcp*fcp*(4.d0-celas*(un+celas)))
            if (lc .gt. lc0) then
                if (materf(9, 2) .lt. zero) then
                    call utmess('A', 'ALGORITH_44')
                else
                    call tecael(iadzi, iazk24)
                    nomail = zk24(iazk24-1+3) (1:8)
                    call utmess('A', 'ALGORITH_45', sk=nomail)
                end if
            end if
            if (pc .lt. ke) then
                fc = fcp*(celas+deux*(un-celas)*pc/ke+(celas-un)*pc*pc/(ke*ke))
                dfcdlc = fcp*deux*(un-celas)*(un-pc/ke)/ke
            else if (pc .lt. ku) then
                fc = fcp*(un-(pc-ku)*(pc-ku)/((ke-ku)*(ke-ku)))
                dfcdlc = -deux*fcp*(pc-ku)/((ke-ku)*(ke-ku))
            else
                fc = fcp*epsi
!               FC = ZERO
                dfcdlc = zero
            end if
        end if
        kuc = ku
    end if
!
! --- VALEUR ET DERIVEE DE LA COURBE D'ADOUCISSEMENT EN TRACTION
!
    if (nseuil .eq. 2 .or. nseuil .eq. 3 .or. nseuil .eq. 22 .or. nseuil .eq. 33) then
!
! -      COURBE POST PIC EN TRACTION LINEAIRE
!
        if (typtra .eq. 0) then
            ku = deux*gt/(lc*ftp)
            lc0 = (deux*e*gt)/(ftp*ftp)
            if (lc .gt. lc0) then
                if (materf(9, 2) .lt. zero) then
                    call utmess('A', 'ALGORITH_46')
                else
                    call tecael(iadzi, iazk24)
                    nomail = zk24(iazk24-1+3) (1:8)
                    call utmess('A', 'ALGORITH_47', sk=nomail)
                end if
            end if
            if (pt .lt. ku) then
                ft = ftp*(un-pt/ku)
                dftdlt = -ftp/ku
            else
                ft = ftp*epsi
                dftdlt = zero
            end if
!
! -      COURBE POST PIC EN TRACTION NON LINEAIRE
!
        else
            ku = 1.d06
            lc0 = (e*gt)/(ftp*ftp)
            if (lc .gt. lc0) then
                if (materf(9, 2) .lt. zero) then
                    call utmess('A', 'ALGORITH_46')
                else
                    call tecael(iadzi, iazk24)
                    nomail = zk24(iazk24-1+3) (1:8)
                    call utmess('A', 'ALGORITH_47', sk=nomail)
                end if
            end if
            ft = ftp*exp(-lc*ftp*pt/gt)
            dftdlt = -ftp*ftp*lc/gt*exp(-lc*ftp*pt/gt)
        end if
        kut = ku
!
    end if
!
end subroutine
