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
subroutine lcpitg(compor, df, line, dp, dvbe, &
                  dtaudf)
!
!
    implicit none
!
! ----------------------------------------------------------------------
!        INTEGRATION DE LA LOI SIMO MIEHE ECROUISSAGE ISOTROPE
!  DERIVEE DE TAU PAR RAPPORT A DF * DFT = TAUDF(IJ,K,P)*DF(L,P)
!  TAU = TAU(DVTAU - TRTAU)
!  TRTAU = TRTAU(DF)
!  DVTAU=DVTAU(DVBE)
!  DVBE = DVBE(BETR)
!  BETR =  ETR(DFB)
!  DFB=DVB(DF)
! ----------------------------------------------------------------------
! IN  COMPOR : COMPORTEMENT
! IN  DF     : INCREMENT DU TENSEUR DE DEFORMATION
! IN  LINE : REGIME DE LA SOLUTION (ELASTIQUE, PLASTIQUE)
! IN  DP     : INCREMENT DE DEFORMATION PLASTIQUE
! IN  DVBE   : PARTIE DEVIATORIQUE DE LA DEFORMATION
! OUT DTAUDF : DERIVEE DE TAU PAR RAPPORT A DF * DFT
! ----------------------------------------------------------------------
#include "blas/ddot.h"
    character(len=16) :: compor
    integer(kind=8) :: line
    real(kind=8) :: df(3, 3), dp, dvbe(6), dtaudf(6, 3, 3)
!
! COMMON GRANDES DEFORMATIONS SIMO - MIEHE
!
    integer(kind=8) :: ind(3, 3), ind1(6), ind2(6)
    real(kind=8) :: kr(6), rac2, rc(6), id(6, 6)
    real(kind=8) :: bem(6), betr(6), dvbetr(6), eqbetr, trbetr
    real(kind=8) :: jp, dj, jm, dfb(3, 3)
    real(kind=8) :: djdf(3, 3), dbtrdf(6, 3, 3)
!
    common/gdsmc/&
     &            bem, betr, dvbetr, eqbetr, trbetr,&
     &            jp, dj, jm, dfb,&
     &            djdf, dbtrdf,&
     &            kr, id, rac2, rc, ind, ind1, ind2
! ----------------------------------------------------------------------
!  COMMON MATERIAU POUR VON MISES
!
    integer(kind=8) :: jprol, jvale, nbval
    real(kind=8) :: pm, young, nu, mu, unk, troisk, cother
    real(kind=8) :: sigm0, epsi0, dt, coefm, rpm, pente, apui, npui, sigy
!
    common/lcpim/&
     &          pm, young, nu, mu, unk, troisk, cother,&
     &          sigm0, epsi0, dt, coefm, rpm, pente,&
     &          apui, npui, sigy, jprol, jvale, nbval
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: ij, kl, k, l
    real(kind=8) :: a2, a3, a4, rb, arg, eqbe
    real(kind=8) :: dtaudj, dtaudb
    real(kind=8) :: dvbbtr(6, 6), dvbedf(6, 3, 3)
    blas_int :: b_incx, b_incy, b_n
! ----------------------------------------------------------------------
!
! 1 - DEFINITION DES COEFFICIENTS UTILES
!
!
    if (line .eq. 0) then
        a2 = 1
        a3 = 0
        a4 = 0
    else
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        eqbe = sqrt(1.5d0*ddot(b_n, dvbe, b_incx, dvbe, b_incy))
        a2 = 1/(1+trbetr*dp/eqbe)
!
        rb = pente
        if (compor(1:4) .eq. 'VISC' .and. dp .ne. 0.d0) then
            arg = (dp/(dt*epsi0))**(1.d0/coefm)
            rb = rb+sigm0*arg/(sqrt(arg**2+1.d0)*coefm*dp)
        end if
!
        a3 = (a2*dp/eqbe)-mu/(rb+mu*trbetr)
        a3 = 3.d0*trbetr*a3/(2.d0*eqbe*eqbe)
!
        a4 = dp/eqbe*((mu*trbetr/(rb+mu*trbetr))-1.d0)
    end if
!
!
! 2 - DERIVEE DE DVBE PAR RAPPORT A BETR = DVBBTR
!
    do ij = 1, 6
        do kl = 1, 6
            dvbbtr(ij, kl) = a2*( &
                             id(ij, kl)-kr(ij)*kr(kl)/3.d0)+a3*dvbe(ij)*dvbe(kl)+a4*dvbe(ij)*kr(&
                             &kl &
                             )
        end do
    end do
!
!
! 3 - DERIVEE DE DVBE PAR RAPPORT A DF = DVBEDF
!
    do ij = 1, 6
        do k = 1, 3
            do l = 1, 3
                b_n = to_blas_int(6)
                b_incx = to_blas_int(6)
                b_incy = to_blas_int(1)
                dvbedf(ij, k, l) = ddot(b_n, dvbbtr(ij, 1), b_incx, dbtrdf(1, k, l), b_incy)
            end do
        end do
    end do
!
!
! 4 - MATRICE TANGENTE = DTAUDF
!
!    DERIVEE PARTIELLE DE TAU PAR RAPPORT A B ET J
    dtaudb = mu
    dtaudj = 0.5d0*(2.d0*unk*jp-cother*(1.d0-1.d0/(jp**2.d0)))
!
    do ij = 1, 6
        do k = 1, 3
            do l = 1, 3
                dtaudf(ij, k, l) = dtaudb*dvbedf(ij, k, l)+dtaudj*kr(ij)*djdf(k, l)
            end do
        end do
    end do
!
end subroutine
