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
subroutine lcrkin(ndim, opt, rela_comp, materf, nbcomm, &
                  cpmono, nmat, mod, nvi, sigd, &
                  sigf, vind, vinf, nbphas, iret)
    implicit none
!     INITIALISATIONS POUR RUNGE-KUTTA
!     ----------------------------------------------------------------
!     IN
!          NDIM   :  2 OU 3
!          OPT    :  OPTION DE CALCUL : RIGI_MECA, FULL_MECA, RAPH_MECA
!          COMP   :  NOM MODELE DE COMPORTEMENT
!          MATERF :  COEF MATERIAU
!          NBCOMM :  INDICES DES COEF MATERIAU
!          NMAT   :  DIMENSION MATER
!          MOD    :  TYPE DE MODELISATION
!          NVI    :  NOMBRE DE VARIABLES INTERNES
!          VIND   :  VARIABLES INTERNES A T
!          SIGD   :  CONTRAINTES A T
!     VAR  NVI    :  NOMBRE DE VARIABLES INTERNES
!          SIGF   :  CONTRAINTES A T+DT
!          VINF   :  VARIABLES INTERNES A T+DT
!     OUT  IRET   :  CODE RETOUR
!     ----------------------------------------------------------------
#include "asterfort/assert.h"
#include "asterfort/lcopli.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
    integer(kind=8) :: ndt, nvi, nmat, ndi, nbcomm(nmat, 3), icp, ndim, iret, ifl
    integer(kind=8) :: indfa, nuecou
    integer(kind=8) :: nbphas, ifa, indpha, iphas, nbsys, nsfv, nbfsys
    real(kind=8) :: materf(nmat, 2), vind(*), vinf(*), id(3, 3), sigd(*)
    real(kind=8) :: sigf(*)
    real(kind=8) :: dsde(6, 6), maxdom, endoc, fp(3, 3)
    character(len=16) :: rela_comp, opt, necoul
    character(len=24) :: cpmono(5*nmat+1)
    character(len=8) :: mod
    common/tdim/ndt, ndi
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    blas_int :: b_incx, b_incy, b_n
    common/polycr/irr, decirr, nbsyst, decal, gdef
    parameter(maxdom=0.99d0)
    data id/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!     ----------------------------------------------------------------
!
    iret = 0
    dsde(:, :) = 0.d0
!
    if (materf(nmat, 1) .eq. 0) then
        call lcopli('ISOTROPE', mod, materf(1, 1), dsde)
    else if (materf(nmat, 1) .eq. 1) then
        call lcopli('ORTHOTRO', mod, materf(1, 1), dsde)
    end if
!
! --    DEBUT TRAITEMENT DE VENDOCHAB --
!     ROUTINE DE DECROISSANCE DES CONTRAINTES QUAND D>MAXDOM
    if (rela_comp .eq. 'VENDOCHAB') then
!
        if (opt .eq. 'RIGI_MECA_TANG') then
            ASSERT(ASTER_FALSE)
        end if
        if (vind(9) .ge. maxdom) then
!
            if (vind(9) .eq. 1.0d0) then
                do icp = 1, 2*ndim
                    sigf(icp) = sigd(icp)*(0.01d0)
                end do
                materf(1, 1) = 0.01d0*materf(1, 1)
                call lcopli('ISOTROPE', mod, materf(1, 1), dsde)
            else
                do icp = 1, 2*ndim
                    sigf(icp) = sigd(icp)*(0.1d0)
                end do
                endoc = (1.0d0-max(maxdom, vind(9)))*0.1d0
                materf(1, 1) = endoc*materf(1, 1)
                call lcopli('ISOTROPE', mod, materf(1, 1), dsde)
                materf(1, 1) = materf(1, 1)/endoc
            end if
            do icp = 1, nvi
                vinf(icp) = vind(icp)
            end do
            vinf(9) = 1.0d0
            iret = 9
            goto 999
        end if
    end if
! --  FIN   TRAITEMENT DE VENDOCHAB --
!
    b_n = to_blas_int(nvi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vind, b_incx, vinf, b_incy)
!
    if (rela_comp .eq. 'VENDOCHAB') then
!        INITIALISATION DE VINF(8) A UNE VALEUR NON NULLE
!        POUR EVITER LES 1/0 DANS RKDVEC
        if (vinf(8) .le. (1.0d-8)) then
            vinf(8) = 1.0d-8
        end if
    end if
!
!     COMPTAGE
!     irr=0
    decirr = 0
    nbsyst = 0
    decal = 0
    if (rela_comp .eq. 'MONOCRISTAL') then
        if (gdef .eq. 1) then
            if (opt .ne. 'RAPH_MECA') then
                ASSERT(ASTER_FALSE)
            end if
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vind(nvi-3-18+1), b_incx, fp, b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, id, b_incx, fp, &
                       b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, fp, b_incx, vinf(nvi-3-18+1), b_incy)
            nvi = nvi-9
        end if
        if (materf(nbcomm(1, 1), 2) .ge. 4) then
!           UNE SEULE FAMILLE
            ASSERT(nbcomm(nmat, 2) .eq. 1)
            necoul = cpmono(3) (1:16)
            irr = 0
            if (necoul .eq. 'MONO_DD_CC_IRRA') then
                irr = 1
                decirr = 6+3*12
            end if
            if (necoul .eq. 'MONO_DD_CFC_IRRA') then
                irr = 1
                decirr = 6+3*12
            end if
        end if
        nvi = nvi-3
!        NE PAS DIMINUER NVI DAVANTAGE A CAUSE DES GDEF
    end if
!      POUR POLYCRISTAL
!     INITIALISATION DE NBPHAS
    nbphas = nbcomm(1, 1)
    if (rela_comp .eq. 'POLYCRISTAL') then
!        RECUPERATION DU NOMBRE DE PHASES
        nbphas = nbcomm(1, 1)
        nsfv = 7+6*nbphas
        do iphas = 1, nbphas
            indpha = nbcomm(1+iphas, 1)
            nbfsys = nbcomm(indpha, 1)
            do ifa = 1, nbfsys
!              indice de la famille IFA
                indfa = indpha+ifa
                ifl = nbcomm(indfa, 1)
                nuecou = nint(materf(ifl, 2))
                nbsys = 12
                nsfv = nsfv+nbsys*3
            end do
        end do
        decirr = nsfv
    end if
!
!
999 continue
end subroutine
