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

subroutine lckimp(ndim, typmod, option, mat, eps, &
                  phitot, vim, sig, forc_endo, vip, dsde_1, dsde_2, dsde_3)

    implicit none

#include "asterf_types.h"
#include "asterfort/rcvala.h"

    integer(kind=8)           :: ndim
    character(len=8)  :: typmod
    character(len=16) :: option
    integer(kind=8)           :: mat
    real(kind=8)      :: eps(:)
    real(kind=8)      :: phitot
    real(kind=8)      :: vim(:)
    real(kind=8)      :: vip(:)
    real(kind=8)      :: sig(:)
    real(kind=8)      :: forc_endo
    real(kind=8)      :: dsde_1(:, :)
    real(kind=8)      :: dsde_2(:)
    real(kind=8)      :: dsde_3
!
!     -----------------------------------------------------------------
!     ENDOMMAGEMENT FRAGILE ENDO_CARRE POUR GVNO
!     -----------------------------------------------------------------
!     IN  NDIM    DIMENSION DE L'ESPACE
!     IN  TYPMOD  TYPE DE MODELISATION
!     IN  OPTION  OPTION DE CALCUL
!     RIGI_MECA_TANG, RIGI_MECA_ELAS
!     RAPH_MECA
!     FULL_MECA, FULL_MECA_ELAS
!     IN  MAT     NATURE DU MATERIAU
!     IN  EPSM    CHAMP DE DEFORMATION EN T- ET PHIM=EPSM(7)
!     IN  DEPS    INCREMENT DU CHAMP DE DEFORMATION ET DPHI=DEPS(7)
!     IN  VIM     VARIABLES INTERNES EN T-
!     OUT VIP     DENSITE DE FISSURATION
!     OUT SIG     CONTRAINTE
!     OUT DSIDEP  MATRICES TANGENTES
!     -----------------------------------------------------------------
!
    aster_logical :: cplan, resi
    integer(kind=8) :: ndimsi, ij, kl
    real(kind=8) :: val(3), nu, lambda, deuxmu
    real(kind=8) :: e, phi
    real(kind=8) :: coplan, w, treps, epseps, sigel(2*ndim)
    real(kind=8) :: fd
    real(kind=8) :: kk, epsd(2*ndim)
    real(kind=8) :: sigm, w0
    integer(kind=8) :: k2(4)
    character(len=8) :: nom(3)
    real(kind=8) :: kron(6), rigmin, kr(2*ndim)
    parameter(rigmin=1.d-5)
    data kron/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/
!     -----------------------------------------------------------------
!
!
!     -----------------------------------------------------------------
!     INITIALISATIONS
!     -----------------------------------------------------------------
!
!
!     -- OPTIONS DE CALCUL
!
    cplan = typmod .eq. 'C_PLAN  '
    resi = option(1:9) .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA'
!      RIGI=OPTION(1:9).EQ.'RIGI_MECA'.OR.OPTION(1:9).EQ.'FULL_MECA'
!      ELAS=OPTION.EQ.'RIGI_MECA_ELAS'.OR.OPTION.EQ.'FULL_MECA_ELAS'
    ndimsi = 2*ndim
    kr = kron(1:ndimsi)

    phi = min(phitot, 1.d0)
!
!
!     -- LECTURE DES CARACTERISTIQUES MATERIAU
!
    nom(1) = 'E'
    nom(2) = 'NU'
    nom(3) = 'SY'
!
    call rcvala(mat, ' ', 'ELAS', 0, ' ', &
                [0.d0], 2, nom(1), val(1), k2, &
                2)
!
    nu = val(2)
    e = val(1)
    lambda = val(1)*val(2)/(1-2*val(2))/(1+val(2))
    deuxmu = val(1)/(1.d0+val(2))
    kk = lambda+deuxmu/(3.0d0)
!
    call rcvala(mat, ' ', 'ECRO_LINE', 0, ' ', &
                [0.d0], 1, nom(3), val(3), k2, &
                2)
!
    sigm = val(3)
!
    w0 = sigm**2/(2*e)

!
!     DEFORMATION HORS PLAN POUR LES CONTRAINTES PLANES
    if (cplan) then
        coplan = -nu/(1.d0-nu)
        eps(3) = coplan*(eps(1)+eps(2))
    end if
!
!
!     -- ENERGIE DE DEFORMATION ET CONTRAINTE ELASTIQUE
!
    treps = eps(1)+eps(2)+eps(3)
!
!     -- DEVIATEUR DES DEFORMATIONS
!
    epsd = eps-treps*kr/3.0d0
!
    epseps = dot_product(eps, eps)
    w = 0.5d0*(lambda*treps**2+deuxmu*epseps)
    sigel = lambda*treps*kr+deuxmu*eps
!
!     CORRECTION 1 DE LA DERIVEE PAR RAPPORT A D EN COMPRESSION
!
    if (treps .lt. 0.d0) then
        w = 0.5d0*deuxmu*dot_product(epsd, epsd)
    end if
!
!     -----------------------------------------------------------------
!     CALCUL DE L'ENDOMMAGEMENT
!     -----------------------------------------------------------------
!
    fd = (1.d0-phi)**2+rigmin
!
    if (.not. resi) goto 500
!
    vip(1) = phi
    vip(2) = merge(1.d0, 0.d0, vip(1) .gt. vim(1))
!
!
!     STOCKAGE DES CONTRAINTES ET DES VARIABLES INTERNES
!
!
!     FORMULATION LOI DE COMPORTMENT AVEC CORRECTION EN COMPRESSION
!
    if (treps .lt. 0.d0) then
        sig = kk*treps*kr+deuxmu*epsd*fd
    else
        sig = sigel*fd
    end if
!
500 continue
!
!     -----------------------------------------------------------------
!     CALCUL DES MATRICES TANGENTES
!     -----------------------------------------------------------------
!
!
    dsde_1 = 0
    dsde_2 = 0
    dsde_3 = 0

    fd = (1.d0-phi)**2+rigmin
!
!     -- CONTRIBUTION ELASTIQUE
!
    if (treps .lt. 0.d0) then
        dsde_1(1:3, 1:3) = lambda+deuxmu/(3.0d0)*(1.d0-fd)
    else
        dsde_1(1:3, 1:3) = fd*lambda
    end if
    do ij = 1, ndimsi
        dsde_1(ij, ij) = dsde_1(ij, ij)+fd*deuxmu
    end do
!
!     -- CORRECTION POUR LES CONTRAINTES PLANES
!
    if (cplan) then
        do 130 ij = 1, ndimsi
            if (ij .eq. 3) goto 130
            do 140 kl = 1, ndimsi
                if (kl .eq. 3) goto 140
                dsde_1(ij, kl) = dsde_1(ij, kl)-1.d0/dsde_1(3, 3)*dsde_1(ij, 3)*dsde_1(3, kl)
140             continue
130             continue
                end if
!
!     -- CORRECTION DISSIPATIVE
!
!     CORRECTION 2 DE LA DERIVEE PAR RAPPORT A D EN COMPRESSION
!
                if (treps .lt. 0.d0) then
                    sigel = deuxmu*epsd
                end if
!
!     DERIVEES CROISEES
!
                dsde_2 = -2.d0*(1.0d0-phi)*sigel
!
!     DERIVEE SECONDE /ENDO
!
                dsde_3 = 2.0d0*w
!
!    DERIVEE PREMIERE /ENDO
!
                forc_endo = 2.0d0*(w0-(1.0d0-phi)*w)
!
!
!
                end subroutine
