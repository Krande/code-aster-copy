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
subroutine nmholi(ndim, axi, nno, npg, ipoids, &
                  ivf, idfde, imate, inst, geom, &
                  depl, chlim)
    implicit none
#include "asterf_types.h"
#include "asterc/matfpe.h"
#include "asterfort/nmgeom.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "blas/dnrm2.h"
    aster_logical :: axi
    integer(kind=8) :: ndim, nno, npg, imate, ipoids, ivf, idfde
    real(kind=8) :: geom(ndim, nno), depl(ndim, nno), inst, chlim(3)
!
! --------------------------------------------------------------------
!        CALCUL DES TERMES POUR LE POST TRAITEMENT CHARGE_LIMITE
! -------------------------------------------------------------------
! IN  NDIM   DIMENSION
! IN  AXI    .TRUE. SI AXISYMETRIQUE
! IN  NNO    NOMBRE DE NOEUDS PORTANT LE DEPLACEMENT
! IN  NPG    NOMBRE DE POINTS DE GAUSS DE MECANIQUE
! IN  VFF    VALEUR DES FOCNTIONS DE FORME
! IN  DFDE   DERIVEES DES FONCTIONS DE FORME (REFERENCE)
! IN  DFDN   DERIVEES DES FONCTIONS DE FORME (REFERENCE)
! IN  DFDK   DERIVEES DES FONCTIONS DE FORME (REFERENCE)
! IN  POIDSG POIDS DES POINTS DE GAUSS       (REFERENCE)
! IN  IMATE  ADRESSE DU MATERIAU
! IN  INST   INSTANT COURANT
! IN  GEOM   COORDONNEES DES NOEUDS
! IN  DEPL   DEPLACEMENTS NODAUX
! OUT CHLIM  TERMES CALCULES :
!             1 - SOMME( SY *EPSEQ )
!             2 - SOMME( A(M)/M * EPSNO**M )
!             3 - MAX( SIEQ/SY )
! -------------------------------------------------------------------
!
    integer(kind=8) :: kpg, ndimsi, spt
    character(len=8) :: fami, poum
    real(kind=8) :: eps(6), poids, epsno, sy, m, am, epsh
    real(kind=8) :: dfdi(27, 3), fbid(3, 3), r
    real(kind=8) :: rac23, val(1)
    integer(kind=8) :: cod(1)
    blas_int :: b_incx, b_n
! ------------------------------------------------------------------
!
!
    call matfpe(-1)
!
! -- INITIALISATION
!
    ndimsi = 2*ndim
    rac23 = sqrt(2.d0/3.d0)
    call r8inir(3, 0.d0, chlim, 1)
!
!
! -- CARACTERISTIQUES
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
    call rcvalb(fami, kpg, spt, poum, imate, &
                ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                1, 'SY', val, cod, 2)
    sy = val(1)
    m = 1+10**(1-inst)
    am = sy*rac23**m
!
    do kpg = 1, npg
!
! -- DEFORMATION
!
        call nmgeom(ndim, nno, axi, .false._1, geom, &
                    kpg, ipoids, ivf, idfde, depl, &
                    .true._1, poids, dfdi, fbid, eps, &
                    r)
        epsh = (eps(1)+eps(2)+eps(3))/3
        eps(1) = eps(1)-epsh
        eps(2) = eps(2)-epsh
        eps(3) = eps(3)-epsh
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        epsno = dnrm2(b_n, eps, b_incx)
!
! - CALCUL DES TERME ELEMENTAIRES
!
        chlim(1) = chlim(1)+poids*sy*rac23*epsno
        chlim(2) = chlim(2)+poids*am/m*epsno**m
        chlim(3) = max(chlim(3), (rac23*epsno)**(m-1))
!
    end do
    call matfpe(1)
!
end subroutine
