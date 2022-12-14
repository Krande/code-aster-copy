! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine cabr2g(kpi, ipoids, ipoid2, ivf, ivf2,&
                  idfde, idfde2, geom, dimdef, dimuel,&
                  ndim, nddls, nddlm, nno, nnos,&
                  nnom, axi, regula, b, poids,&
                  poids2)
! aslint: disable=W1306,W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/utmess.h"
!
    aster_logical :: axi
    integer :: kpi, ipoids, ipoid2, idfde, idfde2, ndim, regula(6), dimdef, ivf
    integer :: ivf2, nno, nnos, nnom, nddls, nddlm, dimuel
    real(kind=8) :: geom(ndim, *), poids, poids2, b(dimdef, dimuel)
! ======================================================================
! --- BUT : CALCUL DE L'OPERATEUR B ------------------------------------
! ======================================================================
! --- IN ---------------------------------------------------------------
! --- NBDDL  : VECTEUR DIMENSION DU NOMBRE DE DDLS ---------------------
! --- NBNO   : VECTEUR DIMENSION DU NOMBRE DE NOEUDS -------------------
! --- KPI    : INDICE DU POINT D'INTEGRATION ---------------------------
! --- IPOIDS : ADRESSE DES FONCTIONS POIDS D'ORDRE 2 -------------------
! --- IPOID2 : ADRESSE DES FONCTIONS POIDS D'ORDRE 1 -------------------
! --- IVF2   : ADRESSE DES FONCTIONS DE FORME D'ORDRE 1 ----------------
! --- IDFDE  : ADRESSE DES DERIVEES DES FONCTIONS DE FORME D'ORDRE 2 ---
! --- IDFDE2 : ADRESSE DES DERIVEES DES FONCTIONS DE FORME D'ORDRE 1 ---
! --- GEOM   : CARACTERISTIQUES GEOMETRIQUES DE L'ELEMENT REEL ---------
! --- DIMDEF : DIMENSION DU VECTEUR DES DEFORMATIONS GENERALISEES ------
! --- NDIM   : DIMENSION DU PROBLEME -----------------------------------
! --- OUT --------------------------------------------------------------
! --- B      : OPERATEUR B DEFINI TEL QUE E=B.U ------------------------
! --- POIDS  : POIDS ASSOCIE AUX FONCTIONS DE FORME D'ORDRE 2 ----------
! --- POIDS2 : POIDS ASSOCIE AUX FONCTIONS DE FORME D'ORDRE 1 ----------
! ======================================================================
! ======================================================================
! --- VARIABLES LOCALES ------------------------------------------------
! ======================================================================
    integer :: i, j, k, n, adder1, adder2, adder3
    real(kind=8) :: dfdi(nno, 3), dfdi2(nnos, 3)
! ======================================================================
! --- INITIALISATION DE LA MATRICE B -----------------------------------
! ======================================================================
    b(:,:) = 0.d0
    adder1 = regula(1)
    adder2 = regula(2)
    adder3 = regula(3)
! ======================================================================
! --- CAS 2D -----------------------------------------------------------
! ======================================================================
    if (ndim .eq. 2) then
! ======================================================================
! --- CAS QUADRATIQUES -------------------------------------------------
! ======================================================================
        call dfdm2d(nno, kpi, ipoids, idfde, geom,&
                    poids, dfdi(1, 1), dfdi(1, 2))
! ======================================================================
! --- CAS LINEAIRES ----------------------------------------------------
! ======================================================================
        call dfdm2d(nnos, kpi, ipoid2, idfde2, geom,&
                    poids2, dfdi2(1, 1), dfdi2(1, 2))
    else if (ndim.eq.3) then
! ======================================================================
! --- CAS QUADRATIQUES -------------------------------------------------
! ======================================================================
        call dfdm3d(nno, kpi, ipoids, idfde, geom,&
                    poids, dfdi(1, 1), dfdi(1, 2), dfdi(1, 3))
! ======================================================================
! --- CAS LINEAIRES ----------------------------------------------------
! ======================================================================
        call dfdm3d(nnos, kpi, ipoid2, idfde2, geom,&
                    poids2, dfdi2(1, 1), dfdi2(1, 2), dfdi2(1, 3))
    else
        call utmess('F', 'ALGORITH6_13')
    endif
! ======================================================================
! --- REMPLISSAGE DE L OPERATEUR B -------------------------------------
! ======================================================================
! --- TERMES -(DUJ/DXI-VJI) --------------------------------------------
! ======================================================================
! --- SUR LES NOEUDS SOMMETS -------------------------------------------
! ======================================================================
    do n = 1, nnos
        do j = 1, ndim
            do i = 1, ndim
                b(adder1-1+(j-1)*ndim+i,(n-1)*nddls+j) = b(&
                                                         adder1-1+(j-1)*ndim+i,&
                                                         (n-1)*nddls+j ) - dfdi(n, i&
                                                         )
                b(adder1-1+(j-1)*ndim+i,(n-1)*nddls+ndim+(j-1)*ndim+i)&
                = b(adder1-1+(j-1)*ndim+i,(n-1)*nddls+ndim+(j-1)*ndim+&
                i) + zr(ivf2+n+(kpi-1)*nnos-1)
            end do
        end do
    end do
    do n = 1, nnom
        do j = 1, ndim
            do i = 1, ndim
                b(adder1-1+(j-1)*ndim+i,nnos*nddls+(n-1)*nddlm+j) =&
                b(adder1-1+(j-1)*ndim+i,nnos*nddls+(n-1)*nddlm+j) -&
                dfdi(n+nnos,i)
            end do
        end do
    end do
! ======================================================================
! --- POUR LES GRADIENTS DE VARIATIONS VOLUMIQUE -----------------------
! --- ON UTILISE LES FONCTIONS DE FORME D'ORDRE 1 ----------------------
! ======================================================================
! --- SUR LES NOEUDS SOMMETS -------------------------------------------
! ======================================================================
    do n = 1, nnos
        do k = 1, ndim
            do j = 1, ndim
                do i = 1, ndim
                    b(adder2-1+(k-1)*ndim*ndim+(j-1)*ndim+i, (n-1)*&
                    nddls+ndim+(k-1)*ndim+j)= b(adder2-1+(k-1)*ndim*&
                    ndim+(j-1)*ndim+i, (n-1)*nddls+ndim+(k-1)*ndim+j)+&
                    dfdi2(n,i)
                end do
            end do
        end do
    end do
! ======================================================================
! --- POUR LE MULTIPLICATEUR DE LAGRANGE -------------------------------
! --- (PRES) -----------------------------------------------------------
! --- ON UTILISE LES FONCTIONS DE FORME D'ORDRE 0 ----------------------
! ======================================================================
! --- SUR LES NOEUDS CENTRAUX ------------------------------------------
! ======================================================================
    do i = 1, ndim
        do j = 1, ndim
            b(adder3-1+(i-1)*ndim+j,nnos*nddls+nnom*nddlm+(i-1)*ndim+&
            j)= 1.0d0
        end do
    end do
! ======================================================================
end subroutine
