! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine xmmbp3(ndim, nno, nnos, nnol, pla, &
                  ffc, ffp, jac, knp, nfh, &
                  seuil, tau1, tau2, mu, singu, &
                  fk, lact, ddls, ddlm, mmat)
!
    implicit none
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, nnos, nnol
    integer :: nfh, ddls, ddlm
    integer :: singu, pla(27), lact(8)
    real(kind=8) :: mmat(216, 216)
    real(kind=8) :: ffc(8), ffp(27), jac, tau1(3), tau2(3)
    real(kind=8) :: seuil, knp(3, 3), mu
    real(kind=8) :: fk(27, 3, 3)
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DE B, BT
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
! IN  NNOS   : NOMBRE DE NOEUDS SOMMET DE L'ELEMENT DE REF PARENT
! IN  NNOL   : NOMBRE DE NOEUDS PORTEURS DE DDLC
! IN  NNOF   : NOMBRE DE NOEUDS DE LA FACETTE DE CONTACT
! IN  PLA    : PLACE DES LAMBDAS DANS LA MATRICE
! IN  IPGF   : NUMÉRO DU POINTS DE GAUSS
! IN  IVFF   : ADRESSE DANS ZR DU TABLEAU FF(INO,IPG)
! IN  FFC    : FONCTIONS DE FORME DE L'ELEMENT DE CONTACT
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  KNP    : PRODUIT KN.P
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  NOEUD  : INDICATEUR FORMULATION (T=NOEUDS , F=ARETE)
! IN  SEUIL  : SEUIL
! IN  TAU1   : TANGENTE A LA FACETTE AU POINT DE GAUSS
! IN  TAU2   : TANGENTE A LA FACETTE AU POINT DE GAUSS
! IN  MU     : COEFFICIENT DE COULOMB
! IN  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! IN  RR     : DISTANCE AU FOND DE FISSURE
! IN  IFA    : INDICE DE LA FACETTE COURANTE
! IN  CFACE  : CONNECTIVITÉ DES NOEUDS DES FACETTES
! IN  LACT   : LISTE DES LAGRANGES ACTIFS
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  DDLM   : NOMBRE DE DDL A CHAQUE NOEUD MILIEU
! IN  LPENAF : INDICATEUR DE PENALISATION DU FROTTEMENT
! I/O MMAT   : MATRICE ELEMENTAITRE DE CONTACT/FROTTEMENT
!
!
    integer :: i, j, k, l, jn, nli, pli
    integer :: alpj
    real(kind=8) :: ffi, tauknp(2, 3), coefj
!
! ----------------------------------------------------------------------
!
    tauknp(:, :) = 0.d0
    coefj = xcalc_saut(1, 0, 1)
!
!     II.3.1. CALCUL DE B ET DE BT
!
    do i = 1, nnol
        pli = pla(i)
        ffi = ffc(i)
        nli = lact(i)
        if (nli .eq. 0) cycle
!
!     CALCUL DE TAU.KN.P
        do j = 1, ndim
            tauknp(1, j) = 0.d0
            do k = 1, ndim
                tauknp(1, j) = tauknp(1, j)+tau1(k)*knp(k, j)
            end do
        end do
!
        if (ndim .eq. 3) then
            do j = 1, ndim
                tauknp(2, j) = 0.d0
                do k = 1, ndim
                    tauknp(2, j) = tauknp(2, j)+tau2(k)*knp(k, j)
                end do
            end do
        end if
!
        do j = 1, nno
            call indent(j, ddls, ddlm, nnos, jn)
            do k = 1, ndim-1
                do l = 1, nfh*ndim
                    mmat(pli+k, jn+ndim+l) = mmat(pli+k, jn+ndim+l)+ &
                                             coefj*mu*seuil*ffi*ffp(j)*tauknp(k, l)*jac
                end do
                do l = 1, singu*ndim
                    do alpj = 1, ndim
                        mmat(pli+k, jn+ndim*(1+nfh)+alpj) = mmat(pli+k, jn+ndim*(1+nfh)+alpj)+ &
                                                   2.d0*fk(j, alpj, l)*mu*seuil*ffi*tauknp(k, l)*jac
                    end do
                end do
            end do
        end do
    end do
!
end subroutine
