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
!
subroutine xmmsa1(algofr, ndim, nno, nnos, nnol, &
                  pla, ffc, ffp, idepd, idepm, &
                  nfh, nd, tau1, tau2, singu, &
                  fk, lact, ddls, ddlm, coeffr, &
                  coeffp, p, adher, knp, ptknp, &
                  ik)
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/promat.h"
#include "asterfort/xadher.h"
#include "asterfort/xcalc_saut.h"
#include "asterfort/xmafr1.h"
    integer :: algofr, ndim, nno, nnos, nnol
    integer :: nfh, ddls, ddlm
    integer :: singu, pla(27), lact(8), idepd, idepm
    real(kind=8) :: coeffr, coeffp, p(3, 3), ik(3, 3)
    real(kind=8) :: ffc(8), ffp(27), tau1(3), tau2(3), ptknp(3, 3)
    real(kind=8) :: knp(3, 3), nd(3)
    aster_logical :: adher
    real(kind=8) :: fk(27, 3, 3)
!
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES INCREMENTS - DÉPLACEMENTS&
! --- SEMI-MULTIPLICATEUR DE FROTTEMENT
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
! IN  IDEPD  :
! IN  IDEPM  :
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  NOEUD  : INDICATEUR FORMULATION (T=NOEUDS , F=ARETE)
! IN  TAU1   : TANGENTE A LA FACETTE AU POINT DE GAUSS
! IN  TAU2   : TANGENTE A LA FACETTE AU POINT DE GAUSS
! IN  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! IN  IFA    : INDICE DE LA FACETTE COURANTE
! IN  CFACE  : CONNECTIVITÉ DES NOEUDS DES FACETTES
! IN  LACT   : LISTE DES LAGRANGES ACTIFS
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  DDLM   : NOMBRE DE DDL A CHAQUE NOEUD MILIEU
! IN  COEFFR  :
! IN  COEFFR : COEFFICIENTS DE STABILISATION DU FROTTEMENT
! IN  COEFFP : COEFFICIENTS DE PENALISATION DU FROTTEMENT
! IN  LPENAF : INDICATEUR DE PENALISATION DU FROTTEMENT
! IN  P      :
! OUT ADHER  :
! OUT KNP    : PRODUIT KN.P
! OUT PTKNP  : MATRICE PT.KN.P
! OUT IK     :
!
!
!
!
    integer :: i, j, nli, ino, in, pli, ig
    integer :: alpi
    real(kind=8) :: ffi, lamb1(3), r3(3), vitang(3), kn(3, 3), saut(3), coefj
!
! ----------------------------------------------------------------------
!
! --- INITIALISATION
    saut(:) = 0.d0
    lamb1(:) = 0.d0
    ptknp(:, :) = 0.d0
    p(:, :) = 0.d0
    knp(:, :) = 0.d0
    kn(:, :) = 0.d0
!
    coefj = xcalc_saut(1, 0, 1)
!
    call xmafr1(ndim, nd, p)
!
    do ino = 1, nno
        call indent(ino, ddls, ddlm, nnos, in)
!
        do j = 1, ndim
            do ig = 1, nfh
                saut(j) = saut(j)-coefj*ffp(ino)*zr(idepd-1+in+ndim*(1+ig-1)+j)
            end do
        end do
        do j = 1, ndim*singu
            do alpi = 1, ndim
                saut(j) = saut(j)-2.d0*fk(ino, alpi, j)*zr(idepd-1+in+ndim*(1+nfh)+alpi)
            end do
        end do
    end do
!
    do i = 1, nnol
        pli = pla(i)
        ffi = ffc(i)
        nli = lact(i)
        if (nli .eq. 0) goto 158
!
        do j = 1, ndim
            lamb1(j) = lamb1(j)+ffi*tau1(j)*(zr(idepd-1+pli+1)+zr( &
                                             idepm-1+pli+1))
!
            if (ndim .eq. 3) lamb1(j) = lamb1(j)+ffi*tau2(j)*(zr(idepd-1+pli+2)+zr(idepm-1+&
                             &pli+2))
        end do
158     continue
    end do
!
!
! --- TEST DE L'ADHERENCE ET CALCUL DES MATRICES DE FROTTEMENT UTILES
!
    call xadher(p, saut, lamb1, coeffr, coeffp, &
                algofr, vitang, r3, kn, ptknp, &
                ik, adher)
!
    call promat(kn, 3, ndim, ndim, p, &
                3, ndim, ndim, knp)
!
end subroutine
