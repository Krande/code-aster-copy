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
subroutine xmvec3(nnol, pla, ffc, reac, jac, &
                  cstaco, vtmp)
!
    implicit none
    integer :: nnol
    integer :: pla(27)
    real(kind=8) :: vtmp(400), cstaco
    real(kind=8) :: ffc(8), jac, reac
!
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DU VECTEUR LN2
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
! IN  RR     : DISTANCE AU FOND DE FISSURE
! IN  IFA    : INDICE DE LA FACETTE COURANTE
! IN  CFACE  : CONNECTIVITÉ DES NOEUDS DES FACETTES
! IN  LACT   : LISTE DES LAGRANGES ACTIFS
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  DDLM   : NOMBRE DE DDL A CHAQUE NOEUD MILIEU
! IN  RHOTK  :
! IN  CSTAFR : COEFFICIENTS DE STABILISATION DU FROTTEMENT
! IN  CPENFR : COEFFICIENTS DE PENALISATION DU FROTTEMENT
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
!
    integer :: i, pli
    real(kind=8) :: ffi
!
! ----------------------------------------------------------------------
!
    do i = 1, nnol
!
        pli = pla(i)
        ffi = ffc(i)
!
        vtmp(pli) = vtmp(pli)-reac*ffi*jac/cstaco
    end do
!
end subroutine
