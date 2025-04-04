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
subroutine xmvec2(ndim, nno, nnos, nnol, pla, &
                  ffc, ffp, reac, jac, nfh, &
                  saut, singu, fk, nd, coeffr, &
                  ddls, ddlm, jheavn, ncompn, nfiss, &
                  ifiss, jheafa, ncomph, ifa, vtmp)
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, nnos, nnol
    integer :: pla(27), nfh
    integer :: singu, ddls, ddlm, nfiss, ifiss, jheafa, ncomph, ifa, jheavn, ncompn
    real(kind=8) :: vtmp(400), coeffr, saut(3), nd(3)
    real(kind=8) :: ffc(8), ffp(27), jac, reac
    real(kind=8) :: fk(27, 3, 3)
!
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DU VECTEUR LN1 & LN2
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
! IN  RHOTK  :
! IN  COEFFR : COEFFICIENTS DE AUGMENTATION DU FROTTEMENT
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
    integer :: i, j, in, pli, ifh, hea_fa(2)
    integer :: alpi
    real(kind=8) :: ffi, dn, coefi
    aster_logical :: lmultc
!
! ----------------------------------------------------------------------
!
    coefi = xcalc_saut(1, 0, 1)
    lmultc = nfiss .gt. 1
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    end if
    dn = 0.d0
    do j = 1, ndim
        dn = dn+saut(j)*nd(j)
    end do
!
! --- TERME LN1
!
    do i = 1, nno
        call indent(i, ddls, ddlm, nnos, in)
        do ifh = 1, nfh
            if (lmultc) then
                coefi = xcalc_saut( &
                        zi(jheavn-1+ncompn*(i-1)+ifh), zi(jheafa-1+ncomph*(ifiss-1)+2*ifa-1), &
                        zi(jheafa-1+ncomph*(ifiss-1)+2*ifa), zi(jheavn-1+ncompn*(i-1)+ncompn) &
                        )
            else
                coefi = xcalc_saut( &
                        zi(jheavn-1+ncompn*(i-1)+ifh), hea_fa(1), hea_fa(2), &
                        zi(jheavn-1+ncompn*(i-1)+ncompn) &
                        )
            end if
            do j = 1, ndim
                vtmp(in+ndim*ifh+j) = vtmp(in+ndim*ifh+j)+(reac-coeffr*dn)*coefi*ffp(i)*nd(j)*&
                                      &jac
            end do
        end do
        do j = 1, singu*ndim
            do alpi = 1, ndim
                vtmp(in+ndim*(1+nfh)+alpi) = vtmp( &
                                             in+ndim*(1+nfh)+alpi)+(reac-coeffr*dn)*2.d0*fk(i, &
                                                                                     alpi, j)*nd(j &
                                                                                               )*jac
            end do
        end do
    end do
!
! --- TERME LN2
!
    do i = 1, nnol
        pli = pla(i)
        ffi = ffc(i)
!
        vtmp(pli) = vtmp(pli)-dn*ffi*jac
!
    end do
end subroutine
