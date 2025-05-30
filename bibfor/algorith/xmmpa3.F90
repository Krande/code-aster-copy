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
subroutine xmmpa3(ndim, nno, nnos, nnol, pla, &
                  ffc, ffp, jac, nfh, nd, &
                  cpenco, singu, fk, ddls, ddlm, &
                  jheavn, ncompn, nfiss, ifiss, jheafa, &
                  ncomph, ifa, mmat)
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, nnos, nnol
    integer :: nfh, ddls, ddlm
    integer :: singu, pla(27), nfiss, ifiss, jheafa, ncomph, ifa, jheavn, ncompn
    real(kind=8) :: mmat(216, 216), nd(3)
    real(kind=8) :: ffc(8), ffp(27), jac
    real(kind=8) :: cpenco
    real(kind=8) :: fk(27, 3, 3)
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES MATRICES A, AT, AU - CAS DU CONTACT
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
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  NOEUD  : INDICATEUR FORMULATION (T=NOEUDS , F=ARETE)
! IN  ND     : NORMALE À LA FACETTE ORIENTÉE DE ESCL -> MAIT
!                 AU POINT DE GAUSS
! IN  CPENCO : COEFFICIENT DE PENALISATION DU CONTACT
! IN  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! IN  RR     : DISTANCE AU FOND DE FISSURE
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  DDLM   : NOMBRE DE DDL A CHAQUE NOEUD MILIEU
! I/O MMAT   : MATRICE ELEMENTAITRE DE CONTACT/FROTTEMENT
!
!
    integer :: i, j, l, jn, jfh
    integer :: pli, plj, hea_fa(2)
    integer :: alpj
    real(kind=8) :: ffi, ffj, coefj
    aster_logical :: lmultc
!
! ----------------------------------------------------------------------
!
    coefj = xcalc_saut(1, 0, 1)
    lmultc = nfiss .gt. 1
    coefj = xcalc_saut(1, 0, 1, -88)
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    end if
! I.1 CALCUL DE A
    do i = 1, nnol
!
        pli = pla(i)
        ffi = ffc(i)
!
        do j = 1, nno
            call indent(j, ddls, ddlm, nnos, jn)
            do jfh = 1, nfh
                if (lmultc) then
                    coefj = xcalc_saut( &
                            zi(jheavn-1+ncompn*(j-1)+jfh), zi(jheafa-1+ncomph*(ifiss-1)+2*ifa-1), &
                            zi(jheafa-1+ncomph*(ifiss-1)+2*ifa), &
                            zi(jheavn-1+ncompn*(j-1)+ncompn) &
                            )
                else
                    coefj = xcalc_saut( &
                            zi(jheavn-1+ncompn*(j-1)+jfh), hea_fa(1), hea_fa(2), &
                            zi(jheavn-1+ncompn*(j-1)+ncompn) &
                            )
                end if
                do l = 1, ndim
                    mmat(pli, jn+ndim*jfh+l) = mmat(pli, jn+ndim*jfh+l)+coefj*ffi*ffp(j)*nd&
                                              &(l)*jac
!
! LBB : ON PREND AUSSI LE TERME EN PENALISATION
!
                    mmat(jn+ndim*jfh+l, pli) = mmat(jn+ndim*jfh+l, pli)+coefj*ffi*ffp(j)*nd&
                                              &(l)*jac
                end do
!
            end do
            do l = 1, singu*ndim
                do alpj = 1, ndim
                    mmat(pli, jn+ndim*(1+nfh)+alpj) = mmat( &
                                                      pli, &
                                                      jn+ndim*(1+nfh)+alpj)+2.d0*ffi*fk(j, &
                                                                                     alpj, l)*nd(l &
                                                                                               )*jac
!
                    mmat(jn+ndim*(1+nfh)+alpj, pli) = mmat(jn+ndim*(1+nfh)+alpj, &
                                                           pli)+2.d0*ffi*fk(j, alpj, l)*nd(l)*jac
                end do
            end do
!
        end do
!
    end do
!
!     CALCUL DE C
!
    do i = 1, nnol
!
        pli = pla(i)
        ffi = ffc(i)
!
        do j = 1, nnol
!
            plj = pla(j)
            ffj = ffc(j)
!
            mmat(pli, plj) = mmat(pli, plj)-ffj*ffi*jac/cpenco
!
        end do
    end do
!
end subroutine
