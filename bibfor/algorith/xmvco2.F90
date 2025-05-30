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
subroutine xmvco2(ndim, nno, nnol, nnos, lamb, &
                  am, delta, pla, lact, nfh, &
                  ddls, ddlm, nfiss, ifiss, jheafa, &
                  ifa, ncomph, jheavn, ncompn, jac, &
                  ffc, ffp, singu, r, fk, &
                  vtmp, p)
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/prmave.h"
#include "asterfort/transp.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, nnol, jheavn, ncompn
    integer :: nfh, ddls, pla(27), lact(8)
    integer :: singu
    real(kind=8) :: vtmp(400), delta(6)
    real(kind=8) :: ffp(27), jac
    real(kind=8) :: ffc(8)
    real(kind=8) :: fk(27, 3, 3)
!
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES SECONDS MEMBRES DE COHESION
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
! IN  SIGMA  : VECTEUR CONTRAINTE EN REPERE LOCAL
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! IN  RR     : DISTANCE AU FOND DE FISSURE
! I/O VTMP   : VECTEUR ELEMENTAIRE DE CONTACT/FROTTEMENT
!
!
    integer :: i, j, k, pli, nli, ddlm, ier, ifa, ifh, ifiss
    integer :: in, jheafa, ncomph, nfiss, nnos, hea_fa(2), alp
    real(kind=8) :: ffi, am(3), coefi, hfix(3), h(3), lamb(3), r
    real(kind=8) :: p(3, 3), ptr(3, 3)
    aster_logical :: lmultc
!
! ---------------------------------------------------------------------
!
! INITIALISATIONS
!
    lmultc = nfiss .gt. 1
    h(:) = 0.d0
    hfix(:) = 0.d0
    ptr(:, :) = 0.d0
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    end if
!
! CALCUL DE H = R*DELTA - FORCE COHESIVE AUGMENTEE EN BASE COVARIANTE
! RAPPEL : AM INVERSE PAR RAPPORT AUX CONVENTIONS X-FEM
!
    do i = 1, ndim
        h(i) = -lamb(i)-r*am(i)+r*delta(i)
    end do
!
! CONVERSION DE H EN BASE FIXE : {HFIX} = [P]T {H}
! RAPPEL : P MATRICE DE PASSAGE BASE FIXE --> BASE COVARIANTE
!
    call transp(p, 3, ndim, ndim, ptr, &
                3)
    call prmave(0, ptr, 3, ndim, ndim, &
                h, ndim, hfix, ndim, ier)
!
! ON STOCKE DANS LE VECTEUR SECOND MEMBRE ELEMENTAIRE DE L EQUILIBRE
! ! IL Y A DEJA UN MOINS DU AUX CONVENTIONS POUR LE SAUT
!
    coefi = xcalc_saut(1, 0, 1)
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
                vtmp(in+ndim*ifh+j) = vtmp(in+ndim*ifh+j)-coefi*ffp(i)*hfix(j)*jac
            end do
        end do
        do alp = 1, singu*ndim
            do j = 1, ndim
                vtmp(in+ndim*(1+nfh)+alp) = vtmp(in+ndim*(1+nfh)+alp)-2.d0*fk(i, alp, j)*hfix(j &
                                                                                              )*jac
            end do
        end do
    end do
!
! SECOND MEMBRE DE L EQUATION D INTERFACE: EXPRESSION DIRECTE
! ATTENTION INVERSION DE CONVENTIONS
!
    do i = 1, nnol
        pli = pla(i)
        ffi = ffc(i)
        nli = lact(i)
        if (nli .eq. 0) goto 20
        do k = 1, ndim
            vtmp(pli-1+k) = vtmp(pli-1+k)+(am(k)-delta(k))*ffi*jac
        end do
20      continue
    end do
!
end subroutine
