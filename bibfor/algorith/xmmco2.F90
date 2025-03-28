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
subroutine xmmco2(ndim, nno, nnos, nnol, ddls, &
                  ddlm, dsidep, p, r, nfh, &
                  jac, ffp, ffc, pla, singu, &
                  nfiss, jheafa, jheavn, ncompn, ifa, &
                  ncomph, ifiss, fk, mmat)
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/indent.h"
#include "asterfort/promat.h"
#include "asterfort/transp.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_saut.h"
    integer :: ndim, nno, nfh, ddls, singu, jheavn, ncompn
    real(kind=8) :: mmat(216, 216), dsidep(6, 6)
    real(kind=8) :: ffp(27), jac
    real(kind=8) :: p(3, 3)
    real(kind=8) :: fk(27, 3, 3)
!
! ROUTINE CONTACT (METHODE XFEM HPP - CALCUL ELEM.)
!
! --- CALCUL DES MATRICES DE COHESION
! ELEMENT COHESIF MIXTE
!
! ----------------------------------------------------------------------
!
! IN  NDIM   : DIMENSION DE L'ESPACE
! IN  NNO    : NOMBRE DE NOEUDS DE L'ELEMENT DE REF PARENT
! IN  DSIDEP : MATRICE TANGENTE BASE LOCALE
! IN  PP     :
! IN  P      : MATRICE PROJECTION PLAN TANGENT
! IN  ND     : DIRECTION NORMALE
! IN  NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  DDLS   : NOMBRE DE DDL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
! IN  JAC    : PRODUIT DU JACOBIEN ET DU POIDS
! IN  FFP    : FONCTIONS DE FORME DE L'ELEMENT PARENT
! IN  SINGU  : 1 SI ELEMENT SINGULIER, 0 SINON
! IN  TAU1   : PREMIERE DIRECTION TANGENTE
! IN  AM     :
! I/O MMAT   : MATRICE ELEMENTAITRE DE CONTACT/FROTTEMENT
!
!
    integer :: i, j, ddlm, ifa, ifh, ifiss, in, jfh, jheafa, jn, hea_fa(2)
    integer :: k, l, ncomph, nfiss, nnol, nnos, pla(27), pli, plj
    integer :: alpj, alpi
    real(kind=8) :: au(3, 3), coefi, coefj, dside2(3, 3), ffc(8), pdotal(3, 3)
    real(kind=8) :: ffi, ffj, r, temp(3, 3), unity(3, 3), ptr(3, 3)
    real(kind=8) :: alocal(3, 3)
    aster_logical :: lmultc
!
!
! ----------------------------------------------------------------------
!
!     INITIALISATIONS
!
    lmultc = nfiss .gt. 1
    unity(:, :) = 0.d0
    alocal(:, :) = 0.d0
    ptr(:, :) = 0.d0
    pdotal(:, :) = 0.d0
    au(:, :) = 0.d0
    dside2(:, :) = 0.d0
    temp(:, :) = 0.d0
!
!   MATRICE -ID+R DSIDEP
!
    do i = 1, ndim
        unity(i, i) = 1.d0
    end do
!
    do i = 1, ndim
        do j = 1, ndim
            dside2(i, j) = dsidep(i, j)
            alocal(i, j) = -unity(i, j)+r*dside2(i, j)
        end do
    end do
!
! MATRICE [P]T[ALOCAL]
!
    call transp(p, 3, ndim, ndim, ptr, &
                3)
!
    call promat(ptr, 3, ndim, ndim, alocal, &
                3, ndim, ndim, pdotal)
!
! MATRICE TANGENTE EN BASE FIXE [P]T [DSIDEP] [P]
!
    call promat(ptr, 3, ndim, ndim, alocal, &
                3, ndim, ndim, temp)
    call promat(temp, 3, ndim, ndim, p, &
                3, ndim, ndim, au)
!
! ON STOCKE DANS LA MATRICE ELEMENTAIRE
!
    coefi = xcalc_saut(1, 0, 1)
    coefj = xcalc_saut(1, 0, 1)
    if (.not. lmultc) then
        hea_fa(1) = xcalc_code(1, he_inte=[-1])
        hea_fa(2) = xcalc_code(1, he_inte=[+1])
    end if
!
    do i = 1, nnol
!
        pli = pla(i)
        ffi = ffc(i)
!
        do k = 1, ndim
!
            do j = 1, nno
                call indent(j, ddls, ddlm, nnos, jn)
                do jfh = 1, nfh
                    if (lmultc) then
                        coefj = xcalc_saut( &
                                zi(jheavn-1+ncompn*(j-1)+jfh), &
                                zi(jheafa-1+ncomph*(ifiss-1)+2*ifa-1), &
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
!
! INDICES INVERSES MATRICE INTERFACE
!
                        mmat(pli-1+k, jn+ndim*jfh+l) = mmat(pli-1+k, jn+ &
                                                            ndim*jfh+l)-coefj*ffi*ffp(j)*pdotal(l, &
                                                                                              k)*jac
!
! INDICES MEME ORDRE MATRICE EQUILIBRE
!
                        mmat(jn+ndim*jfh+l, pli-1+k) = mmat(jn+ndim*jfh+ &
                                                            l, pli-1+k)-coefj*ffi*ffp(j)*pdotal(l, &
                                                                                              k)*jac
                    end do
!
                end do
                do alpj = 1, singu*ndim
                    do l = 1, ndim
                        mmat(pli-1+k, jn+ndim*(1+nfh)+alpj) = mmat( &
                                                              pli-1+k, &
                                                              jn+ndim*(1+nfh)+alpj &
                                                              )-2.d0*ffi*fk(j, &
                                                                            alpj, l)*pdotal(l, k &
                                                                                            )*jac
!
                        mmat(jn+ndim*(1+nfh)+alpj, pli-1+k) = mmat(jn+ndim*(1+ &
                                                     nfh)+alpj, pli-1+k)-coefj*ffi*fk(j, alpj, l)* &
                                                              pdotal(l, k)*jac
                    end do
                end do
!
            end do
        end do
    end do
!
! -- MATRICE VENANT S AJOUTER A LA RAIDEUR
!
    do i = 1, nno
        call indent(i, ddls, ddlm, nnos, in)
        do j = 1, nno
            call indent(j, ddls, ddlm, nnos, jn)
            do ifh = 1, nfh
                if (lmultc) then
                    coefi = xcalc_saut( &
                            zi(jheavn-1+ncompn*(i-1)+ifh), zi(jheafa-1+ncomph*(ifiss-1)+2*ifa-1), &
                            zi(jheafa-1+ncomph*(ifiss-1)+2*ifa), &
                            zi(jheavn-1+ncompn*(i-1)+ncompn) &
                            )
                else
                    coefi = xcalc_saut( &
                            zi(jheavn-1+ncompn*(i-1)+ifh), hea_fa(1), hea_fa(2), &
                            zi(jheavn-1+ncompn*(i-1)+ncompn) &
                            )
                end if
                do jfh = 1, nfh
                    if (lmultc) then
                        coefj = xcalc_saut( &
                                zi(jheavn-1+ncompn*(j-1)+jfh), &
                                zi(jheafa-1+ncomph*(ifiss-1)+2*ifa-1), &
                                zi(jheafa-1+ncomph*(ifiss-1)+2*ifa), &
                                zi(jheavn-1+ncompn*(j-1)+ncompn) &
                                )
                    else
                        coefj = xcalc_saut( &
                                zi(jheavn-1+ncompn*(j-1)+jfh), hea_fa(1), hea_fa(2), &
                                zi(jheavn-1+ncompn*(j-1)+ncompn) &
                                )
                    end if
                    do k = 1, ndim
                        do l = 1, ndim
                            mmat(in+ndim*ifh+k, jn+ndim*jfh+l) = &
                                mmat(in+ndim*ifh+k, jn+ndim*jfh+l)- &
                                coefi*coefj*r*au(k, l)*ffp(i)*ffp(j)*jac
                        end do
!
                        do alpj = 1, singu*ndim
                            do l = 1, ndim
                                mmat(in+ndim+k, jn+ndim*(1+nfh)+alpj) = &
                                    mmat(in+ndim+k, jn+ndim*(1+nfh)+alpj)- &
                                    coefi*2.d0*ffp(i)*fk(j, alpj, l)*r*au(k, l)*jac
                            end do
                        end do
                    end do
                end do
            end do
!
            do alpi = 1, singu*ndim
                do k = 1, ndim
                    do jfh = 1, nfh
                        if (lmultc) then
                            coefj = xcalc_saut( &
                                    zi(jheavn-1+ncompn*(j-1)+jfh), &
                                    zi(jheafa-1+ncomph*(ifiss-1)+2*ifa-1), &
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
                            mmat(in+ndim*(1+nfh)+alpi, jn+ndim*jfh+l) = mmat( &
                                                                       in+ndim*(1+nfh)+alpi, &
                                                                       jn+ndim*jfh+l)-coefj*2.d&
                                                                       &0*fk(i, &
                                                                       alpi, k)*ffp(j)*r*au(k, l &
                                                                       )*jac
                        end do
                    end do
!
                    do alpj = 1, singu*ndim
                        do l = 1, ndim
                            mmat(in+ndim*(1+nfh)+alpi, jn+ndim*(1+nfh)+alpj) = &
                                mmat(in+ndim*(1+nfh)+alpi, jn+ndim*(1+nfh)+alpj)- &
                                4.d0*fk(i, alpi, k)*fk(j, alpj, l)*r*au(k, l)*jac
                        end do
                    end do
                end do
            end do
!
        end do
    end do
!
! -- MATRICE D INTERFACE : EXPRESSION DIRECTE
!
    do i = 1, nnol
!
        pli = pla(i)
        ffi = ffc(i)
        do k = 1, ndim
!
            do j = 1, nnol
!
                plj = pla(j)
                ffj = ffc(j)
                do l = 1, ndim
!
                    mmat(pli-1+k, plj-1+l) = mmat(pli-1+k, plj-1+l)-ffj*dside2(k, l)*ffi*jac
                end do
            end do
        end do
    end do
end subroutine
