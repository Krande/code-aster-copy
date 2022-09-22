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
subroutine xadher(p, saut, lamb1, cstafr, cpenfr,&
                  algofr, vitang, pboul, kn, ptknp,&
                  ik, adher)
!
! person_in_charge: samuel.geniaut at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/matini.h"
#include "asterfort/utbtab.h"
    integer :: algofr
    real(kind=8) :: p(3, 3), saut(3), lamb1(3), cstafr, cpenfr, pboul(3)
    real(kind=8) :: vitang(3), ptknp(3, 3), ik(3, 3), kn(3, 3)
    aster_logical :: adher
!
! ----------------------------------------------------------------------
!                      TEST DE L'ADHÉRENCE AVEC X-FEM
!                ET CALCUL DES MATRICES DE FROTTEMENT UTILES
!
! IN    P           : OPÉRATEUR DE PROJECTION
! IN    SAUT        : SAUT DES INCRÉMENTS DE DÉPLACEMENTS
!                     DEPUIS L'ÉQUILIBRE PRÉCÉDENT
! IN    LAMB1       : INCRÉMENTS DU SEMI-MULTIPLICATEUR DE FROTTEMENT
!                     DEPUIS L'ÉQUILIBRE PRÉCÉDENT DANS LA BASE GLOBALE
! IN    CSTAFR       : COEFFICIENT DE REGULARISATION DE FROTTEMENT
!                       DIVISÉ PAR L'INCRÉMENT DE TEMPS
! IN    CSTAFR      : COEFFICIENT DE STABILISTAION DE FROTTEMENT
!                       DIVISÉ PAR L'INCRÉMENT DE TEMPS
! IN    CPENFR      : COEFFICIENT DE PENALISATION DU FROTTEMENT
!                       DIVISÉ PAR L'INCRÉMENT DE TEMPS
! IN    ALGOFR      : FORMULATION POUR LE FROTTEMENT
!
! OUT   VITANG      : PROJECTION TANGENTE DU SAUT
! OUT   PBOUL       : PROJECTION SUR LA BOULE B(0,1)
! OUT   KN          : KN(LAMDBA + CSTAFR [[DX]]/DELTAT )
! OUT   PTKNP       : MATRICE Pt.KN.P
! OUT   IK          : MATRICE Id-KN
! OUT   ADHER       : STATUT D'ADHERENCE
!
!
!
!
    integer :: ndim, i, j, k
    real(kind=8) :: prec, norme, xab(3, 3), gt(3)
    real(kind=8) :: p2(2, 2), ptknp2(2, 2), kn2(2, 2), xab2(2, 2)
    real(kind=8) :: gt2(3), norme2
    aster_logical :: lpenaf
    parameter  (prec=1.d-12)
!
!-----------------------------------------------------------------------
!     CALCUL DE GT = LAMDBA + RHO [[DX]]/DELTAT ET DE SA PROJECTION
!
    call elrefe_info(fami='RIGI', ndim=ndim)
    lpenaf = (algofr.eq.2)
    do i = 1, ndim
        vitang(i)=0.d0
!       "VITESSE TANGENTE" : PROJECTION DU SAUT
        do k = 1, ndim
            vitang(i)=vitang(i)+p(i,k)*saut(k)
        end do
        if (lpenaf) then
!         PENALISATION SEULE
            gt(i)=cpenfr * vitang(i)
        else
            gt(i)=lamb1(i)+cstafr*vitang(i)
        endif
    end do
    if (ndim .eq. 3) then
        norme=sqrt(gt(1)*gt(1)+gt(2)*gt(2)+gt(3)*gt(3))
    else
        norme=sqrt(gt(1)*gt(1)+gt(2)*gt(2))
    endif
!
    if (lpenaf) then
!       PENALISATION SEULE
        do i = 1, ndim
            gt2(i)=cpenfr * vitang(i)
        end do
    else
        do i = 1, ndim
            gt2(i)=lamb1(i)+cstafr * vitang(i)
        end do
    endif
    if (ndim .eq. 3) then
        norme2=sqrt(gt2(1)*gt2(1)+gt2(2)*gt2(2)+gt2(3)*gt2(3))
    else
        norme2=sqrt(gt2(1)*gt2(1)+gt2(2)*gt2(2))
    endif
!
!     ADHER : TRUE SI ADHÉRENCE, FALSE SI GLISSEMENT
    if (norme .le. (1.d0+prec)) then
        adher = .true.
        do j = 1, ndim
            pboul(j)=gt2(j)
        end do
    else
        adher = .false.
        do j = 1, ndim
            pboul(j)=gt2(j)/norme2
        end do
    endif
!
!-----------------------------------------------------------------------
!     CALCUL DE KN(LAMDBA + CSTA [[DX]]/DELTAT )
!
!     ADHERENT
!     OU GLISSANT
!        ET LAMBDA + CSTA [[DX]]/DELTAT EST DANS LA BOULE UNITE
    if (adher .or. ((.not.adher).and.norme2.le.(1.d0+prec))) then
!
        call matini(3, 3, 0.d0, kn)
        do i = 1, ndim
            kn(i,i)=1.d0
        end do
!
!     GLISSANT
!       ET LAMBDA + CSTA [[DX]]/DELTAT N'EST PAS DANS LA BOULE UNITE
    else
!
        do i = 1, ndim
            do j = 1, ndim
                kn(i,j)=-gt2(i)*gt2(j)/(norme2*norme2)
            end do
        end do
!
        do i = 1, ndim
            kn(i,i)= kn(i,i) + 1.d0
        end do
!
        do i = 1, ndim
            do j = 1, ndim
                kn(i,j)= kn(i,j)/norme2
            end do
        end do
!
    endif
!
!-----------------------------------------------------------------------
!
!     CALCUL DE PT.KN.P
    if (ndim .eq. 3) then
        call utbtab('ZERO', ndim, ndim, kn, p,&
                    xab, ptknp)
    else
        do i = 1, ndim
            do j = 1, ndim
                p2(i,j)=p(i,j)
                kn2(i,j)=kn(i,j)
            end do
        end do
        call utbtab('ZERO', ndim, ndim, kn2, p2,&
                    xab2, ptknp2)
        do i = 1, ndim
            do j = 1, ndim
                ptknp(i,j)=ptknp2(i,j)
            end do
        end do
    endif
!
!     CALCUL DE Id-KN
    do i = 1, ndim
        do j = 1, ndim
            ik(i,j)= -1.d0 * kn(i,j)
        end do
    end do
    do i = 1, ndim
        ik(i,i)= 1.d0 + ik(i,i)
    end do
!
!-----------------------------------------------------------------------
end subroutine
