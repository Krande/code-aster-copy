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
subroutine mefpre(ndim, alpha, z, cf, dh, &
                  vit, rho, pstat, dpstat, dvit, &
                  itypg, zg, hg, axg, pm, &
                  xig, afluid, cdg, cfg, vitg, &
                  rhog)
! aslint: disable=W1504
    implicit none
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
    integer(kind=8) :: ndim(14)
    real(kind=8) :: alpha, z(*), cf(*), dh, vit(*), rho(*), pstat(*)
    real(kind=8) :: dpstat(*), dvit(*)
!
    integer(kind=8) :: itypg(*)
    real(kind=8) :: zg(*), hg(*), axg(*), xig(*), afluid, pm
    real(kind=8) :: cdg(*), cfg(*), vitg(*), rhog(*)
!     CALCUL DE LA PRESSION ET DU GRADIENT DE PRESSION STATIONNAIRE
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NDIM   : TABLEAU DES DIMENSIONS
! IN  : ALPHA  : COEFFICIENT DE PROPORTIONALITE DE LA PESENTEUR PAR
!                RAPPORT A LA VALEUR STANDARD (9.81). LA PROJECTION DU
!                VECTEUR V SUIVANT Z VAUT 9.81*ALPHA.
! IN  : Z      : COORDONNEES 'Z'  DES DES POINTS DE DISCRETISATION DANS
!                LE REPERE AXIAL
! IN  : CF     : COEFFICIENT DE TRAINEE VISQUEUSE DU FLUIDE LE LONG DES
!                PAROIS, AUX POINTS DE DISCRETISATION
! IN  : DH     : DIAMETRE HYDRAULIQUE
! IN  : VIT    : VITESSE D ECOULEMENT DU FLUIDE AUX POINTS DE
!                DISCRETISATION
! IN  : RHO    : MASSE VOLUMIQUE DU FLUIDE AUX POINTS DE DISCRETISATION
! OUT : PSTAT  : PROFIL DE PRESSION STATIONNAIRE
! OUT : DPSTAT : PROFIL DE GRADIENT DE PRESSION STATIONNAIRE
! --  : DVIT   : TABLEAU DE TRAVAIL, GRADIENT DE VITESSE D ECOULEMENT DU
!                FLUIDE
!
! IN  : ITYPG  : VECTEUR DES TYPES DE GRILLES
! IN  : ZG     : COORDONNEES 'Z' DES POSITIONS DES GRILLES DANS LE
!                 REPERE AXIAL
! IN  : HG     :  VECTEUR DES HAUTEURS DE GRILLE
! IN  : AXG    : VECTEUR DES SECTIONS SOLIDE DES TYPES DE GRILLES
! IN  : XIG    : VECTEUR DES PERIMETRES MOUILLES DES TYPES DE GRILLES
! IN  : AFLUID: SECTION FLUIDE DE L'ECOULEMENT EN L'ABSENCE DE GRILLES
! IN  : PM     : PERIMETRE MOUILLE DE L'ECOULEMENT EN L'ABSENCE
!                DE GRILLES
! IN  : CDG    : VECTEUR DES COEFF DE TRAINEE DES TYPES DE GRILLES
! IN  : CFG    : VECTEUR DES COEEF DE FROTTEMENT DES TYPES DE GRILLES
! IN  : VITG   : VITESSE D'ECOULEMENT DU  FLUIDE AUX POINTS DE
!                POSITIONNEMENT DES GRILLES
! IN  : RHOG   : MASSE VOLUMIQUE DU FLUIDE AUX MEMES POINTS
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k, n, nbz, nbgtot, ntypg
    real(kind=8) :: ecart, g, pi
    real(kind=8), pointer :: cfnew(:) => null()
    real(kind=8), pointer :: deltap(:) => null()
! ----------------------------------------------------------------------
    call jemarq()
!
! --- LECTURE DES DIMENSIONS
    nbz = ndim(1)
    ntypg = ndim(13)
    nbgtot = ndim(14)
!
! --- CREATION DES OBJETS DE TRAVAIL
    if (ntypg .ne. 0) then
        AS_ALLOCATE(vr=deltap, size=nbgtot)
        AS_ALLOCATE(vr=cfnew, size=nbgtot)
    end if
!
    pi = r8pi()
!
! --- ACCELERATION DE LA PESANTEUR
    g = 9.81d0*alpha
!
! --- VITESSE MOYENNE D ECOULEMENT ET MASSE VOLUMIQUE MOYENNE
!
!
! --- CALCUL DE VIT'(Z) -> G(Z)
! --- MINIMISATION QUADRATIQUE DES RESTES DES
! --- DEVELOPPEMENTS DE TAYLOR DE VIT(Z)
! --- A GAUCHE ET A DROITE
!
    dvit(1) = (vit(2)-vit(1))/(z(2)-z(1))
!
    do n = 2, nbz-1
        dvit(n) = ( &
                  ( &
                  vit(n+1)-vit(n))*(z(n+1)-z(n))+(vit(n-1)-vit(n))*(z(n-1)-z(n)))/((z(n+1)-z&
                  &(n))*(z(n+1)-z(n))+(z(n-1)-z(n))*(z(n-1)-z(n) &
                  ) &
                  )
    end do
!
    dvit(nbz) = (vit(nbz)-vit(nbz-1))/(z(nbz)-z(nbz-1))
!
! --- CALCUL DU PROFIL DE GRADIENT DE PRESSION STATIONNAIRE
!
    do n = 1, nbz
        dpstat(n) = -rho(n)*vit(n)*dvit(n)+rho(n)*g-2.d0*rho(n)*cf(n)*abs(vit(n))*vit(n)/pi/d&
                    &h
    end do
!
! --- CALCUL DU PROFIL DE PRESSION STATIONNAIRE
!
    pstat(1) = 0.d0
    do n = 2, nbz
        pstat(n) = pstat(n-1)+(dpstat(n-1)+dpstat(n))*(z(n)-z(n-1))/2.d0
    end do
!
!--- CALCUL DU SAUT DE PRESSION AU PASSAGE DE CHAQUE GRILLE
!
    if (ntypg .ne. 0) then
!
        do j = 1, nbgtot
            deltap(j) = 0.d0
        end do
!
        do i = 2, nbz
            do j = 1, nbgtot
                ecart = (z(i)-zg(j))*(z(i-1)-zg(j))
!
                if (ecart .le. 0.d0) then
                    cfnew(j) = (cf(i-1)*(z(i)-zg(j))+cf(i)*(zg( &
                                                            j)-z(i-1)))/(z(i)-z(i-1))
                end if
            end do
        end do
!
        do j = 1, nbgtot
            do k = 1, ntypg
                if (itypg(j) .eq. k) then
                    deltap(j) = 0.5d0*rhog(j)*abs(vitg(j))*vitg(j)*(axg(k)*cdg(k)+xig(k)*hg(k)*&
                                &cfg(j))/afluid+0.5d0*rhog(j)*abs(vitg(j))*vitg(j)*(1.d0-(1.d0&
                                &-axg(k)/afluid)**2)*pm*hg(k)*cfnew(j)/afluid
                end if
            end do
        end do
!
        do n = 2, nbz
            do j = 1, nbgtot
                ecart = (z(n)-zg(j))*(z(n-1)-zg(j))
                if (ecart .le. 0.d0) then
                    do k = n, nbz
                        pstat(k) = pstat(k)-deltap(j)
                    end do
                end if
            end do
        end do
!
    end if
!
    if (ntypg .ne. 0) then
        AS_DEALLOCATE(vr=deltap)
        AS_DEALLOCATE(vr=cfnew)
    end if
    call jedema()
end subroutine
