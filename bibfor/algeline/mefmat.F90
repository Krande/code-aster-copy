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
subroutine mefmat(ndim, numgrp, nbz, nbgrp, nbmod, &
                  matma, dcent, cp, cf, vit, &
                  rho, pstat, dpstat, rint, phix, &
                  phiy, z, matm, matr, mata, &
                  itypg, axg, zg, rhog, vitg, &
                  cdg, cpg)
! aslint: disable=W1504
    implicit none
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mefin1.h"
#include "asterfort/mefin2.h"
#include "asterfort/mefin3.h"
#include "asterfort/mefin4.h"
#include "asterfort/mefin5.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: ndim(14), numgrp(*), nbmod, nbz
    real(kind=8) :: matma(*), dcent(*)
    real(kind=8) :: cp(*), cf(*), vit(0:*), rho(0:*), pstat(*), dpstat(*)
    integer(kind=8) :: nbgrp
    real(kind=8) :: matm(nbmod, nbmod), rint(*), phix(nbz*nbgrp, nbmod)
    real(kind=8) :: phiy(nbz*nbgrp, nbmod), z(*)
    real(kind=8) :: matr(nbmod, nbmod), mata(nbmod, nbmod)
!
    integer(kind=8) :: itypg(*)
    real(kind=8) :: zg(*), axg(*), rhog(*), vitg(*), cdg(*), cpg(*)
!     CALCUL DES MATRICES DE MASSE, DE RAIDEUR, D AMORTISSEMENT SOUS
!     ECOULEMENT
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NDIM   : TABLEAU DES DIMENSIONS
! IN  : NUMGRP : INDICES DES GROUPES D EQUIVALENCE
! IN  : NBZ     : NOMBRE DE POINTS DE DISCRETISATION
! IN  : NBGRP  : NOMBRE DE GROUPES D'EQUIVALENCE DE CYLINDRES
! IN  : NBMOD  : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
! IN  : MATMA  : VECTEUR CONTENANT LES MATRICES MODALES, MASSE,RIGIDITE,
!                AMORTISSEMENT
! IN  : DCENT  : VECTEUR CONTENANT LES TABLEAUX DE COEFFICIENTS ET
!                LES MATRICES EN AIR
! IN  : CP     : COEFFICIENT DE PORTANCE CP DU FLUIDE AUTOUR D UN
!                CYLINDRE INCLINE, AUX POINTS DE DISCRETISATION
! IN  : CF     : COEFFICIENT DE TRAINEE VISQUEUSE DU FLUIDE LE LONG DES
!                PAROIS, AUX POINTS DE DISCRETISATION
! IN  : VIT    : VITESSE D ECOULEMENT DU FLUIDE AUX POINTS DE
!                DISCRETISATION
! IN  : RHO    : MASSE VOLUMIQUE DU FLUIDE AUX POINTS DE DISCRETISATION
! IN  : PSTAT  : PROFIL DE PRESSION STATIONNAIRE
! IN  : DPSTAT : PROFIL DE GRADIENT DE PRESSION STATIONNAIRE
! IN  : RINT   : RAYONS DES CYLINDRES
! IN  : PHIX   : DEFORMEES MODALES INTERPOLEES DANS LE REPERE AXIAL
! IN  : PHIY   : DEFORMEES MODALES INTERPOLEES DANS LE REPERE AXIAL
! OUT : MATM   : MATRICE DE MASSE AJOUTEE REPRESENTANT LA PROJECTION DES
!                EFFORTS FLUIDES INERTIELS DANS LA BASE DES DEFORMEES
!                MODALES DES CYLINDRES
! OUT : MATR   : MATRICE DE RAIDEUR  AJOUTEE REPRESENTANT LA PROJECTION
!                DES EFFORTS FLUIDES DE RAIDEUR DANS LA BASE DES
!                DEFORMEES MODALES DES CYLINDRES
! OUT : MATA   : MATRICE D AMORTISSEMENT AJOUTEE REPRESENTANT LA
!                PROJECTION DES EFFORTS FLUIDES D AMORTISSEMENT DANS LA
!                BASE DES DEFORMEES MODALES DES CYLINDRES
!
! IN  : ITYPG  : VECTEUR DES GROUPES D'APPARTENANCE DES GRILLES
! IN  : ZG     : POINTS DE DISCRETISATION DES GRILLES (EN LEURMILIEU)
! IN  : AXG    : SECTION SOLIDE DES TYPES DE GRILLES
! IN  : RHOG   : PROFIL DE MASSE VOLUMIQUE DE L'ECOULEMENT AU NIVEAU
!                         DES GRILLES
! IN  : VITG   : PROFIL DE VITESSE DE L'ECOULEMENT AU NIVEAU DES GRILLES
! IN  : CDG    : COEFF DE TRAINEE POUR CHAQUE TYPE DE GRILLES
! IN  : CPG    : PENTE DU COEFF DE PORTANCE POUR CHAQUE TYPE DE GRILLES
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j
    integer(kind=8) :: imod, igrp, jmod, jgrp
    integer(kind=8) :: ncyl
    real(kind=8) :: rayo
!
    integer(kind=8) :: kg, k, ngz1, ngz2
    real(kind=8) :: ecart
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: idphxg, idphyg, ig, ih, imataa, imatra
    integer(kind=8) :: iphixg, iphiyg, ippxx, ippxy, ippyx, ippyy, ivnxx
    integer(kind=8) :: ivnxy, ivnyx, ivnyy, nbcyl, nbgtot
    integer(kind=8) :: ntypg
    real(kind=8) :: aire, pi, rho0, vit0
    real(kind=8), pointer :: aireg(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
! --- LECTURE DES DIMENSIONS
    nbcyl = ndim(3)
    ntypg = ndim(13)
    nbgtot = ndim(14)
!
! --- CREATION DES OBJETS DE TRAVAIL
    call wkvect('&&MEFMAT.TMP.GH', 'V V R', nbz*2, ig)
    ih = ig+nbz
    if (ntypg .ne. 0) then
        AS_ALLOCATE(vr=aireg, size=ntypg)
        call wkvect('&&MEFMAT.PHIXG', 'V V R', nbgrp*nbgtot*nbmod, iphixg)
        call wkvect('&&MEFMAT.PHIYG', 'V V R', nbgrp*nbgtot*nbmod, iphiyg)
        call wkvect('&&MEFMAT.DPHIXG', 'V V R', nbgrp*nbgtot*nbmod, idphxg)
        call wkvect('&&MEFMAT.DPHIYG', 'V V R', nbgrp*nbgtot*nbmod, idphyg)
    end if
!
    pi = r8pi()
!
! --- VITESSE MOYENNE D ECOULEMENT ET MASSE VOLUMIQUE MOYENNE
!
    rho0 = rho(0)
    vit0 = vit(0)
!
! --- DECALAGES DES TABLEAUX DE COEFFICIENTS ET DES MATRICES EN AIR
! --- DANS LE VECTEUR DCENT
!
    ippxx = nbcyl+nbcyl+nbcyl*nbcyl+nbcyl*nbcyl
    ippxy = ippxx+nbcyl*nbgrp
    ippyx = ippxy+nbcyl*nbgrp
    ippyy = ippyx+nbcyl*nbgrp
    ivnxx = ippyy+nbcyl*nbgrp
    ivnxy = ivnxx+nbcyl*nbgrp
    ivnyx = ivnxy+nbcyl*nbgrp
    ivnyy = ivnyx+nbcyl*nbgrp
!
! --- DECALAGES DES MATRICES MODALES DANS LE VECTEUR MATMA
!
    imatra = nbmod
    imataa = imatra+nbmod
!
! --- INITIALISATIONS
!
    matm(:, :) = 0.d0
    mata(:, :) = 0.d0
    matr(:, :) = 0.d0
!
!---  SECTION SOLIDE D'UNE CELLULE ELEMENTAIRE DE GRILLE
!
    if (ntypg .ne. 0) then
!
        do k = 1, ntypg
            aireg(k) = axg(k)/dble(nbcyl)
        end do
!
!---  CALCUL (PAR INTERPOLATION  LINEAIRE) DES DEFORMEES MODALES
!---  ET DE LEUR GRADIENT AUX COTES ZG DES GRILLES
!
        do imod = 1, nbmod
            do igrp = 1, nbgrp
                do i = 2, nbz
                    do j = 1, nbgtot
                        ecart = (z(i)-zg(j))*(z(i-1)-zg(j))
                        if (ecart .le. 0.d0) then
                            zr(iphixg+(imod-1)*nbgrp*nbgtot+(igrp-1)* &
                               nbgtot+j-1) = phix((igrp-1)*nbz+i-1, imod)* &
                                          (z(i)-zg(j))/(z(i)-z(i-1))+ &
                                          phix((igrp-1)*nbz+i, imod)* &
                                          (zg(j)-z(i-1))/(z(i)-z(i-1))
!
                            zr(iphiyg+(imod-1)*nbgrp*nbgtot+(igrp-1)* &
                               nbgtot+j-1) = phiy((igrp-1)*nbz+i-1, imod)* &
                                          (z(i)-zg(j))/(z(i)-z(i-1))+ &
                                          phiy((igrp-1)*nbz+i, imod)* &
                                          (zg(j)-z(i-1))/(z(i)-z(i-1))
!
                            zr(idphxg+(imod-1)*nbgrp*nbgtot+(igrp-1)* &
                               nbgtot+j-1) = (phix((igrp-1)*nbz+i, imod)- &
                                              phix((igrp-1)*nbz+i-1, imod))/ &
                                          (z(i)-z(i-1))
!
                            zr(idphyg+(imod-1)*nbgrp*nbgtot+(igrp-1)* &
                               nbgtot+j-1) = (phiy((igrp-1)*nbz+i, imod)- &
                                              phiy((igrp-1)*nbz+i-1, imod))/ &
                                          (z(i)-z(i-1))
!
                        end if
                    end do
                end do
            end do
        end do
!
    end if
!
    do jmod = 1, nbmod
        do imod = 1, nbmod
!
            do jgrp = 1, nbgrp
                do igrp = 1, nbgrp
!
                    do i = 1, nbcyl
                        if (numgrp(i) .eq. igrp) then
                            rayo = rint(i)
                        end if
                    end do
                    aire = pi*rayo*rayo
!
! --- CONTRIBUTION DES EFFORTS NORMAUX DE FROTTEMENT VISQUEUX
! --- -> TERMES D'AMORTISSEMENT ET DE RAIDEUR AJOUTES
!
! --- AMORTISSEMENT AJOUTE
!
                    ncyl = 0
                    if (igrp .eq. jgrp) then
                        do i = 1, nbcyl
                            if (numgrp(i) .eq. igrp) ncyl = ncyl-1
                        end do
                    end if
!
                    mata(imod, jmod) = mata(imod, jmod)- &
                                       rho0*abs(vit0)*rayo* &
                                       ((dcent(ivnxx+nbcyl*(jgrp-1)+igrp)+dble(ncyl))* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phix, cf)+ &
                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phiy, cf)+ &
                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phix, cf)+ &
                                        (dcent(ivnyy+nbcyl*(jgrp-1)+igrp)+dble(ncyl))* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phiy, cf))
!
!
                    mata(imod, jmod) = mata(imod, jmod)- &
                                       rho0*abs(vit0)*rayo* &
                                       ((dcent(ivnxx+nbcyl*(jgrp-1)+igrp)+dble(ncyl))* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phix, cp)+ &
                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phiy, cp)+ &
                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phix, cp)+ &
                                        (dcent(ivnyy+nbcyl*(jgrp-1)+igrp)+dble(ncyl))* &
                                        mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phiy, cp))
!
!
! ---  RAIDEUR AJOUTEE
!
                    matr(imod, jmod) = matr(imod, jmod)- &
                                       rho0*abs(vit0)*rayo* &
                                       (dcent(ivnxx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phix, vit, cf, zr(ig))+ &
                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phiy, vit, cf, zr(ig))+ &
                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phix, vit, cf, zr(ig))+ &
                                        dcent(ivnyy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phiy, vit, cf, zr(ig)))
!
!
                    matr(imod, jmod) = matr(imod, jmod)- &
                                       rho0*abs(vit0)*rayo* &
                                       ((dcent(ivnxx+nbcyl*(jgrp-1)+igrp)+dble(ncyl))* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phix, vit, cp, zr(ig))+ &
                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phiy, vit, cp, zr(ig))+ &
                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phix, vit, cp, zr(ig))+ &
                                        (dcent(ivnyy+nbcyl*(jgrp-1)+igrp)+dble(ncyl))* &
                                        mefin4(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phiy, vit, cp, zr(ig)))
!
!
!
! --- CONTRIBUTION DES EFFORTS DE PRESSION PERTURBEE
! --- -> TERMES DE MASSE, AMORTISSEMENT ET RAIDEUR AJOUTES
!
! --- MASSE AJOUTEE
!
                    matm(imod, jmod) = matm(imod, jmod)- &
                                       aire*(dcent(ippxx+nbcyl*(jgrp-1)+igrp)* &
                                             mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                                    jgrp, z, phix, phix, rho)+ &
                                             dcent(ippxy+nbcyl*(jgrp-1)+igrp)* &
                                             mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                                    jgrp, z, phix, phiy, rho)+ &
                                             dcent(ippyx+nbcyl*(jgrp-1)+igrp)* &
                                             mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                                    jgrp, z, phiy, phix, rho)+ &
                                             dcent(ippyy+nbcyl*(jgrp-1)+igrp)* &
                                             mefin1(nbz, nbgrp, imod, igrp, jmod, &
                                                    jgrp, z, phiy, phiy, rho))
!
! --- AMORTISSEMENT AJOUTE
!
                    mata(imod, jmod) = mata(imod, jmod)- &
                                       2.d0*rho0*vit0*aire* &
                                       (dcent(ippxx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin2(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phix, zr(ig))+ &
                                        dcent(ippxy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin2(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phiy, zr(ig))+ &
                                        dcent(ippyx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin2(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phix, zr(ig))+ &
                                        dcent(ippyy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin2(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phiy, zr(ig)))
!
! --- RAIDEUR AJOUTEE
!
                    matr(imod, jmod) = matr(imod, jmod)- &
                                       rho0*vit0*aire* &
                                       (dcent(ippxx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin3(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phix, vit, zr(ig), zr(ih))+ &
                                        dcent(ippxy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin3(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phix, phiy, vit, zr(ig), zr(ih))+ &
                                        dcent(ippyx+nbcyl*(jgrp-1)+igrp)* &
                                        mefin3(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phix, vit, zr(ig), zr(ih))+ &
                                        dcent(ippyy+nbcyl*(jgrp-1)+igrp)* &
                                        mefin3(nbz, nbgrp, imod, igrp, jmod, &
                                               jgrp, z, phiy, phiy, vit, zr(ig), zr(ih)))
!
                end do
            end do
!
            do igrp = 1, nbgrp
!
                ncyl = 0
                do i = 1, nbcyl
                    if (numgrp(i) .eq. igrp) ncyl = ncyl+1
                end do
!
                do i = 1, nbcyl
                    if (numgrp(i) .eq. igrp) then
                        rayo = rint(i)
                    end if
                end do
                aire = pi*rayo*rayo
!C
! ---    CONTRIBUTION DES EFFORTS DE PRESSION STATIONNAIRE
! ---    -> TERMES DE RAIDEUR AJOUTEE
!
! ---    RAIDEUR AJOUTEE
!
                matr(imod, jmod) = matr(imod, jmod)- &
                                   aire*ncyl* &
                                   (mefin3(nbz, nbgrp, imod, igrp, jmod, igrp, &
                                           z, phix, phix, pstat, zr(ig), zr(ih))+ &
                                    mefin3(nbz, nbgrp, imod, igrp, jmod, igrp, z, &
                                           phiy, phiy, pstat, zr(ig), zr(ih))+ &
                                    mefin5(nbz, nbgrp, imod, igrp, jmod, igrp, &
                                           z, phix, phix, dpstat, zr(ig))+ &
                                    mefin5(nbz, nbgrp, imod, igrp, jmod, &
                                           igrp, z, phiy, phiy, dpstat, zr(ig)))
!
            end do
!
!---    CONTRIBUTION DES EFFORTS DE CONTRAINTES SUR LES GRILLES
!           ---> TERMES D'AMORTISSEMENT ET DE RAIDEUR AJOUTES
!
            if (ntypg .ne. 0) then
!
                do kg = 1, nbgtot
!
                    do jgrp = 1, nbgrp
                        do igrp = 1, nbgrp
!
                            ncyl = 0
                            if (igrp .eq. jgrp) then
                                do i = 1, nbcyl
                                    if (numgrp(i) .eq. igrp) ncyl = ncyl-1
                                end do
                            end if
!
                            ngz1 = (igrp-1)*nbgtot+kg
                            ngz2 = (jgrp-1)*nbgtot+kg
!
!---   AMORTISSEMENT AJOUTE
!
                            do k = 1, ntypg
                                if (itypg(kg) .eq. k) then
                                    mata(imod, jmod) = mata(imod, jmod)- &
                                                       0.5d0*rhog(kg)* &
                                                       abs(vitg(kg))*aireg(k)*cpg(k)* &
                                                       ((dcent(ivnxx+nbcyl*(jgrp-1)+igrp)+ &
                                                         dble(ncyl))* &
                                                        zr(iphixg+(imod-1)*nbgtot*nbgrp+ &
                                                           ngz1-1) &
                                                        *zr(iphixg+(jmod-1)*nbgtot*nbgrp+ &
                                                            ngz2-1)+ &
                                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphixg+(imod-1)*nbgtot*nbgrp+ &
                                                           ngz1-1)* &
                                                        zr(iphiyg+(jmod-1)*nbgtot*nbgrp+ &
                                                           ngz2-1)+ &
                                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphiyg+(imod-1)*nbgtot*nbgrp+ &
                                                           ngz1-1)* &
                                                        zr(iphixg+(jmod-1)*nbgtot*nbgrp+ &
                                                           ngz2-1)+ &
                                                        (dcent(ivnyy+nbcyl*(jgrp-1)+igrp)+ &
                                                         dble(ncyl))* &
                                                        zr(iphiyg+(imod-1)*nbgtot*nbgrp+ &
                                                           ngz1-1)* &
                                                        zr(iphiyg+(jmod-1)*nbgtot*nbgrp+ &
                                                           ngz2-1))
!
                                    mata(imod, jmod) = mata(imod, jmod)- &
                                                       0.5d0*rhog(kg)* &
                                                       abs(vitg(kg))*aireg(k)*cdg(k)* &
                                                       ((dcent(ivnxx+nbcyl*(jgrp-1)+igrp)+ &
                                                         dble(ncyl))* &
                                                        zr(iphixg+(imod-1)*nbgtot*nbgrp+ngz1-1)* &
                                                        zr(iphixg+(jmod-1)*nbgtot*nbgrp+ngz2-1)+ &
                                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphixg+(imod-1)*nbgtot*nbgrp+ngz1-1)* &
                                                        zr(iphiyg+(jmod-1)*nbgtot*nbgrp+ngz2-1)+ &
                                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphiyg+(imod-1)*nbgtot*nbgrp+ngz1-1)* &
                                                        zr(iphixg+(jmod-1)*nbgtot*nbgrp+ngz2-1)+ &
                                                        (dcent(ivnyy+nbcyl*(jgrp-1)+igrp)+ &
                                                         dble(ncyl))* &
                                                        zr(iphiyg+(imod-1)*nbgtot*nbgrp+ngz1-1)* &
                                                        zr(iphiyg+(jmod-1)*nbgtot*nbgrp+ngz2-1))
                                end if
                            end do
!
! ---  RAIDEUR AJOUTEE
!
                            do k = 1, ntypg
                                if (itypg(kg) .eq. k) then
                                    matr(imod, jmod) = matr(imod, jmod)- &
                                                       0.5d0*rhog(kg)*abs(vitg(kg))* &
                                                       aireg(k)*cdg(k)*vitg(kg)* &
                                                       (dcent(ivnxx+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphixg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphxg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1)+ &
                                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphixg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphyg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1)+ &
                                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphiyg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphxg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1)+ &
                                                        dcent(ivnyy+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphiyg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphyg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1))
!
                                    matr(imod, jmod) = matr(imod, jmod)- &
                                                       0.5d0*rhog(kg)*abs(vitg(kg))* &
                                                       aireg(k)*cpg(k)*vitg(kg)* &
                                                       ((dcent(ivnxx+nbcyl*(jgrp-1)+igrp)+ &
                                                         dble(ncyl))* &
                                                        zr(iphixg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphxg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1)+ &
                                                        dcent(ivnxy+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphixg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphyg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1)+ &
                                                        dcent(ivnyx+nbcyl*(jgrp-1)+igrp)* &
                                                        zr(iphiyg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphxg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1)+ &
                                                        (dcent(ivnyy+nbcyl*(jgrp-1)+igrp)+ &
                                                         dble(ncyl))* &
                                                        zr(iphiyg+(imod-1)*nbgrp*nbgtot+ &
                                                           ngz1-1)* &
                                                        zr(idphyg+(jmod-1)*nbgrp*nbgtot+ &
                                                           ngz2-1))
!
                                end if
                            end do
!
                        end do
                    end do
                end do
!
            end if
!
! ---    TERMES DE MASSE, AMORTISSEMENT ET RAIDEUR DE STRUCTURE
!
            if (imod .eq. jmod) then
                matm(imod, jmod) = matm(imod, jmod)+matma(imod)
                mata(imod, jmod) = mata(imod, jmod)+matma(imataa+imod)
                matr(imod, jmod) = matr(imod, jmod)+matma(imatra+imod)
            end if
!
        end do
    end do
!
! --- MENAGE
!
    call jedetr('&&MEFMAT.TMP.GH')
    AS_DEALLOCATE(vr=aireg)
    call jedetr('&&MEFMAT.PHIXG')
    call jedetr('&&MEFMAT.PHIYG')
    call jedetr('&&MEFMAT.DPHIXG')
    call jedetr('&&MEFMAT.DPHIYG')
!
    call jedema()
end subroutine
