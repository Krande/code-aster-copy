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

subroutine lcmmop(fami, kpg, ksp, rela_comp, nbcomm, &
                  cpmono, nmat, nvi, vini, x, &
                  dtime, mod, coeft, epsd, detot, &
                  coel, nbphas, nfs, nsg, toutms, &
                  dvin, nhsr, numhsr, hsr, itmax, &
                  toler, iret)
! aslint: disable=W1306,W1504
    implicit none
#include "asterfort/calsig.h"
#include "asterfort/lcloca.h"
#include "asterfort/lcmmec.h"
#include "asterfort/lcmmfe.h"
#include "asterfort/lcmmfi.h"
#include "asterfort/lcnrts.h"
    integer(kind=8) :: kpg, ksp, nmat, nbcomm(nmat, 3), nvi, nbphas, itmax, iret
    integer(kind=8) :: nfs, nsg, nhsr, numhsr(*)
    real(kind=8) :: vini(*), dvin(*), x, dtime, coeft(nmat), coel(nmat)
    real(kind=8) :: sigi(6), epsd(6), detot(6)
!     POUR GAGNER EN TEMPS CPU
    real(kind=8) :: toutms(nbphas, nfs, nsg, 7)
    real(kind=8) :: toler, hsr(nsg, nsg, nhsr)
    character(len=*) :: fami
    character(len=16) :: rela_comp

! ======================================================================
! ======================================================================
!       IN   FAMI   : FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!            KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!           COMP    :  NOM DU MODELE DE COMPORTEMENT
!           MOD     :  TYPE DE MODELISATION
!           IMAT    :  ADRESSE DU MATERIAU CODE
!         NBCOMM :  NOMBRE DE COEF MATERIAU PAR FAMILLE
!         CPMONO :  NOMS DES LOIS MATERIAU PAR FAMILLE
!           PGL   : MATRICE DE PASSAGE GLOBAL LOCAL
!           NVI     :  NOMBRE DE VARIABLES INTERNES
!           VINI    :  VARIABLES INTERNES A T
!           X       :  INTERVALE DE TEMPS ADAPTATIF
!           DTIME   :  INTERVALE DE TEMPS
!           COEFT   :  COEFFICIENTS MATERIAU INELASTIQUE A T
!           EPSD    :  DEFORMATION TOTALE A T
!           DETOT   :  INCREMENT DE DEFORMATION TOTALE
!     OUT:
!           DVIN    :  DERIVEES DES VARIABLES INTERNES A T
! INTEGRATION DES LOIS POLYCRISTALLINES PAR UNE METHODE DE RUNGE KUTTA
!
!     CETTE ROUTINE FOURNIT LA DERIVEE DE L ENSEMBLE DES VARIABLES
!     INTERNES DU MODELE
!
!       OBJETS DE STOCKAGE DES COMPORTEMENTS:
!           COEFT(*) = Fractions Volumiques et angles de chaque phase
!                      + COEFFICIENT DE CHAQUE COMPORTEMENT MONOCRSITAL
!                        pour chaque famille de systèmes de glissement
!                        à la température actuelle (COEFTF)
!                       et à la température précédente (COEFTD)
!           NBCOMM = indices des coefficents de chaque comportement
!                    dans COEFT(*,2)
!           CPMONO = noms des différentes "briques" de comportement
!
!      STRUCTURE DES OBJETS CREES
!
!           COEFT(*) : Nombre de monocristaux
!                        indice debut premier monocristal
!                        indice debut deuxième monocristal
!..............
!                        indice debut denier monocristal
!                        indice des paramètes localisation
!                        Fv et 3 angles par phase
!           pour chaque monocristal différent
!                 par famille de système de glissement
!                    nb coef écoulement
!                       numéro de la loi d'écoulement
!                       + coef,
!                    nb coef écrou isot + num_loi + coef,
!                    nb coef ecou cine + num_loi + coef
!                        puis 2 (ou plus) paramètres localisation
!
!
!           CPMONO(*) : nom de la methode de localisation
!                 puis, pour chaque matériau différent
!                 nom du monocristal, nombre de familles SG, et,
!                    par famille de système de glissement
!                       Nom de la famille
!                       Nom du materiau
!                       Nom de la loi d'écoulement
!                       Nom de la loi d'écrouissage isotrope
!                       Nom de la loi d'écrouissage cinématique
!                       Nom de la loi d'élasticité (ELAS ou ELAS_ORTH)
!
!           NBCOMM(*,3) :
!                        Colonne 1      Colonne 2      Colonne3
!                    _____________________________________________
!
!            Ligne 1     Nb phases      Nb var.int.   Nb monocristaux
!                                                     différents
!   pour chaque phase g  Num ligne g    Ind CPMONO    ind frac vol
!   ..................
!   ...................
!   pour chaque phase
!   pour la localisation  indice coef   nb param      0
!   phase g              nb fam g       0            NVIg
!                ... et pour chaque famille de système de glissement :
!             famille 1  ind coef       ind coef      ind coef
!                        ecoulement     ecr iso       ecr cin
!    .....
!         (ind signifie l'indice dans COEFT(*)
!                    _____________________________________________
!     7 variables : tenseur EVP + Norme(EVP)
!    description des variables internes :
!    pour chaque phase
!        6 variables : beta ou epsilonp par phase
!    pour chaque phase
!        pour chaque systeme de glissement
!              3 variables Alpha, Gamma, P
!   1 variable : indic
!     ----------------------------------------------------------------
    character(len=8) :: mod
    character(len=16) :: necoul, necris, necrci
    character(len=24) :: cpmono(5*nmat+1)
    character(len=16) :: loca
    real(kind=8) :: vis(3), dt, evg(6), dl, da, gammas
    real(kind=8) :: evi(6), sigg(6), rp, devg(6), fv
    real(kind=8) :: devi(6), ms(6), taus, dgamma, dalpha, dp
    real(kind=8) :: devgeq, dbeta, beta, dvineq, granb(6)
    real(kind=8) :: crit, sgns, expbp(nsg), dy(nsg)
    integer(kind=8) :: itens, nbfsys, i, nuvi, ifa, nbsys, is, iv
    integer(kind=8) :: indpha, indfv, iphas, indcp, indfa, iexp
    integer(kind=8) :: ifl, nuecou, ihsr
    integer(kind=8) :: irr, decirr, nbsyst, decal, gdef
    common/polycr/irr, decirr, nbsyst, decal, gdef
!     ----------------------------------------------------------------
! --  VARIABLES INTERNES
!
    do itens = 1, 6
        evi(itens) = vini(itens)
        devi(itens) = 0.d0
    end do
    do i = 1, nsg
        dy(i) = 0.d0
    end do
    iret = 0
!
    call calsig(fami, kpg, ksp, evi, mod, &
                rela_comp, vini, x, dtime, epsd, &
                detot, nmat, coel, sigi)
!
! LOCALISATION
!   RECUPERATION DU NOMBRE DE PHASES
!      NBPHAS=NBCOMM(1,1)
    loca = cpmono(1) (1:16)
!     CALCUL DE  B
    do i = 1, 6
        granb(i) = 0.d0
    end do
    do i = 1, 6
        do iphas = 1, nbphas
            indfv = nbcomm(1+iphas, 3)
            fv = coeft(indfv)
            granb(i) = granb(i)+fv*vini(7+6*(iphas-1)+i)
        end do
    end do
!
!     DEBUT DES VARIABLES INTERNES DES SYSTEMES DE GLISSEMENT
    nuvi = 7+6*nbphas
    decal = nuvi
!
    nbsyst = 0
!
    do iphas = 1, nbphas
!        INDPHA indice debut phase IPHAS dans NBCOMM
        indpha = nbcomm(1+iphas, 1)
        indfv = nbcomm(1+iphas, 3)
!
!         recuperer l'orientation de la phase et la proportion
!         INDORI=INDFV+1
        fv = coeft(indfv)
        call lcloca(coeft, nmat, nbcomm, nbphas, sigi, &
                    vini, iphas, granb, loca, sigg)
        nbfsys = nbcomm(indpha, 1)
        indcp = nbcomm(1+iphas, 2)
!        Nombre de variables internes de la phase (=monocristal)
!         NVIG=NBCOMM(INDPHA,3)
        do itens = 1, 6
            devg(itens) = 0.d0
            evg(itens) = 0.d0
        end do
        ihsr = numhsr(iphas)
!
        do ifa = 1, nbfsys
!
            necoul = cpmono(indcp+5*(ifa-1)+3) (1:16)
            necris = cpmono(indcp+5*(ifa-1)+4) (1:16)
            necrci = cpmono(indcp+5*(ifa-1)+5) (1:16)
!
            nbsys = nint(toutms(iphas, ifa, 1, 7))
!
!           indice de la famille IFA
            indfa = indpha+ifa
!
            ifl = nbcomm(indfa, 1)
            nuecou = nint(coeft(ifl))
!
            do is = 1, nbsys
!              VARIABLES INTERNES DU SYST GLIS
                do iv = 1, 3
                    nuvi = nuvi+1
                    vis(iv) = vini(nuvi)
                end do
                dvin(nuvi-2) = 0.d0
                dvin(nuvi-1) = 0.d0
                dvin(nuvi) = 0.d0
!
!              CALCUL DE LA SCISSION REDUITE
!              PROJECTION DE SIG SUR LE SYSTEME DE GLISSEMENT
!              TAU      : SCISSION REDUITE TAU=SIG:MS
                do i = 1, 6
                    ms(i) = toutms(iphas, ifa, is, i)
                end do
!
                taus = 0.d0
                do i = 1, 6
                    taus = taus+sigg(i)*ms(i)
                end do
!
!              ECROUISSAGE ISOTROPE
!
!              DECAL est le début des systemes de glissement de la
!              phase en cours
!              NVIG est le nombre de variables internes dela phase G
!
!               IF (NECOUL.NE.'KOCKS_RAUCH') THEN
                if (nuecou .ne. 4) then
!
                    iexp = 0
                    if (is .eq. 1) iexp = 1
                    call lcmmfi(coeft, indfa, nmat, nbcomm, necris, &
                                is, nbsys, vini, decal, dy(1), &
                                nfs, nsg, hsr(1, 1, ihsr), iexp, expbp, &
                                rp)
!
                end if
!
!              ECOULEMENT VISCOPLASTIQUE:
!              ROUTINE COMMUNE A L'IMPLICITE (PLASTI-LCPLNL)
!              ET L'EXPLICITE (NMVPRK-GERPAS-RK21CO-RDIF01)
!              CAS IMPLCITE : IL FAUT PRENDRE EN COMPTE DTIME
!              CAS EXPLICITE : IL NE LE FAUT PAS (VITESSES)
!              D'OU :
                dt = 1.d0
!
                call lcmmfe(taus, coeft, coel, indfa, nmat, &
                            nbcomm, necoul, is, nbsys, vini, &
                            dy(1), rp, vis(1), vis(2), dt, &
                            dalpha, dgamma, dp, crit, sgns, &
                            nfs, nsg, hsr(1, 1, ihsr), iret)
                if (dp .gt. 0.d0) then
!
!                 ECROUISSAGE CINEMATIQUE
!
                    if (nuecou .lt. 4) then
                        call lcmmec(coeft, indfa, nmat, nbcomm, necrci, &
                                    itmax, toler, vis(1), dgamma, dalpha, &
                                    iret)
                        if (iret .ne. 0) goto 999
                    end if
!                 DEVG designe ici DEPSVPG
                    do itens = 1, 6
                        devg(itens) = devg(itens)+ms(itens)*dgamma
                    end do
!
!                 EVG designe ici EPSVPG
                    if (loca .eq. 'BETA') then
                        gammas = vis(2)+dgamma
                        do itens = 1, 6
                            evg(itens) = evg(itens)+ms(itens)*gammas
                        end do
                    end if
!
                    dvin(nuvi-2) = dalpha
                    dvin(nuvi-1) = dgamma
                    dvin(nuvi) = dp
                else
                    dvin(nuvi-2) = 0.d0
                    dvin(nuvi-1) = 0.d0
                    dvin(nuvi) = 0.d0
                end if
                nbsyst = nbsyst+1
            end do
!
        end do
!
        decal = nuvi
!
!         "homogenesisation" des déformations viscoplastiques
        do i = 1, 6
            devi(i) = devi(i)+fv*devg(i)
        end do
        devgeq = lcnrts(devg)/1.5d0
!         localisation BETA
        if (loca .eq. 'BETA') then
            dl = coeft(nbcomm((nbphas+2), 1)+1)
            da = coeft(nbcomm((nbphas+2), 1)+2)
            do i = 1, 6
                beta = vini(7+6*(iphas-1)+i)
                dbeta = devg(i)-dl*(beta-da*evg(i))*devgeq
                dvin(7+6*(iphas-1)+i) = dbeta
            end do
        else
            do i = 1, 6
                dvin(7+6*(iphas-1)+i) = devg(i)
            end do
!
        end if
!
! fin boucle sur nombre de phases
    end do
!
! --    DERIVEES DES VARIABLES INTERNES
!
    do itens = 1, 6
        dvin(itens) = devi(itens)
    end do
!     Norme de DEVP cumulée
    dvineq = lcnrts(devi)/1.5d0
!
    dvin(7) = dvineq
    do itens = 1, 6*nbphas
        dvin(nuvi+i) = 0.d0
    end do
    dvin(nvi) = 0
999 continue
end subroutine
