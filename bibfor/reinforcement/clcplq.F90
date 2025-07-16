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

subroutine clcplq(typcmb, typco, nb, precs, &
                  ferrsyme, slsyme, ferrcomp, epucisa, &
                  ferrmin, rholmin, rhotmin, compress, cequi, &
                  enrobi, enrobs, sigs, sigci, sigcs, &
                  alphacc, gammas, gammac, facier, eys, typdiag, &
                  fbeton, clacier, uc, um, &
                  wmaxi, wmaxs, sigelsqp, kt, phixi, phixs, phiyi, phiys, &
                  ht, effrts, dnsits, ierrl, ierrt)

!_____________________________________________________________________
!
!      CLCPLQ
!
!      CALCUL DES ARMATURES EN ACIER DANS LES ELEMENTS DE PLAQUE - CAPRA MAURY
!
!      I TYPCMB        TYPE DE COMBINAISON (0 = ELU, 1 = ELS)
!      I TYPCO         CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I NB            NOMBRE DE FACETTES
!      I PRECS         PRECISION ITERATION
!      I FERRSYME      FERRAILLAGE SYMETRIQUE?
!                            (0 = NON, 1 = OUI)
!      I SLSYME        SECTION SEUIL DE TOLERANCE POUR UN FERRAILLAGE SYMETRIQUE
!      I FERRCOMP      FERRAILLAGE DE COMPRESSION ADMIS?
!                      (0 = NON, 1 = OUI)
!      I EPUCISA       IMPACT DE L'EFFORT TRANCHANT ET DE LA TORSION SUR LE
!                      FERRAILLAGE LONGITUDINAL?
!                      (0 = NON, 1 = OUI)
!      I FERRMIN       PRISE EN COMPTE DU FERRA MINI (0 = NON, 1 = OUI, 2 = CODE)
!      I RHOLMIN       RATIO DE FERRAILLAGE LONGI MINI (A RENSEIGNER SI FERMIN='OUI')
!      I RHOTMIN       RATIO DE FERRAILLAGE TRNSV MINI (A RENSEIGNER SI FERMIN='OUI')
!      I COMPRESS      VALORISATION DE LA COMPRESSION POUR LES ACIERS TRANSVERSAUX
!                      (0 = NON, 1 = OUI)
!      I CEQUI         COEFFICIENT D'EQUIVALENCE ACIER/BETON
!      I ENROBI        ENROBAGE DES ARMATURES INFERIEURES
!      I ENROBS        ENROBAGE DES ARMATURES SUPERIEURES
!      I SIGS          CONTRAINTE ADMISSIBLE DANS L'ACIER
!      I SIGCI         CONTRAINTE DE COMPRESSION MAXI DU BETON EN FACE INF
!      I SIGCS         CONTRAINTE DE COMPRESSION MAXI DU BETON EN FACE SUP
!      I ALPHACC       COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DU BETON EN COMPRESSION
!      I GAMMAS        COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DES ACIERS
!      I GAMMAC        COEFFICIENT DE SECURITE SUR LA RESISTANCE
!                      DE CALCUL DU BETON
!      I FACIER        LIMITE D'ELASTICITE DES ACIERS (CONTRAINTE)
!      I EYS           MODULE D'YOUNG DE L'ACIER
!      I TYPDIAG       TYPE DE DIAGRAMME UTILISÉ POUR L'ACIER
!                            TYPDIAG = 1 ("B1" ==> PALIER INCLINÉ)
!                            TYPDIAG = 2 ("B2" ==> PALIER HORIZONTAL)
!      I FBETON        RESISTANCE EN COMPRESSION DU BETON (CONTRAINTE)
!      I CLACIER       CLASSE DE DUCTILITE DES ACIERS (UTILISE POUR EC2) :
!                            CLACIER = 0 ACIER PEU DUCTILE (CLASSE A)
!                            CLACIER = 1 ACIER MOYENNEMENT DUCTILE (CLASSE B)
!                            CLACIER = 2 ACIER FORTEMENT DUCTILE (CLASSE C)
!      I UC            UNITE DES CONTRAINTES :
!                            UC = 0 CONTRAINTES EN Pa
!                            UC = 1 CONTRAINTES EN MPa
!      I UM            UNITE DES DIMENSIONS :
!                            UM = 0 DIMENSIONS EN m
!                            UM = 1 DIMENSIONS EN mm
!      I WMAXI         OUVERTURE MAXIMALE DES FISSURES EN FACE INFÉRIEURE
!      I WMAXS         OUVERTURE MAXIMALE DES FISSURES EN FACE SUPÉRIEURE
!      I SIGELSQP      CONTRAINTE ADMISSIBLE DANS LE BETON À L'ELS QP
!      I KT            COEFFICIENT DE DURÉE DE CHARGEMENT
!      I PHIXI         DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT X
!      I PHIXS         DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT X
!      I PHIYI         DIAMÈTRE APPROXIMATIF DES ARMATURES INFÉRIEURES SUIVANT Y
!      I PHIYS         DIAMÈTRE APPROXIMATIF DES ARMATURES SUPÉRIEURES SUIVANT Y
!      I HT            EPAISSEUR DE LA COQUE
!      I EFFRTS        (DIM 8) TORSEUR DES EFFORTS, MOMENTS, ...
!
!      O DNSITS        (DIM 6) DENSITES
!                            1..4 : SURFACES D'ACIER LONGITUDINAL
!                            5..6 : TRANSVERSAL
!      O IERRL         CODE RETOUR LONGI (0 = OK)
!      O IERRT         CODE RETOUR TRNSV (0 = OK)
!
!_____________________________________________________________________
!
!
    implicit none
!
!
#include "asterfort/cafelu.h"
#include "asterfort/cafels.h"
#include "asterfort/cafelsqp.h"
#include "asterfort/clcopt.h"
#include "asterfort/cftelu.h"
#include "asterfort/cftels.h"
#include "asterfort/trgfct.h"
!
!
    integer(kind=8) :: typcmb
    integer(kind=8) :: typco
    integer(kind=8) :: nb
    integer(kind=8) :: precs
    integer(kind=8) :: ferrsyme
    real(kind=8) :: slsyme
    integer(kind=8) :: ferrcomp
    integer(kind=8) :: epucisa
    integer(kind=8) :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer(kind=8) :: compress
    real(kind=8) :: cequi
    real(kind=8) :: enrobi
    real(kind=8) :: enrobs
    real(kind=8) :: sigs
    real(kind=8) :: sigci
    real(kind=8) :: sigcs
    real(kind=8) :: alphacc
    real(kind=8) :: gammas
    real(kind=8) :: gammac
    real(kind=8) :: facier
    real(kind=8) :: eys
    integer(kind=8) :: typdiag
    real(kind=8) :: fbeton
    integer(kind=8) :: clacier
    integer(kind=8) :: uc
    integer(kind=8) :: um
    real(kind=8) :: wmaxi
    real(kind=8) :: wmaxs
    real(kind=8) :: sigelsqp
    real(kind=8) :: kt
    real(kind=8) :: phixi
    real(kind=8) :: phixs
    real(kind=8) :: phiyi
    real(kind=8) :: phiys
    real(kind=8) :: ht
    real(kind=8) :: effrts(8)
    real(kind=8) :: dnsits(6)
    integer(kind=8) :: ierrl
    integer(kind=8) :: ierrt
!
!       NOMBRE DE DIVISIONS ENTRE -PI/2 ET +PI/2
    real(kind=8) :: fcttab(nb, 6)
!       NOMBRE DE FACETTES COMPRIMEES EN PIVOT C (ELU ET ELS)
    real(kind=8) :: nb_fac_comp_elu, nb_fac_comp_els
!       EFFORT NORMAL DANS CETTE DIRECTION
    real(kind=8) :: effn
!       EFFORT TRANCHANT DANS CETTE DIRECTION
    real(kind=8) :: efft
!       MOMENT DE FLEXION DANS CETTE DIRECTION
    real(kind=8) :: effm
!       DIAMETRE DES BARRES SUPÉRIEURES DANS CETTE DIRECTION
    real(kind=8) :: phisup
!       DIAMETRE DES BARRES INFÉRIEURES DANS CETTE DIRECTION
    real(kind=8) :: phiinf
!       DENSITE DE FERRAILLAGE TRANSVERSAL
    real(kind=8) :: dnstra(nb)
!       SECTIONS DES ACIERS INFERIEURS SUIVANT LES NB FACETTES
    real(kind=8) :: ai(nb)
!       SECTIONS DES ACIERS SUPERIEURS SUIVANT LES NB FACETTES
    real(kind=8) :: as(nb)
!       CONTRAINTE DANS L'ACIER DE TRACTION A L'ELS QP
    real(kind=8) :: Sacier
!       VARIABLE D'ITERATION
    integer(kind=8) :: i
!       AUTRES VARIABLES
    real(kind=8) :: fctm, unite_m, unite_pa, d, ak, uk, alpha, thetab
    real(kind=8) :: Asl, ecinf, ecsup, kvarf, sigmci, sigmcs, sigmsi, sigmss
    real(kind=8) :: wfini, wfins
    integer(kind=8) :: etat, pivot, ierr, ierrlG, ierrtG
!
!   INITIALISATION DES VARIABLES
!
    ierrlG = 0
    ierrtG = 0
    nb_fac_comp_elu = 0
    nb_fac_comp_els = 0
!
!   INITIALISATION DES FACETTES
!
    call trgfct(nb, fcttab)
    do 5 i = 1, 6
        dnsits(i) = -1.d0
5       continue
!
!   BOUCLE SUR LES FACETTES DE CAPRA ET MAURY
!   DETERMINATION DU FERRAILLAGE POUR CHACUNE DES FACETTES

        do 10 i = 1, nb

            ierrl = 0
            ierrt = 0

            effn = fcttab(i, 1)*effrts(1)+fcttab(i, 2)*effrts(2)+fcttab(i, 3)*effrts(3)
            effm = fcttab(i, 1)*effrts(4)+fcttab(i, 2)*effrts(5)+fcttab(i, 3)*effrts(6)
            effn = -effn
            effm = -effm
            efft = abs(effrts(7)*fcttab(i, 4)+effrts(8)*fcttab(i, 5))
            phisup = fcttab(i, 1)*phixs+fcttab(i, 2)*phiys
            phiinf = fcttab(i, 1)*phixi+fcttab(i, 2)*phiyi

!       CALCUL DU FERRAILLAGE A L'ELU

            if (typcmb .eq. 0) then

!           CALCUL DES ACIERS DE FLEXION A L'ELU

                call cafelu(typco, alphacc, effm, effn, ht, 1.0, &
                            enrobi, enrobs, facier, fbeton, gammas, gammac, &
                            clacier, eys, typdiag, ferrcomp, precs, ferrsyme, slsyme, &
                            uc, um, &
                            ai(i), as(i), sigmsi, sigmss, ecinf, ecsup, &
                            alpha, pivot, etat, ierr)
!
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELU

                if (ierr .eq. 1) then

!               Facette en pivot B ou C trop comprimée !
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 1001

                end if

                if (ierr .eq. 2) then

!               Ferraillage symétrique non possible!
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 10011

                end if
!
!           CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELU

                if (ierrl .eq. 1) then
                    ai(i) = 0
                    as(i) = 0
                    alpha = -1
                    sigmsi = -facier/gammas
                    sigmss = -sigmsi
                end if

!               Calcul si aucune alarne émise pour les aciers de flexion
                call cftelu(typco, 0, effrts, effm, effn, efft, 0.0, &
                            ai(i), as(i), sigmsi, sigmss, alpha, &
                            ht, 1.0, enrobi, enrobs, facier, fbeton, &
                            alphacc, gammac, gammas, uc, um, &
                            compress, dnstra(i), thetab, ak, uk, ierr)
!
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELU

                if (ierr .eq. 1) then

!                   Béton trop cisaillé !
!                   Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'élément
                    ierrt = 1
                    ierrtG = 1002

                end if

!               Ajout de l'epure (impact de l'effort tranchant sur le ferr longi)

                if ((ierrl .eq. 0) .and. (ierrt .eq. 0) &
                    & .and. (dnstra(i) .gt. 0) .and. (epucisa .eq. 1)) then

                    Asl = abs(efft)/(tan(thetab))

                    if (effm .ge. 0) then

                        if (sigmsi .ne. (-1)) then
                            Asl = Asl/(-sigmsi)
                        else
                            Asl = Asl/(facier/gammas)
                        end if
                        ai(i) = max(ai(i)+Asl, 0.0)

                    else

                        if (sigmss .ne. (-1)) then
                            Asl = Asl/(-sigmss)
                        else
                            Asl = Asl/(facier/gammas)
                        end if
                        as(i) = max(as(i)+Asl, 0.0)

                    end if

                end if

!       CALCUL DU FERRAILLAGE A L'ELS

            elseif (typcmb .eq. 1) then

!           CALCUL DES ACIERS DE FLEXION A L'ELS

                call cafels(cequi, effm, effn, ht, 1.0, &
                            enrobi, enrobs, sigci, sigcs, sigs, &
                            ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                            ai(i), as(i), sigmsi, sigmss, &
                            sigmci, sigmcs, &
                            alpha, pivot, etat, ierr)

!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELS

                if (ierr .eq. 1) then

!               Facette en pivot B trop comprimée !
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 1003

                end if

                if (ierr .eq. 2) then

!               Ferraillage symétrique non possible!
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 10011

                end if
!
!           CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELS

!               Calcul si aucune alarne émise pour les aciers de flexion

                if (ierrl .eq. 1) then

                    ai(i) = 0
                    as(i) = 0
                    alpha = -1
                    sigmsi = -sigs
                    sigmss = -sigs
                    sigmci = sigci
                    sigmcs = sigcs

                end if

                call cftels(typco, 0, effrts, effm, effn, efft, 0.0, &
                            ai(i), as(i), &
                            sigmsi, sigmss, sigmci, sigmcs, alpha, &
                            ht, 1.0, enrobi, enrobs, facier, fbeton, &
                            sigci, sigcs, sigs, uc, um, &
                            compress, dnstra(i), thetab, ak, uk, ierr)
!
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELS

                if (ierr .eq. 1) then
!                   Béton trop cisaillé !
!                   Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'élément
                    ierrt = 1
                    ierrtG = 1004
                end if

!               Ajout de l'epure (impact de l'effort tranchant sur le ferr longi)

                if ((ierrl .eq. 0) .and. (ierrt .eq. 0) &
                    & .and. (dnstra(i) .gt. 0) .and. (epucisa .eq. 1)) then

                    Asl = abs(efft)/(tan(thetab))

                    if (effm .ge. 0) then

                        if (sigmsi .ne. (-1)) then
                            Asl = Asl/(-sigmsi)
                        else
                            Asl = Asl/sigs
                        end if
                        ai(i) = max(ai(i)+Asl, 0.0)

                    else

                        if (sigmss .ne. (-1)) then
                            Asl = Asl/(-sigmss)
                        else
                            Asl = Asl/sigs
                        end if
                        as(i) = max(as(i)+Asl, 0.0)

                    end if

                end if

!       CALCUL DU FERRAILLAGE A L'ELS QP

            elseif (typcmb .eq. 2) then

!           CALCUL DES ACIERS DE FLEXION A L'ELS QP

                call cafelsqp(cequi, effm, effn, ht, 1.0, &
                              enrobi, enrobs, wmaxi, wmaxs, &
                              ferrcomp, precs, ferrsyme, slsyme, uc, um, &
                              kt, facier, fbeton, eys, sigelsqp, phiinf, phisup, &
                              ai(i), as(i), sigmsi, sigmss, sigmci, sigmcs, &
                              alpha, pivot, etat, &
                              wfini, wfins, kvarf, ierr)

!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELS QP

                if (ierr .eq. 1) then

!               Facette en pivot B trop comprimée !
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 1005

                end if

                if (ierr .eq. 2) then

!               Ferraillage symétrique non possible!
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 10011

                end if

                if (ierr .eq. 3) then

!               Résolution itérative impossible à l'els qp !
!               Alarme dans te0146 + on sort de la boucle + densité = -1 pour l'élément
                    ierrl = 1
                    ierrlG = 1006

                end if

!           CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELS QP

!               Calcul si aucune alarne émise pour les aciers de flexion
                if (ierrl .eq. 1) then
                    ai(i) = 0
                    as(i) = 0
                    alpha = -1
                    sigmsi = -facier
                    sigmss = -facier
                    sigmci = sigelsqp
                    sigmcs = sigelsqp
                    kvarf = 1
                end if

                Sacier = kvarf*facier
                call cftels(typco, 0, effrts, effm, effn, efft, 0.0, &
                            ai(i), as(i), &
                            sigmsi, sigmss, sigmci, sigmcs, alpha, &
                            ht, 1.0, enrobi, enrobs, facier, fbeton, &
                            sigelsqp, sigelsqp, Sacier, uc, um, &
                            compress, dnstra(i), thetab, ak, uk, ierr)

!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELS QP

                if (ierr .eq. 1) then

!                   Béton trop cisaillé !
!                   Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'élément
                    ierrt = 1
                    ierrtG = 1007

                end if

!               Ajout de l'epure (impact de l'effort tranchant sur le ferr longi)

                if ((ierrl .eq. 0) .and. (ierrt .eq. 0) &
                    & .and. (dnstra(i) .gt. 0) .and. (epucisa .eq. 1)) then

                    Asl = abs(efft)/(tan(thetab))

                    if (effm .ge. 0) then

                        if (sigmsi .ne. (-1)) then
                            Asl = Asl/(-sigmsi)
                        else
                            Asl = Asl/Sacier
                        end if
                        ai(i) = max(ai(i)+Asl, 0.0)

                    else

                        if (sigmss .ne. (-1)) then
                            Asl = Asl/(-sigmss)
                        else
                            Asl = Asl/Sacier
                        end if
                        as(i) = max(as(i)+Asl, 0.0)

                    end if

                end if

            end if

!  -- VERIFICATION DU FERRAILLAGE MINIMUM :
!  ----------------------------------------

            if ((ferrmin .eq. 1) .or. (ferrmin .eq. 2)) then

                if (uc .eq. 0) then
                    unite_pa = 1.e6
                    unite_m = 1.
                elseif (uc .eq. 1) then
                    unite_pa = 1.
                    unite_m = 1.e-3
                end if

                if (fbeton .le. (50*unite_pa)) then
                    fctm = 0.30*((fbeton/unite_pa)**(2.0/3.0))
                else
                    fctm = 2.12*LOG(1.0+((fbeton/unite_pa)+8.0)/10.0)
                end if

                if (ferrmin .eq. 2) then
                    rholmin = max(0.26*(fctm/facier), 0.0013)
                    rhotmin = 0
                end if

                !ferraillage inferieur
                d = ht-enrobi
                if ((ai(i) .lt. (rholmin*d)) .and. (ierrl .eq. 0)) then
                    ai(i) = rholmin*d
                end if

                !ferraillage supérieur
                d = ht-enrobs
                if ((as(i) .lt. (rholmin*d)) .and. (ierrl .eq. 0)) then
                    as(i) = rholmin*d
                end if

                if ((dnstra(i) .lt. (rhotmin*ht)) .and. (ierrt .eq. 0)) then
                    dnstra(i) = rhotmin*ht
                end if

            end if

10          continue
!
!   OPTIMISATION DES FERRAILLAGES

!   FER2_R =  DNSXI DNSXS DNSYI DNSYS DNSXT DNSYT DNSVOL CONSTRUC
!               1     2     3     4     5     6      7       8

            if (ierrlG .eq. 0) then
                call clcopt(nb, fcttab, ai, dnsits(1), dnsits(2))
                call clcopt(nb, fcttab, as, dnsits(3), dnsits(4))
            else
                dnsits(1) = -1.d0
                dnsits(2) = -1.d0
                dnsits(3) = -1.d0
                dnsits(4) = -1.d0
            end if

            ierrl = ierrlG

            if (ierrtG .eq. 0) then
                call clcopt(nb, fcttab, dnstra, dnsits(5), dnsits(6))
            else
                dnsits(5) = -1.d0
                dnsits(6) = -1.d0
            end if

            ierrt = ierrtG

            end subroutine
