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
subroutine clcplq(typcmb, typco, ferrsyme, slsyme, ferrcomp,&
                  epucisa, ferrmin, rholmin, rhotmin, compress,&
                  cequi, enrobi, enrobs, sigs, sigci,&
                  sigcs, alphacc, gammas, gammac, facier,&
                  eys, typdiag, fbeton, clacier, uc,&
                  um, wmaxi, wmaxs, sigelsqp, kt,&
                  phixi, phixs, phiyi, phiys, ht,&
                  effrts, dnsits, ierr)
!
!_____________________________________________________________________
!
!      CLCPLQ
!
!      CALCUL DES ARMATURES EN ACIER DANS LES ELEMENTS DE PLAQUE
!
!      I TYPCMB        TYPE DE COMBINAISON (0 = ELU, 1 = ELS)
!      I TYPCO         CODIFICATION UTILISEE (1 = BAEL91, 2 = EC2)
!      I FERRSYME      FERRAILLAGE SYMETRIQUE?
!                      (0 = NON, 1 = OUI)
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
!      I TYPDIAG       TYPE DE DIAGRAMME UTILIS?? POUR L'ACIER
!                            TYPDIAG = 1 ("B1" ==> PALIER INCLIN??)
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
!      I WMAXI         OUVERTURE MAXIMALE DES FISSURES EN FACE INF??RIEURE
!      I WMAXS         OUVERTURE MAXIMALE DES FISSURES EN FACE SUP??RIEURE
!      I SIGELSQP      CONTRAINTE ADMISSIBLE DANS LE BETON ?? L'ELS QP
!      I KT            COEFFICIENT DE DUR??E DE CHARGEMENT
!      I PHIXI         DIAM??TRE APPROXIMATIF DES ARMATURES INF??RIEURES SUIVANT X
!      I PHIXS         DIAM??TRE APPROXIMATIF DES ARMATURES SUP??RIEURES SUIVANT X
!      I PHIYI         DIAM??TRE APPROXIMATIF DES ARMATURES INF??RIEURES SUIVANT Y
!      I PHIYS         DIAM??TRE APPROXIMATIF DES ARMATURES SUP??RIEURES SUIVANT Y
!      I HT            EPAISSEUR DE LA COQUE
!      I EFFRTS        (DIM 8) TORSEUR DES EFFORTS, MOMENTS, ...
!
!      O DNSITS        (DIM 6) DENSITES
!                            1..4 : SURFACES D'ACIER LONGITUDINAL
!                            5..6 : TRANSVERSAL
!      O IERR          CODE RETOUR (0 = OK)
!
!_____________________________________________________________________
!
!
    implicit none
!
!
#include "asterfort/cafels.h"
#include "asterfort/cafelsqp.h"
#include "asterfort/cafelu.h"
#include "asterfort/cftels.h"
#include "asterfort/cftelu.h"
#include "asterfort/clcopt.h"
#include "asterfort/trgfct.h"
#include "asterfort/utmess.h"
!
!
    integer :: typcmb
    integer :: typco
    integer :: ferrsyme
    real(kind=8) :: slsyme
    integer :: ferrcomp
    integer :: epucisa
    integer :: ferrmin
    real(kind=8) :: rholmin
    real(kind=8) :: rhotmin
    integer :: compress
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
    integer :: typdiag
    real(kind=8) :: fbeton
    integer :: clacier
    integer :: uc
    integer :: um
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
    integer :: ierr
!
!
!       NOMBRE DE DIVISIONS ENTRE -PI/2 ET +PI/2
    real(kind=8) :: fcttab(36, 5)
!       NOMBRE DE FACETTES COMPRIMEES EN PIVOT C (ELU ET ELS)
    real(kind=8) :: nb_fac_comp_elu, nb_fac_comp_els
!       EFFORT NORMAL DANS CETTE DIRECTION
    real(kind=8) :: effn
!       EFFORT TRANCHANT DANS CETTE DIRECTION
    real(kind=8) :: efft
!       MOMENT DE FLEXION DANS CETTE DIRECTION
    real(kind=8) :: effm
!       DIAMETRE DES BARRES SUP??RIEURES DANS CETTE DIRECTION
    real(kind=8) :: phisup
!       DIAMETRE DES BARRES INF??RIEURES DANS CETTE DIRECTION
    real(kind=8) :: phiinf
!       DENSITE DE FERRAILLAGE TRANSVERSAL
    real(kind=8) :: dnstra(36)
!       SECTIONS DES ACIERS INFERIEURS SUIVANT LES 36 FACETTES
    real(kind=8) :: ai(36)
!       SECTIONS DES ACIERS SUPERIEURS SUIVANT LES 36 FACETTES
    real(kind=8) :: as(36)
!       CONTRAINTE DANS L'ACIER DE TRACTION A L'ELS QP
    real(kind=8) :: Sacier
!       VARIABLE D'ITERATION
    integer :: i
!       AUTRES VARIABLES
    real(kind=8) :: fctm, unite_m, unite_pa, d, ak, uk, alpha, thetab
    real(kind=8) :: Asl, ecinf, ecsup, kvarf, sigmci, sigmcs, sigmsi, sigmss
    real(kind=8) :: wfini, wfins
    integer :: etat, pivot
!
!
!   INITIALISATION DES VARIABLES
!
    ierr = 0
    nb_fac_comp_elu = 0
    nb_fac_comp_els = 0
!
!   INITIALISATION DES FACETTES
!
    call trgfct(fcttab)
    do i = 1, 6
        dnsits(i) = -1.d0
    end do
!
!   BOUCLE SUR LES FACETTES DE CAPRA ET MAURY
!   DETERMINATION DU FERRAILLAGE POUR CHACUNE DES FACETTES
!
!print *,"effrts = ",effrts
    do i = 1, 36
!print *,"i = ",i
        effn = fcttab(i,1) * effrts(1) + fcttab(i,2) * effrts(2) + fcttab(i,3) * effrts(3)
        effm = fcttab(i,1) * effrts(4) + fcttab(i,2) * effrts(5) + fcttab(i,3) * effrts(6)
        effn = -effn
        effm = -effm
!
        efft = abs(effrts(7)*fcttab(i,4) + effrts(8)*fcttab(i,5))
        phisup = fcttab(i,1) * phixs + fcttab(i,2) * phiys
        phiinf = fcttab(i,1) * phixi + fcttab(i,2) * phiyi
!
!       CALCUL DU FERRAILLAGE A L'ELU
!
        if (typcmb .eq. 0) then
!
!           CALCUL DES ACIERS DE FLEXION A L'ELU
            call cafelu(typco, alphacc, effm, effn, ht,&
                        1.0, enrobi, enrobs, facier, fbeton,&
                        gammas, gammac, clacier, eys, typdiag,&
                        ferrcomp, ferrsyme, slsyme, uc, ai(i),&
                        as(i), sigmsi, sigmss, ecinf, ecsup,&
                        alpha, pivot, etat, ierr)
!
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELU
            if (ierr .eq. 1) then
!               Facette en pivot B ou C trop comprim??e !
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 1001
                goto 999
            endif
!
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELU
            if (ierr .eq. 2) then
!               Ferraillage sym??trique non possible!
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 10011
                goto 999
            endif
!
!           CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELU
            if (ierr .eq. 0) then
!               Calcul si aucune alarne ??mise pour les aciers de flexion
                call cftelu(typco, 0, effrts, effm, effn,&
                            efft, 0.0, ai(i), as(i), sigmsi,&
                            sigmss, alpha, ht, 1.0, enrobi,&
                            enrobs, facier, fbeton, alphacc, gammac,&
                            gammas, uc, um, compress, dnstra(i),&
                            thetab, ak, uk, ierr)
!
!
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELU
                if (ierr .eq. 1) then
!                   B??ton trop cisaill?? !
!                   Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'??l??ment
                    ierr = 1002
                    goto 999
                endif
!
!               Ajout de l'epure (impact de l'effort tranchant sur le ferr longi)
                if ((ierr.eq.0) .and. (epucisa.eq.1)) then
                    Asl = abs(efft)/(tan(thetab))
                    if (effm .ge. 0) then
                        if (sigmsi .ne. (-1)) then
                            Asl = Asl/(-sigmsi)
                        else
                            Asl = Asl/(facier/gammas)
                        endif
                        ai(i) = max(ai(i) + Asl,0.0)
                    else
                        if (sigmss .ne. (-1)) then
                            Asl = Asl/(-sigmss)
                        else
                            Asl = Asl/(facier/gammas)
                        endif
                        as(i) = max(as(i) + Asl,0.0)
                    endif
                endif
!
            endif
!
!       CALCUL DU FERRAILLAGE A L'ELS
!
        else if (typcmb .eq. 1) then
!
!           CALCUL DES ACIERS DE FLEXION A L'ELS
            call cafels(cequi, effm, effn, ht, 1.0,&
                        enrobi, enrobs, sigci, sigcs, sigs,&
                        ferrcomp, ferrsyme, slsyme, uc, ai(i),&
                        as(i), sigmsi, sigmss, sigmci, sigmcs,&
                        alpha, pivot, etat, ierr)
!
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELS
            if (ierr .eq. 1) then
!               Facette en pivot B trop comprim??e !
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 1003
                goto 999
            endif
!
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELU
            if (ierr .eq. 2) then
!               Ferraillage sym??trique non possible!
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 10011
                goto 999
            endif
!
!           CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELS
            if (ierr .eq. 0) then
!               Calcul si aucune alarne ??mise pour les aciers de flexion
                call cftels(typco, 0, effrts, effm, effn,&
                            efft, 0.0, ai(i), as(i), sigmsi,&
                            sigmss, sigmci, sigmcs, alpha, ht,&
                            1.0, enrobi, enrobs, facier, fbeton,&
                            sigci, sigcs, sigs, uc, um,&
                            compress, dnstra(i), thetab, ak, uk,&
                            ierr)
!
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELS
                if (ierr .eq. 1) then
!                   B??ton trop cisaill?? !
!                   Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'??l??ment
                    ierr = 1004
                    goto 999
                endif
!
!               Ajout de l'epure (impact de l'effort tranchant sur le ferr longi)
                if ((ierr.eq.0) .and. (dnstra(i).gt.0) .and. (epucisa.eq.1)) then
                    Asl = abs(efft)/(tan(thetab))
                    if (effm .ge. 0) then
                        if (sigmsi .ne. (-1)) then
                            Asl = Asl/(-sigmsi)
                        else
                            Asl = Asl/sigs
                        endif
                        ai(i) = max(ai(i) + Asl,0.0)
                    else
                        if (sigmss .ne. (-1)) then
                            Asl = Asl/(-sigmss)
                        else
                            Asl = Asl/sigs
                        endif
                        as(i) = max(as(i) + Asl,0.0)
                    endif
                endif
!
            endif
!
!       CALCUL DU FERRAILLAGE A L'ELS QP
!
        else if (typcmb .eq. 2) then
!
!           CALCUL DES ACIERS DE FLEXION A L'ELS QP
            call cafelsqp(cequi, effm, effn, ht, 1.0,&
                          enrobi, enrobs, wmaxi, wmaxs, ferrcomp,&
                          ferrsyme, slsyme, uc, um, kt,&
                          facier, fbeton, eys, sigelsqp, phiinf,&
                          phisup, ai(i), as(i), sigmsi, sigmss,&
                          sigmci, sigmcs, alpha, pivot, etat,&
                          wfini, wfins, kvarf, ierr)
!
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELS QP
            if (ierr .eq. 1) then
!               Facette en pivot B trop comprim??e !
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 1005
                goto 999
            endif
!           GESTION DES ALARMES EMISES POUR LES ACIERS DE FLEXION A L'ELU
            if (ierr .eq. 2) then
!               Ferraillage sym??trique non possible!
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 10011
                goto 999
            endif
            if (ierr .eq. 3) then
!               R??solution it??rative impossible ?? l'els qp !
!               Alarme dans te0146 + on sort de la boucle + densit?? = -1 pour l'??l??ment
                ierr = 1006
                goto 999
            endif
!
!           CALCUL DU FERRAILLAGE TRANSVERSAL A L'ELS QP
            if (ierr .eq. 0) then
!               Calcul si aucune alarne ??mise pour les aciers de flexion
                Sacier = kvarf*facier
                call cftels(typco, 0, effrts, effm, effn,&
                            efft, 0.0, ai(i), as(i), sigmsi,&
                            sigmss, sigmci, sigmcs, alpha, ht,&
                            1.0, enrobi, enrobs, facier, fbeton,&
                            sigelsqp, sigelsqp, Sacier, uc, um,&
                            compress, dnstra(i), thetab, ak, uk,&
                            ierr)
!
!               GESTION DES ALARMES EMISES POUR LE FERRAILLAGE TRANSVERSAL A L'ELS
                if (ierr .eq. 1) then
!                   B??ton trop cisaill?? !
!                   Alarme dans te0146 + on sort de la boucle + dnstra = -1 pour l'??l??ment
                    ierr = 1007
                    goto 999
                endif
!
!               Ajout de l'epure (impact de l'effort tranchant sur le ferr longi)
                if ((ierr.eq.0) .and. (dnstra(i).gt.0) .and. (epucisa.eq.1)) then
                    Asl = abs(efft)/(tan(thetab))
                    if (effm .ge. 0) then
                        if (sigmsi .ne. (-1)) then
                            Asl = Asl/(-sigmsi)
                        else
                            Asl = Asl/Sacier
                        endif
                        ai(i) = max(ai(i) + Asl,0.0)
                    else
                        if (sigmss .ne. (-1)) then
                            Asl = Asl/(-sigmss)
                        else
                            Asl = Asl/Sacier
                        endif
                        as(i) = max(as(i) + Asl,0.0)
                    endif
                endif
!
            endif
!
        endif
!
!  -- VERIFICATION DU FERRAILLAGE MINIMUM :
!  ----------------------------------------
!
        if ((ferrmin.eq.1) .or. (ferrmin.eq.2)) then
!
            if (uc .eq. 0) then
                unite_pa = 1.e6
                unite_m = 1.
            else if (uc.eq.1) then
                unite_pa = 1.
                unite_m = 1.e-3
            endif
!
            if (fbeton .le. (50*unite_pa)) then
                fctm = 0.30*((fbeton/unite_pa)**(2.0/3.0))
            else
                fctm = 2.12*LOG(1.0+((fbeton/unite_pa)+8.0)/10.0)
            endif
!
            if (ferrmin .eq. 2) then
                rholmin = max(0.26*(fctm/facier),0.0013)
                rhotmin = 0
            endif
!
!ferraillage inferieur
            d = ht - enrobi
            if ((ai(i).lt.(rholmin*d)) .and. (ierr.ne.1001) .and. (ierr.ne.1003) .and.&
                (ierr.ne.1005) .and. (ierr.ne.1006)) then
                ai(i) = rholmin*d
            endif
!
!ferraillage sup??rieur
            d = ht - enrobs
            if ((as(i).lt.(rholmin*d)) .and. (ierr.ne.1001) .and. (ierr.ne.1003) .and.&
                (ierr.ne.1005) .and. (ierr.ne.1006)) then
                as(i) = rholmin*d
            endif
!
            if ((dnstra(i).lt.(rhotmin*ht)) .and. (ierr.ne.1002) .and. (ierr.ne.1004) .and.&
                (ierr.ne.1007)) then
                dnstra(i) = rhotmin*ht
            endif
!
        endif
    end do
!
!   OPTIMISATION DES FERRAILLAGES
!
!   FER2_R =  DNSXI DNSXS DNSYI DNSYS DNSXT DNSYT DNSVOL CONSTRUC
!               1     2     3     4     5     6      7       8
!
    call clcopt(fcttab, ai, dnsits(1), dnsits(2))
    call clcopt(fcttab, as, dnsits(3), dnsits(4))
    call clcopt(fcttab, dnstra, dnsits(5), dnsits(6))
!
999 continue
!
end subroutine
