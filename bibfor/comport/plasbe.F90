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
subroutine plasbe(BEHinteg, &
                  fami, kpg, ksp, typmod, imat, l_epsi_varc, &
                  crit, epsdt, depst, sigd, vind, &
                  opt, sigf, vinf, dsde, &
                  icomp, nvi, irteti)
!
    use Behaviour_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/betcvx.h"
#include "asterfort/betimp.h"
#include "asterfort/betjpl.h"
#include "asterfort/betmat.h"
#include "asterfort/codent.h"
#include "asterfort/lcdedi.h"
#include "asterfort/lcdehy.h"
#include "asterfort/lcelin.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcplbe.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!

!       ================================================================
!       INTEGRATION DE LOIS DE COMPORTEMENT ELASTO PLASTIQUE ET VISCO
!       PLASTIQUE
!               AVEC    . N VARIABLES INTERNES
!                       . UNE FONCTION SEUIL ELASTIQUE
!
!       INTEGRATION DES CONTRAINTES           = SIG(T+DT)
!       INTEGRATION DES VARIABLES INTERNES    = VIN(T+DT)
!       ET CALCUL DU JACOBIEN ASSOCIE         = DS/DE(T+DT) OU DS/DE(T)
!
!       EN CAS DE NON-CONVERGENCE LOCALE ON EFFECTUE UN REDECOUPAGE DU
!       PAS DE TEMPS, L ORDRE D EXECUTION ETANT REMONTE EN ARGUMENT
!       DANS REDECE, APPELE PAR NMCOMP AVANT PLASBE
!       ================================================================
!       ROUTINE CONSTRUITE SUIVANT LE MODELE ET L'ARCHITECTURE DE
!                                PLASTI
!       ================================================================
!       ROUTINES UTILITAIRES DE CALCUL MATRICIEL(6,6) - VECTORIEL (6)
!       COMMON /TDIM/ NDT,NDI    A INCLURE OBLIGATOIREMENT
!       () = DEBUG
!
!       ( LCIMMA  IMPRESSION MATRICE )
!       ( LCIMVE  IMPRESSION VECTEUR )
!       ( LCIMSC  IMPRESSION SCALAIRE )
!       LCPRTE  PRODUIT TENSORIEL DE VECTEURS
!       LCPRSC  PRODUIT SCALAIRE  DE VECTEURS
!       LCPRSM  PRODUIT SCALAIRE * MATRICE
!       LCPTMV  PRODUIT MATRICE TRANSPOSEE * VECTEUR
!       LCPRMV  PRODUIT MATRICE  * VECTEUR
!       LCPRMM  PRODUIT MATRICE  * MATRICE
!       LCDIMA  DIFFERENCE DE MATRICES
!       LCSOMA  SOMME DE MATRICES
!       LCTRMA  TRANSPOSEE DE MATRICE
!       LCEQMA  EGALITE DE MATRICES
!       LCINMA  INITIALISATION DE MATRICE A UNE VALEUR
!       LCPRSV  PRODUIT SCALAIRE * VECTEUR
!       LCSOVE  SOMME DE VECTEUR
!       LCDIVE  DIFFERENCE DE VECTEURS
!       LCINVE  INITIALISATION DE VECTEUR
!       LCEQVE  EGALITE DE VECTEURS
!
!       ----------------------------------------------------------------
!       ROUTINES UTILITAIRES D INTEGRATION D UN MODELE DE COMPORTEMENT
!
!       LCDEVI  PARTIE DEVIATORIQUE D UN TENSEUR
!       LCHYDR  PARTIE SPHERIQUE    D UN TENSEUR
!       LCVS    DERIVEE DE LA CONTRAINTE VON MISES / CONTRAINTE
!       LCVVSS  DERIVEE SECONDE DE LA CONTRAINTE VON MISES / CONTRAINTE
!       LCIV2S  SECOND INVARIANT DU TENSEUR CONTRAINTE
!       LCIV2E  SECOND INVARIANT DU TENSEUR DEFORMATION
!       LCNRTS  'NORME' DU TENSEUR CONTRAINTE
!       LCNRTE  'NORME' DU TENSEUR DEFORMATION
!       LCOPLI  OPERATEUR ELASTIQUE LINEAIRE
!       LCELIN  INTEGRATION  ELASTIQUE LINEAIRE ISOTROPE
!       LCVERR  CALCUL DU VECTEUR D'ERREUR RELATIVE, ABSOLU, NORMEE...
!
!       ================================================================
!       ARGUMENTS
!
!       IN      TYPMOD  TYPE DE MODELISATION
!               IMAT    ADRESSE DU MATERIAU CODE
!               COMP    COMPORTEMENT DE L ELEMENT
!                       COMP(1) = RELATION DE COMPORTEMENT
!                       COMP(2) = NB DE VARIABLES INTERNES
!                       COMP(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!               OPT     OPTION DE CALCUL A FAIRE
!                               'RIGI_MECA_TANG'> DSDE(T)
!                               'FULL_MECA'     > DSDE(T+DT) , SIG(T+DT)
!                               'RAPH_MECA'     > SIG(T+DT)
!               CRIT    CRITERES  LOCAUX
!                       CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                                 (ITER_INTE_MAXI == ITECREL)
!                       CRIT(2) = TYPE DE JACOBIEN A T+DT
!                                 (TYPE_MATR_COMP == MACOMP)
!                                 0 = EN VITESSE     > SYMETRIQUE
!                                 1 = EN INCREMENTAL > NON-SYMETRIQUE
!                       CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                                 (RESI_INTE == RESCREL)
!                       CRIT(5) = NOMBRE D'INCREMENTS POUR LE
!                                 REDECOUPAGE LOCAL DU PAS DE TEMPS
!                                 (RESI_INTE_PAS == ITEDEC )
!                                 0 = PAS DE REDECOUPAGE
!                                 N = NOMBRE DE PALIERS
!               EPSDT   DEFORMATION TOTALE A T
!               DEPST   INCREMENT DE DEFORMATION TOTALE
!               SIGD    CONTRAINTE A T
!               VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!               ICOMP   COMPTEUR POUR LE REDECOUPAGE DU PAS DE TEMPS
!       OUT     SIGF    CONTRAINTE A T+DT
!               VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSDE    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!               IRTETI = 1 CONTROLE DU REDECOUPAGE DU PAS DE TEMPS
!       ----------------------------------------------------------------
!       INFO    MATERD        (*,1) = CARACTERISTIQUES ELASTIQUES A T
!                             (*,2) = CARACTERISTIQUES PLASTIQUES A T
!               MATERF        (*,1) = CARACTERISTIQUES ELASTIQUES A T+DT
!                             (*,2) = CARACTERISTIQUES PLASTIQUES A T+DT
!               MATCST          'OUI' SI MATERIAU CST ENTRE T ET T+DT
!                               'NON' SINON
!               NDT             NB DE COMPOSANTE TOTALES DES TENSEURS
!                                       = 6  3D
!                                       = 4  AXIS  C_PLAN  D_PLAN
!                                       = 1  1D
!               NDI             NB DE COMPOSANTE DIRECTES DES TENSEURS
!               NVI             NB DE VARIABLES INTERNES
!               NR              NB EQUATION SYSTEME INTEGRE A RESOUDRE
!       ----------------------------------------------------------------
!       ROUTINE LC....UTILITAIRES POUR INTEGRATION LOI DE COMPORTEMENT
!       ----------------------------------------------------------------
!       ORDRE DES TENSEURS      3D      XX YY ZZ XY XZ YZ
!                               DP      XX YY ZZ XY
!                               AX      RR ZZ TT RZ
!                               1D      XX YY ZZ
!       ----------------------------------------------------------------
!       ATTENTION
!       SI OPT = 'RIGI_MECA_TANG' NE PAS TOUCHER AUX VARIABLES SIGF,VINF
!       QUI N ONT PAS DE PLACE MEMOIRE ALLOUEE
!
!       SIG EPS DEPS  ONT DEJA LEURS COMPOSANTES DE CISAILLEMENT
!       MULTIPLIES PAR RACINE DE 2 > PRISE EN COMPTE DES DOUBLES
!       PRODUITS TENSORIELS ET CONSERVATION DE LA SYMETRIE
!
!       ----------------------------------------------------------------
    type(Behaviour_Integ), intent(in) :: BEHinteg
    aster_logical, intent(in) :: l_epsi_varc
    integer(kind=8) :: imat, ndt, ndi, nr, nvi
    integer(kind=8) :: itmax, icomp
    integer(kind=8) :: nmat, irtet, irteti, nseui4
    integer(kind=8) :: nseuil, nseui1, nseui2, nseui3
    integer(kind=8) :: iadzi, iazk24
    real(kind=8) :: toler
    real(kind=8) :: epsi
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret, kpg, ksp
    real(kind=8) :: tneg, tref
!-----------------------------------------------------------------------
    parameter(epsi=1.d-15)
    parameter(nmat=90)
    parameter(tneg=-1.d3)
!
    real(kind=8) :: crit(*)
    real(kind=8) :: vind(*), vinf(*)
    real(kind=8) :: tempd, tempf
    real(kind=8) :: epsd(6), deps(6)
    real(kind=8) :: epsdt(6), depst(6)
    real(kind=8) :: sigd(6), sigf(6), sige(6)
!
    real(kind=8) :: dsde(6, 6)
    real(kind=8) :: materd(nmat, 2), materf(nmat, 2), tmpmx
    character(len=7) :: etatd, etatf
    character(len=8) :: mod, typma, typmod(*)
    character(len=16) :: opt
    character(len=3) :: matcst, cnseui
    character(len=8) :: nomail
    character(len=*) :: fami
    aster_logical :: rigi, resi, istemp
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
    common/ecri/nomail
!       ----------------------------------------------------------------
!
! --    INITIALISATION DES PARAMETRES DE CONVERGENCE ET ITERATIONS
!
    irteti = 0
    itmax = int(crit(1))
    toler = crit(3)
!        LOI      = COMP(1)
    mod = typmod(1)
    nseuil = 0
    nseui1 = 0
    nseui2 = 0
    nseui3 = 0
    nseui4 = 0
    nomail = ' '
!
    resi = opt(1:9) .eq. 'FULL_MECA' .or. opt .eq. 'RAPH_MECA'
    rigi = opt(1:9) .eq. 'FULL_MECA' .or. opt(1:9) .eq. 'RIGI_MECA'
    ASSERT((opt(1:9) .eq. 'RIGI_MECA') .or. (opt(1:9) .eq. 'FULL_MECA') .or. (opt .eq. 'RAPH_MECA'))
    call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                ksp, tempd, iret)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, tempf, iret)
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, &
                ksp, tref, iret)
!
! --    C'EST INTERDIT DE MELANGER DEUX MODELISATIONS AVEC OU SANS
! --      DEPENDENCE DES PARAMETRES DE LA TEMPERATURE
    if ((.not. isnan(tempd)) .and. (isnan(tempf))) then
        call utmess('F', 'ALGORITH9_100')
    else if ((isnan(tempd)) .and. (.not. isnan(tempf))) then
        call utmess('F', 'ALGORITH9_100')
    else if ((vind(3) .eq. tneg) .and. (.not. isnan(tempf))) then
        call utmess('F', 'ALGORITH9_100')
    else
        istemp = .not. isnan(tempd) .and. .not. isnan(tempf)
    end if
!
! --    OPTION SUPPRIMEE CAR TYPMA EST IMPOSE SUIVANT QUE L'ON EST EN
! --    PLASTCITE OU VISCOPLASTICITE. TYPMA EST DEFINI DANS LCMATE
! --    POUR LES MODELE VISCO-PLASTIQUES A LA VALEUR 'COHERENT'.
!          IF ( INT(CRIT(2)) .EQ. 0 ) THEN
!          TYPMA = 'VITESSE '
!          ELSE
!          TYPMA = 'COHERENT'
!          ENDIF
!
    typma = 'VITESSE '
!
    if (itmax .le. 0) itmax = -itmax
!
! --    LES PARAMETRES SONT FONCTIONS DE LA TEMPERATURE MAXIMALE
! --    VIND(3) EST LE MAX DES TEMPERATURES DANS L'HISTORIQUE DES TEMP.
!
    if (istemp) then
        tmpmx = vind(3)
        if (tempd .gt. tmpmx) tmpmx = tempd
    else
        tmpmx = r8nnem()
    end if
!
! --    RECUPERATION COEF(TEMP(T))) LOI ELASTO-PLASTIQUE A T ET/OU T+DT
!                    NB DE CMP DIRECTES/CISAILLEMENT + NB VAR. INTERNES
!
!
    call betmat(fami, kpg, ksp, mod, imat, &
                nmat, tmpmx, tempf, materd, materf, &
                matcst, ndt, ndi, nr, nvi)
!
! --    RETRAIT INCREMENT DE DEFORMATION DUE A LA DILATATION THERMIQUE
!
    call lcdedi(fami, kpg, ksp, nmat, materd, &
                materf, tempd, tempf, tref, depst, &
                epsdt, deps, epsd, l_epsi_varc)
!
! --    RETRAIT ENDOGENNE ET RETRAIT DE DESSICCATION
!
    call lcdehy(fami, kpg, ksp, nmat, materd, &
                materf, deps, epsd, l_epsi_varc)
!
! --    SEUIL A T > ETAT ELASTIQUE OU PLASTIQUE A T
!
    if (abs(vind(nvi)) .le. epsi) then
        etatd = 'ELASTIC'
    else
        etatd = 'PLASTIC'
    end if
!
!  -->  REDECOUPAGE IMPOSE
    if (icomp .eq. -1 .and. opt .ne. 'RIGI_MECA_TANG') then
        irteti = 0
        goto 999
    end if
!
!       ----------------------------------------------------------------
!       OPTIONS 'FULL_MECA' ET 'RAPH_MECA' = CALCUL DE SIG(T+DT)
!       ----------------------------------------------------------------
!
    if (resi) then
!
! --    INTEGRATION ELASTIQUE SUR DT
!
        call lcelin(mod, nmat, materd, materf, deps, &
                    sigd, sige)
        vinf(1:nvi-1) = vind(1:nvi-1)
        vinf(nvi) = 0.d0
!
        if (istemp) then
            if (tmpmx .lt. tempf) tmpmx = tempf
            vinf(3) = tmpmx
        else
            vinf(3) = tneg
        end if
!
        sigf(1:ndt) = sige(1:ndt)
!
! --    PREDICTION ETAT ELASTIQUE A T+DT : F(SIG(T+DT),VIN(T)) = 0 ?
!
        call betcvx(BEHinteg, nmat, materf, sigf, vind, vinf, &
                    nvi, nseuil)
!
        if (nseuil .ge. 0) then
!
! --       PREDICTION INCORRECTE > INTEGRATION ELASTO-PLASTIQUE SUR DT
!
            etatf = 'PLASTIC'
!
            nseui1 = nseuil
            call lcplbe(BEHinteg, toler, itmax, nmat, materf, nvi, &
                        vind, sigf, vinf, nseuil, &
                        irtet)
!           GOTO (1), IRTET
!
            call betcvx(BEHinteg, nmat, materf, sigf, vind, vinf, &
                        nvi, nseuil)
            nseui2 = nseuil
!
            if (nseui2 .gt. 0) then
                if (nseui2 .eq. 44) then
                    call utmess('A', 'ALGORITH9_93')
                    goto 1
                end if
                if (nseui2 .eq. 4) then
                    call codent(nseui1, 'G', cnseui)
                    call utmess('A', 'ALGORITH9_94', sk=cnseui)
                    goto 1
                end if
                if (nseui2 .eq. nseui1) then
                    if (nseui2 .ne. 3) then
                        nseui2 = 3
                    else
                        nseui2 = 2
                    end if
                    nseuil = nseui2
                end if
                sigf(1:ndt) = sige(1:ndt)
                call lcplbe(BEHinteg, toler, itmax, nmat, materf, nvi, &
                            vind, sigf, vinf, nseuil, &
                            irtet)
!              GOTO (1), IRTET
!
                call betcvx(BEHinteg, nmat, materf, sigf, vind, vinf, &
                            nvi, nseuil)
                nseui3 = nseuil
            end if
!
            if (nseui3 .gt. 0) then
                if (nseui3 .eq. 44) then
                    call utmess('A', 'ALGORITH9_93')
                    goto 1
                end if
                if (nseui3 .eq. 4) then
                    call codent(nseui2, 'G', cnseui)
                    call utmess('A', 'ALGORITH9_94', sk=cnseui)
                    goto 1
                end if
                if (nseui3 .eq. nseui1 .or. nseui3 .eq. nseui2) then
                    nseui3 = 6-nseui1-nseui2
                    nseuil = nseui3
                end if
                sigf(1:ndt) = sige(1:ndt)
                call lcplbe(BEHinteg, toler, itmax, nmat, materf, nvi, &
                            vind, sigf, vinf, nseuil, &
                            irtet)
!              GOTO (1), IRTET
!
                call betcvx(BEHinteg, nmat, materf, sigf, vind, vinf, &
                            nvi, nseuil)
                nseui4 = nseuil
            end if
!
            if (nseui4 .gt. 0) then
                if (nseui4 .eq. 44) then
                    call utmess('A', 'ALGORITH9_93')
                    goto 1
                end if
                if (nseui4 .eq. 4) then
                    call codent(nseui3, 'G', cnseui)
                    call utmess('A', 'ALGORITH9_94', sk=cnseui)
                    goto 1
                end if
                nseuil = 22
                nseui4 = nseuil
                sigf(1:ndt) = sige(1:ndt)
                call lcplbe(BEHinteg, toler, itmax, nmat, materf, nvi, &
                            vind, sigf, vinf, nseuil, &
                            irtet)
!             GOTO (1), IRTET
!
                call betcvx(BEHinteg, nmat, materf, sigf, vind, vinf, &
                            nvi, nseuil)
            end if
!
            if (nseuil .ge. 0) then
                call codent(nseui4, 'G', cnseui)
                call utmess('A', 'ALGORITH9_94', sk=cnseui)
                goto 1
            end if
!
        else
!
! --       PREDICTION CORRECTE > INTEGRATION ELASTIQUE FAITE
!
            etatf = 'ELASTIC'
        end if
!
    end if
!
!       ----------------------------------------------------------------
!       OPTIONS 'FULL_MECA', 'RIGI_MECA_TANG', 'RIGI_MECA_ELAS' :
!          CALCUL DE DSDE
!       ----------------------------------------------------------------
!       EVALUATION DU JACOBIEN DSDE A (T+DT) POUR 'FULL_MECA'
!       ET CALCUL ELASTIQUE    ET   A (T)    POUR 'RIGI_MECA_TANG'
!       ----------------------------------------------------------------
!
    if (rigi) then
        if (opt(1:9) .eq. 'RIGI_MECA') then
            if ((etatd .eq. 'ELASTIC') .or. (opt(10:14) .eq. '_ELAS')) then
                call lcopli('ISOTROPE', mod, materd(1, 1), dsde)
            else if (etatd .eq. 'PLASTIC') then
!   ------> ELASTOPLASTICITE ==> TYPMA = 'VITESSE '
!   ------> VISCOPLASTICITE  ==> TYPMA = 'COHERENT '==> CALCUL ELASTIQUE
                if (typma .eq. 'COHERENT') then
! PAS UTILISE ICI  CALL LCJELA ( LOI  , MOD ,  NMAT, MATERD,VIND, DSDE)
                else if (typma .eq. 'VITESSE ') then
                    call betjpl(BEHinteg, mod, nmat, materd, sigd, vind, &
                                dsde)
                end if
            end if
!
        else if (opt(1:9) .eq. 'FULL_MECA') then
            if ((etatf .eq. 'ELASTIC') .or. (opt(10:14) .eq. '_ELAS')) then
                call lcopli('ISOTROPE', mod, materf(1, 1), dsde)
            else if (etatf .eq. 'PLASTIC') then
!   ------> ELASTOPLASTICITE ==>  TYPMA = 'VITESSE '
!   ------> VISCOPLASTICITE  ==>  TYPMA = 'COHERENT '
                if (typma .eq. 'COHERENT') then
! PAS UTILISE ICI  CALL LCJPLC ( LOI  , MOD ,  NMAT, MATERD, DSDE)
                else if (typma .eq. 'VITESSE ') then
                    call betjpl(BEHinteg, mod, nmat, materd, sigf, vinf, &
                                dsde)
                end if
            end if
        end if
    end if
!
!       ----------------------------------------------------------------
!
    irteti = 0
    goto 999
1   continue
    irteti = 1
    call betimp(BEHinteg, nmat, materf, sigf, vind, vinf, &
                nseui1, nseui2, nseui3, nseui4, &
                sige, sigd)
!
    call tecael(iadzi, iazk24)
    nomail = zk24(iazk24-1+3) (1:8)
    call utmess('A', 'ALGORITH9_95', sk=nomail)
!
999 continue
!
end subroutine
