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
! aslint: disable=W0413

subroutine nmcjs(typmod, imat, crit, &
                 epsd, &
                 deps, sigd, vind, opt, sigf, &
                 vinf, dsde, iret)
    implicit none
!
!
!       INTEGRATION DE LA LOI DE COMPORTEMENT ELASTO PLASTIQUE CJS
!               AVEC    . N VARIABLES INTERNES
!                       . 2 FONCTIONS SEUIL ELASTIQUE
!
!       INTEGRATION DES CONTRAINTES           = SIG(T+DT)
!       INTEGRATION DES VARIABLES INTERNES    = VIN(T+DT)
!       ET CALCUL DU JACOBIEN ASSOCIE         = DS/DE(T+DT) OU DS/DE(T)
!       ================================================================
!       IN      TYPMOD  TYPE DE MODELISATION
!               IMAT    ADRESSE DU MATERIAU CODE
!               COMP    COMPORTEMENT DE L ELEMENT
!                       COMP(1) = RELATION DE COMPORTEMENT (CHABOCHE...)
!                       COMP(2) = NB DE VARIABLES INTERNES
!                       COMP(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
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
!               INSTAM  INSTANT T
!               INSTAP  INSTANT T+DT
!               EPSD    DEFORMATION TOTALE A T
!               DEPS    INCREMENT DE DEFORMATION TOTALE
!               SIGD    CONTRAINTE A T
!               VIND    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!               OPT     OPTION DE CALCUL A FAIRE
!                               'RIGI_MECA_TANG'> DSDE(T)
!                               'FULL_MECA'     > DSDE(T+DT) , SIG(T+DT)
!                               'RAPH_MECA'     > SIG(T+DT)
!       OUT     SIGF    CONTRAINTE A T+DT
!               VINF    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSDE    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!               IRET    CODE RETOUR DE  L'INTEGRATION DE LA LOI CJS
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ECHEC
!       ----------------------------------------------------------------
!       INFO    MATERD        (*,1) = CARACTERISTIQUES ELASTIQUES A T
!                             (*,2) = CARACTERISTIQUES PLASTIQUES A T
!               MATERF        (*,1) = CARACTERISTIQUES ELASTIQUES A T+DT
!                             (*,2) = CARACTERISTIQUES PLASTIQUES A T+DT
!               NDT             NB DE COMPOSANTE TOTALES DES TENSEURS
!                                       = 6  3D
!                                       = 4  AXIS  C_PLAN  D_PLAN
!                                       = 1  1D
!               NDI             NB DE COMPOSANTE DIRECTES DES TENSEURS
!               NVI             NB DE VARIABLES INTERNES
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
#include "asterf_types.h"
#include "asterfort/cjsela.h"
#include "asterfort/cjsinp.h"
#include "asterfort/cjsmat.h"
#include "asterfort/cjspla.h"
#include "asterfort/cjssmd.h"
#include "asterfort/cjssmi.h"
#include "asterfort/cjstde.h"
#include "asterfort/cjstel.h"
#include "asterfort/cjstid.h"
#include "asterfort/cjstis.h"
#include "asterfort/iunifi.h"
#include "asterfort/utmess.h"
#include "asterfort/get_varc.h"
    integer(kind=8) :: imat, ndt, ndi, nvi, iret
!
    real(kind=8) :: crit(*)
    real(kind=8) :: vind(*), vinf(*)
    real(kind=8) :: tempm, tempf, tref
    real(kind=8) :: epsd(6), deps(6), epsf(6)
    real(kind=8) :: sigd(6), sigf(6)
!
    real(kind=8) :: seuili, seuild, q0, rinit, pa, qinit
!
    real(kind=8) :: dsde(6, 6)
!
    real(kind=8) :: materf(14, 2)
!
!
    character(len=4) :: nivcjs
    character(len=6) :: mecand, mecanf
    character(len=7) :: etatd, etatf
    character(len=8) :: mod, typmod(*)
    character(len=16) :: opt
    integer(kind=8) :: niter, i, ndec
    real(kind=8) :: epscon
    real(kind=8) :: epsthe, epsthm
    real(kind=8) :: depsth(6), epsdth(6), alphaf, alpham
    aster_logical :: trac, lTemp
!
    real(kind=8) :: i1d
!
    integer(kind=8) :: umess
!
!       ----------------------------------------------------------------
    common/tdim/ndt, ndi
!       ----------------------------------------------------------------
!
!
!
! - Get temperatures
!
    call get_varc('RIGI', 1, 1, 'T', &
                  tempm, tempf, tref, l_temp_=lTemp)
!
    umess = iunifi('MESSAGE')
    mod = typmod(1)
    trac = .false.
!
! --    RECUPERATION COEF DE LA LOI CJS (INDEPENDANTS DE LA TEMPERATURE)
!                    NB DE CMP DIRECTES/CISAILLEMENT
!                    NB VARIABLES INTERNES
!                    NIVEAU DE LA LOI CJS: CJS1, CJS2 OU CJS3
    call cjsmat(mod, imat, tempf, materf, ndt, &
                ndi, nvi, nivcjs)
    pa = materf(12, 2)
    qinit = materf(13, 2)
!
!  COEF DE DILATATION LE MEME A TPLUS ET TMOINS
!
    alphaf = materf(3, 1)
    alpham = materf(3, 1)
!
    niter = 0
    ndec = 0
    epscon = 0
!
    i1d = 0.d0
    do i = 1, ndi
        i1d = i1d+sigd(i)
    end do
!
!     --  CALCUL DE DEPSTH ET EPSDTH
!     --------------------------------
!
    if (materf(3, 1) .eq. 0.d0) then
        epsthe = 0.d0
        epsthm = 0.d0
    else
        if (lTemp) then
            epsthe = alphaf*(tempf-tref)-alpham*(tempm-tref)
            epsthm = alpham*(tempm-tref)
        else
            call utmess('F', 'CALCULEL_15')
        end if
    end if
    depsth = 0.d0
    epsdth = 0.d0
    do i = 1, ndi
        depsth(i) = deps(i)-epsthe
        epsdth(i) = epsd(i)-epsthm
    end do
    do i = ndi+1, ndt
        depsth(i) = deps(i)
        epsdth(i) = epsd(i)
    end do
    if (ndt .lt. 6) then
        do i = ndt+1, 6
            depsth(i) = 0.d0
            epsdth(i) = 0.d0
        end do
    end if
!
!  TESTER QUE VIND DE NVI EST 1 2 OU 3
!
    if ((vind(nvi) .ne. 0.d0) .and. (vind(nvi) .ne. 1.d0) .and. (vind(nvi) .ne. 2.d0) .and. &
        (vind(nvi) .ne. 3.d0)) then
        write (umess, *) ' INDICATEUR DE PLASTICITE ERRONE : ', vind(nvi)
        call utmess('F', 'ALGORITH6_80')
    end if
!
!
! --    BLOCAGE DES VARIABLES INTERNES EN FONCTION DU NIVEAU DE LA LOI
!       CJS CHOISI, ET ON PREND DES VALEURS INITIALES DE Q ET R PETITES
!       MAIS NON NULLES
!
    if (nivcjs .eq. 'CJS1') then
        vind(1) = 1.d25*materf(12, 2)
        vind(2) = materf(2, 2)
    end if
!
    if (nivcjs .eq. 'CJS3') then
        vind(2) = materf(2, 2)
    end if
!
!  SI SEUIL ISOTROPE = 0 SEUIL ISOTROPE = PRESSION HYDRO
!
!
!
    if (vind(1) .eq. 0.d0) then
        if (i1d .lt. 0.d0) then
            q0 = (1.d-3*pa+i1d+qinit)/3.d0
        else
            q0 = (1.d-3*pa+qinit)/3.d0
        end if
!
        vind(1) = q0
        if (opt(1:14) .ne. 'RIGI_MECA_TANG') then
            vinf(1) = q0
        end if
    end if
!
! INITIALISATION SEUIL DEVIATOIRE SI NUL
!
    if (vind(2) .eq. 0.d0) then
        if (materf(14, 2) .eq. 0.d0) then
            rinit = 1.d-3
        else
            rinit = materf(14, 2)
        end if
        vind(2) = rinit
        if (opt(1:14) .ne. 'RIGI_MECA_TANG') then
            vinf(2) = rinit
        end if
    end if
!
!
! --    ETAT ELASTIC OU PLASTIC A T
!
    if (vind(nvi) .eq. 0.d0) then
        etatd = 'ELASTIC'
    end if
!
    if (vind(nvi) .eq. 1.d0) then
        etatd = 'PLASTIC'
        mecand = 'ISOTRO'
    end if
!
    if (vind(nvi) .eq. 2.d0) then
        etatd = 'PLASTIC'
        mecand = 'DEVIAT'
    end if
!
    if (vind(nvi) .eq. 3.d0) then
        etatd = 'PLASTIC'
        mecand = 'ISODEV'
    end if
!
!
!       ----------------------------------------------------------------
!       OPTIONS 'FULL_MECA' ET 'RAPH_MECA' = CALCUL DE SIG(T+DT)
!       ----------------------------------------------------------------
!
    if (opt .eq. 'RAPH_MECA' .or. opt .eq. 'FULL_MECA') then
!
!
!
!
! --    INTEGRATION ELASTIQUE SUR DT
!
!
!
        call cjsela(mod, crit, materf, depsth, sigd, &
                    sigf, nvi, vind, vinf, iret)
        if (iret .eq. 1) goto 999
!
!
! --    PREDICTION ETAT ELASTIQUE A T+DT :
!       FI(SIG(T+DT),VIN(T)) = 0 ?     SEUIL MECANISME ISOTROPE
!       FD(SIG(T+DT),VIN(T)) = 0 ?     SEUIL MECANISME DEVIATOIRE
!
!
        call cjssmi(materf, sigf, vind, seuili)
        call cjssmd(materf, sigf, vind, seuild)
!
!
!
        if ((seuili .gt. 0.d0) .or. (seuild .gt. 0.d0)) then
!
! ecriture des contraintes
!
!
! --      PREDICTION INCORRECTE > INTEGRATION ELASTO-PLASTIQUE SUR DT
!
            etatf = 'PLASTIC'
!
            call cjspla(mod, crit, materf, seuili, seuild, &
                        nvi, epsdth, depsth, sigd, vind, &
                        sigf, vinf, mecanf, nivcjs, niter, &
                        ndec, epscon, iret, trac)
            if (iret .eq. 1) goto 999
            if ((trac)) then
                etatf = 'ELASTIC'
            end if
!
        else
!
! --      PREDICTION CORRECTE > INTEGRATION ELASTIQUE FAITE
!
            etatf = 'ELASTIC'
        end if
!
    end if
!
!
!       ----------------------------------------------------------------
!       OPTIONS 'FULL_MECA' ET 'RIGI_MECA_TANG' = CALCUL DE DSDE
!       ----------------------------------------------------------------
!       CALCUL ELASTIQUE ET EVALUATION DE DSDE A (T)
!       POUR 'RIGI_MECA_TANG' ET POUR 'FULL_MECA'
!       ----------------------------------------------------------------
!
    if (opt .eq. 'RIGI_MECA_TANG') then
        dsde(:, :) = 0.d0
!
!
!
! REMARQUE: CALCUL DE DSDE A T AVEC MATERF CAR PARAMETRES CJS
!           INDEPENDANTS DE LA TEMPERATURE
!
!
!
! --      CALCUL MATRICE DE RIGIDITE ELASTIQUE
        if (etatd .eq. 'ELASTIC') then
            call cjstel(mod, materf, sigd, dsde)
        end if
!
! --      CALCUL MATRICE TANGENTE DU PROBLEME CONTINU
        if (etatd .eq. 'PLASTIC') then
!
!          MECANISME ISOTROPE SEUL
            if (mecand .eq. 'ISOTRO') then
                call cjstis(mod, materf, sigd, vind, dsde)
            end if
!
!          MECANISME DEVIATOIRE SEUL
            if ((mecand .eq. 'DEVIAT')) then
                call cjstde(mod, materf, nvi, epsdth, sigd, &
                            vind, dsde)
            end if
!
!          MECANISMES ISOTROPE ET DEVIATOIRE
            if (mecand .eq. 'ISODEV') then
                call cjstid(mod, materf, nvi, epsdth, sigd, &
                            vind, dsde)
            end if
!
!
        end if
!
    end if
!
!
    if (opt .eq. 'FULL_MECA') then
!
        dsde(:, :) = 0.d0
!
!
!
! --      CALCUL MATRICE DE RIGIDITE ELASTIQUE
        if (etatf .eq. 'ELASTIC') then
            call cjstel(mod, materf, sigf, dsde)
        end if
!
!
! --      CALCUL MATRICE TANGENTE DU PROBLEME CONTINU
        if (etatf .eq. 'PLASTIC') then
!
!          MECANISME ISOTROPE SEUL
            if (mecanf .eq. 'ISOTRO') then
                call cjstis(mod, materf, sigf, vinf, dsde)
            end if
!
!          MECANISME DEVIATOIRE SEUL
            if ((mecanf .eq. 'DEVIAT')) then
                epsf(1:ndt) = epsdth(1:ndt)+depsth(1:ndt)
                call cjstde(mod, materf, nvi, epsf, sigf, &
                            vinf, dsde)
            end if
!
!          MECANISMES ISOTROPE ET DEVIATOIRE
            if (mecanf .eq. 'ISODEV') then
                epsf(1:ndt) = epsdth(1:ndt)+depsth(1:ndt)
                call cjstid(mod, materf, nvi, epsf, sigf, &
                            vinf, dsde)
            end if
!
!
        end if
!
    end if
!
!  VARIABLES INTERNES POUR SORTIES
!
    if ((opt .eq. 'FULL_MECA') .or. (opt .eq. 'RAPH_MECA')) then
        call cjsinp(materf, epsdth, depsth, sigf, vinf, &
                    niter, nvi, nivcjs, ndec, epscon)
    end if
!
999 continue
end subroutine
