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

subroutine vpsorc(lmasse, ldynfa, nbeq, nbvect, nfreq, &
                  tolsor, vect, resid, workd, workl, &
                  lonwl, selec, dsor, sigma, vaux, &
                  workv, rwork, ddlexc, ddllag, neqact, &
                  maxitr, ifm, niv, priram, alpha, &
                  nconv, flage, solveu)
!
!     SUBROUTINE ASTER ORCHESTRANT LA METHODE DE SORENSEN: UN ARNOLDI
!     AVEC REDEMARRAGE IMPLICITE VIA QR (VERSION ARPACK 2.8).
!---------------------------------------------------------------------
!     PARTANT DU PROBLEME GENERALISE AUX VALEURS PROPRES
!                           (A)*X = LAMBDA*(B)*X
!     AVEC
!      - LES MATRICES SYMETRIQUES (A) ET (B), CORRESPONDANT AUX
!     MATRICES DE RAIDEUR ET A CELLE DE MASSE (RESP. RAIDEUR GEOMETRIQUE
!     , EN FLAMBEMENT),
!      - LE COMPLEXE, LAMBDA, CORRESPONDANT AU CARRE DE LA PULSATION
!       (RESP. L'OPPOSE DE LA CHARGE CRITIQUE, EN FLAMBEMENT),
!      - LE OU LES VECTEURS PROPRES ASSOCIES, X, (A- ET B- ORTHOGO
!        NAUX ENTRE EUX AINSI QU'AVEC CEUX DES AUTRES VALEURS PROPRES).
!
!     ON RESOUT LE PROBLEME STANDARD
!                             (OP)*X =  MU*X
!     AVEC
!       - L'OPERATEUR SHIFTE (OP) = INV((A)-SIGMA*(B))*(B),
!       - LE 'SHIFT' COMPLEXE SIGMA,
!       - LA VALEUR PROPRE MU = 1/(LAMBDA-SIGMA),
!       - LE OU LES MEMES VECTEURS PROPRES QU'INITIALEMENT, X.
!   ------------------------------------------------------------------
!   CETTE METHODE, PARTANT D'UNE FACTORISATION DE TYPE ARNOLDI D'ORDRE
!   M=K+P DU PROBLEME, PILOTE UN RESTART A L'ORDRE K SUR P NOUVELLES
!   ITERATIONS. CE RESTART PERMET D'AMELIORER LES K PREMIERES VALEURS
!   PROPRES SOUHAITEES, LES P DERNIERES SERVANT UNIQUEMENT AUX CALCULS
!   AUXILIAIRES.
!   ELLE PERMET
!     - DE MINIMISER LA TAILLE DU SOUS-ESPACE DE PROJECTION,
!     - D'EFFECTUER DES RESTARTS DE MANIERE TRANSPARENTE, EFFICACE ET
!       AVEC DES PRE-REQUIS MEMOIRE FIXES,
!     - DE MIEUX PRENDRE EN COMPTE LES MULTIPLICITES,
!     - DE TRAITER AVEC UN BON COMPROMIS LA STRATEGIE DE RE-ORTHONORMA
!       LISATION.
!   ------------------------------------------------------------------
!     PARAMETRES D'APPELS:
!
! IN  LMASSE : IS : DESCRIPTEUR MATRICE DE "MASSE".
! IN  LDYNFA : IS : DESCRIPTEUR MATRICE DE "RAIDEUR"-SHIFT"MASSE"
!                     FACTORISEE.
! IN  NBEQ   : IS : DIMENSION DES VECTEURS.
! IN  NBVECT : IS : DIMENSION DE L'ESPACE DE PROJECTION.
! IN  NFREQ  : IS : NOMBRE DE VALEURS PROPRES DEMANDEES.
! IN  TOLSOR : R8 : NORME D'ERREUR SOUHAITEE (SI 0.D0 ALORS LA VALEUR
!                   PAR DEFAUT EST LA PRECISION MACHINE).
! IN  LONWL  : IS : TAILLE DU VECTEUR DE TRAVAIL WORKL.
! IN  FSHIFT : R8 : VALEUR DU SHIFT SIGMA EN OMEGA2.
! IN  DDLEXC : IS : DDLEXC(1..NBEQ) VECTEUR POSITION DES DDLS BLOQUES.
! IN  DDLLAG : IS : DDLLAG(1..NBEQ) VECTEUR POSITION DES LAGRANGES.
! IN  NEQACT : IS : NOMBRE DE DDLS ACTIFS.
! IN  MAXITR : IS : NOMBRE MAXIMUM DE RESTARTS.
! IN  IFM    : IS : UNITE LOGIQUE D'IMPRESSION DES .MESS
! IN  NIV    : IS : NIVEAU D'IMPRESSION
! IN  PRIRAM : IS : PRIRAM(1..8) VECTEUR NIVEAU D'IMPRESSION ARPACK
! IN  ALPHA  : R8 : PARAMETRE VPGSKP D'ORTHONORMALISATION.
! IN  OMECOR : R8 : OMEGA2 DE CORPS RIGIDE
!
! OUT VECT   : C16 : VECT(1..NBEQ,1..NBVECT) MATRICE DES
!                   VECTEURS D'ARNOLDI.
! OUT RESID  : C16 : RESID(1..NBEQ) VECTEUR RESIDU.
! OUT WORKD  : C16 : WORKD(1..3*NBEQ) VECTEUR DE TRAVAIL PRIMAIRE IRAM
! OUT WORKL  : C16 : WORKL(1..LONWL) VECTEUR DE TRAVAIL SECONDAIRE IRAM
! OUT RWORK  : R8  : RWORK(1..NBEQ) VECTEUR DE TRAVAIL REEL
! OUT SELEC  : LS : SELEC(1..NBVECT) VECTEUR DE TRAVAIL POUR DNEUPD.
! OUT DSOR   : C16 : DSOR(1..NFREQ+1) VECTEUR DES VALEURS PROPRES.
! OUT VAUX   : C16 : VAUX(1..NBEQ) VECTEUR DE TRAVAIL.
! OUT WORKV  : C16: WORKV(1..3*NBVECT) VECTEUR DE TRAVAIL POUR DNEUPD
!                     ET VPGSKP.
! OUT NCONV  : IS : NOMBRE DE MODES CONVERGES.
! OUT FLAGE  : LO : FLAG PERMETTANT DE GERER LES IMPRESSIONS F
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
! aslint: disable=W1504
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/mcmult.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "asterfort/vpordc.h"
#include "asterfort/znaupd.h"
#include "asterfort/zneupd.h"
!
    integer(kind=8) :: lmasse, ldynfa, nbeq, nbvect, nfreq, lonwl
    integer(kind=8) :: ddlexc(nbeq), ddllag(nbeq), neqact
    integer(kind=8) :: maxitr, ifm, niv, priram(8), nconv
    real(kind=8) :: tolsor, alpha, rwork(*)
    complex(kind=8) :: vect(nbeq, *), resid(*), workd(*), workl(*), dsor(*), vaux(*), workv(*)
    complex(kind=8) :: sigma
    aster_logical :: selec(nbvect), flage
    character(len=19) :: solveu
!
!
! DECLARATION VARIABLES LOCALES
!
! POUR LE FONCTIONNEMENT GLOBAL
    integer(kind=8) :: i, j
    integer(kind=8) :: iret
!
! POUR ARPACK
    integer(kind=8) :: ido, info, ishfts, mode, iparam(11), ipntr(14), vali(6)
    aster_logical :: rvec
    character(len=1) :: bmat, kbid
    character(len=2) :: which
    character(len=19) :: k19bid, matass, chcine, criter
!
    integer(kind=8) :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr
    integer(kind=8) :: mneigh, mnapps, mngets, mneupd
    common/debug/&
     &  logfil, ndigit, mgetv0,&
     &  mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
!------------------------------------------------------------------
! INITIALISATION POUR ARPACK
!
! NIVEAU D'IMPRESSION ARPACK
    ndigit = -3
    logfil = ifm
    mgetv0 = priram(1)
    mnaupd = priram(2)
    mnaup2 = priram(3)
    mnaitr = priram(4)
    mneigh = priram(5)
    mnapps = priram(6)
    mngets = priram(7)
    mneupd = priram(8)
!
! FONCTIONNEMENT D'ARPACK
    ido = 0
    info = 0
    ishfts = 1
    mode = 3
    rvec = .true.
    bmat = 'G'
    which = 'LM'
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(5) = 0
    iparam(4) = 1
    iparam(7) = mode
!
! INIT. OBJETS ASTER
    matass = zk24(zi(ldynfa+1)) (1:19)
    chcine = ' '
    criter = ' '
    k19bid = ' '
!------------------------------------------------------------------
! BOUCLE PRINCIPALE
!
20  continue
!
! CALCUL DES VALEURS PROPRES DE (OP)
    call znaupd(ido, bmat, nbeq, which, nfreq, &
                tolsor, resid, nbvect, vect, nbeq, &
                iparam, ipntr, workd, workl, lonwl, &
                rwork, info, neqact, alpha)
!
! NOMBRE DE MODES CONVERGES
    nconv = iparam(5)
!
! GESTION DES FLAGS D'ERREURS
    if ((info .eq. 1) .and. (niv .ge. 1)) then
        vali(1) = maxitr
        call utmess('I', 'ALGELINE6_89', si=vali(1))
    else if (info .eq. 2) then
        call utmess('F', 'ALGELINE3_72')
    else if ((info .eq. 3) .and. (niv .ge. 1)) then
        call utmess('I', 'ALGELINE6_90')
    else if (info .eq. -7) then
        call utmess('F', 'ALGELINE3_73')
    else if (info .eq. -8) then
        call utmess('F', 'ALGELINE3_74')
    else if (info .eq. -9) then
        call utmess('F', 'ALGELINE3_75')
    else if ((info .eq. -9999) .and. (niv .ge. 1)) then
        call utmess('F', 'ALGELINE6_91')
    else if (info .lt. 0) then
        vali(1) = info
        call utmess('F', 'ALGELINE4_82', si=vali(1))
    end if
!
!---------------------------------------------------------------------
! ZONE GERANT LA 'REVERSE COMMUNICATION' VIA IDO
!
    if (ido .eq. -1) then
! CALCUL DU Y = (OP)*X INITIAL
! 1/ CALCUL D'UN ELT. INITIAL X REPONDANT AU C.I. DE LAGRANGE
! 2/ CALCUL DE Y = (OP)* X AVEC DDL CINEMATIQUEMENT BLOQUES
! X <- X*DDL_LAGRANGE
        do j = 1, nbeq
            vaux(j) = workd(ipntr(1)+j-1)*ddllag(j)
        end do
! X <- (INV((A)-SIGMA*(B))*X)*DDL_LAGRANGE
        call resoud(matass, k19bid, solveu, chcine, 1, &
                    k19bid, k19bid, kbid, [0.d0], vaux, &
                    criter, .false._1, 0, iret)
        do j = 1, nbeq
            workd(ipntr(1)+j-1) = vaux(j)*ddllag(j)
        end do
! X <- (OP)*(X*DDL_BLOQUE)
        call mcmult('ZERO', lmasse, workd(ipntr(1)), vaux, 1, &
                    .false._1)
        do j = 1, nbeq
            vaux(j) = workd(ipntr(1)+j-1)*ddlexc(j)
        end do
        call resoud(matass, k19bid, solveu, chcine, 1, &
                    k19bid, k19bid, kbid, [0.d0], vaux, &
                    criter, .false._1, 0, iret)
! RETOUR VERS DNAUPD
        do j = 1, nbeq
            workd(ipntr(2)+j-1) = vaux(j)*ddlexc(j)
        end do
        goto 20
!
    else if (ido .eq. 1) then
! CALCUL DU Y = (OP)*X CONNAISSANT DEJA (B)*X (EN FAIT ON CONNAIT
! SEULEMENT (ID)*X VIA IDO= 2 CAR PRODUIT SCALAIRE= L2)
! X <- (B)*X*DDL_BLOQUE
        call mcmult('ZERO', lmasse, workd(ipntr(3)), vaux, 1, &
                    .false._1)
        do j = 1, nbeq
            vaux(j) = vaux(j)*ddlexc(j)
        end do
! X <- (OP)*X
        call resoud(matass, k19bid, solveu, chcine, 1, &
                    k19bid, k19bid, kbid, [0.d0], vaux, &
                    criter, .false._1, 0, iret)
! RETOUR VERS DNAUPD
        do j = 1, nbeq
            workd(ipntr(2)+j-1) = vaux(j)*ddlexc(j)
        end do
        goto 20
!
    else if (ido .eq. 2) then
! X <- X*DDL_BLOQUE  (PRODUIT SCALAIRE= L2)
        do j = 1, nbeq
            workd(ipntr(2)+j-1) = workd(ipntr(1)+j-1)*ddlexc(j)
        end do
! RETOUR VERS DNAUPD
        goto 20
!
! GESTION DES MODES CONVERGES
    else if (ido .eq. 99) then
        if (nconv .lt. nfreq) then
            vali(1) = nconv
            vali(2) = nfreq
            call utmess('E', 'ALGELINE5_49', ni=2, vali=vali)
            flage = .true.
        else if (nconv .gt. nfreq) then
            vali(1) = nconv
            vali(2) = nfreq
            call utmess('I', 'ALGELINE5_50', ni=2, vali=vali)
            nconv = nfreq
        end if
    else
        ASSERT(.false.)
    end if
!--------------------------------------------------------------------
! CALCUL DES MODES PROPRES APPROCHES DU PB INITIAL
!
    info = 0
    call zneupd(rvec, 'A', selec, dsor, vect, &
                nbeq, sigma, workv, bmat, nbeq, &
                which, nfreq, tolsor, resid, nbvect, &
                vect, nbeq, iparam, ipntr, workd, &
                workl, lonwl, rwork, info)
!
! GESTION DES FLAGS D'ERREURS
    if (info .eq. 1) then
        call utmess('F', 'ALGELINE3_74')
    else if (info .eq. -7) then
        call utmess('F', 'ALGELINE3_73')
    else if (info .eq. -8) then
        call utmess('F', 'ALGELINE3_76')
    else if (info .eq. -9) then
        call utmess('F', 'ALGELINE3_77')
    else if (info .eq. -14) then
        call utmess('F', 'ALGELINE3_78')
    else if (info .lt. 0) then
        vali(1) = info
        call utmess('F', 'ALGELINE4_82', si=vali(1))
    end if
!--------------------------------------------------------------------
! TESTS ET POST-TRAITEMENTS
!
! POUR TEST
!      DO 59 J=1,NCONV
!       RBID=DIMAG(DSOR(J))/DBLE(DSOR(J))
!       WRITE(6,*) '******** VALEUR DE RITZ N ********',J
!       WRITE(6,*) 'RE: LANDAJ/ FJ INIT',DSOR(J)
!   59 CONTINUE
!
! TRI DES MODES PROPRES PAR RAPPORT AU NCONV DSOR(I)
    call vpordc(1, 0, nconv, dsor, vect, &
                nbeq)
!
    do i = 1, nconv
        dsor(i) = dsor(i)+sigma
    end do
!
! FIN DE VPSORC
!
end subroutine
