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
subroutine wpsorn(appr, lmasse, lamor, lmatra, nbeq, &
                  nbvect, nfreq, tolsor, vect, resid, &
                  workd, workl, lonwl, selec, dsor, &
                  vpr, vpi, sigma, vaux, workv, &
                  ddlexc, ddllag, neqact, maxitr, ifm, &
                  niv, priram, alpha, nconv, flage, &
                  vaur, vauc, vaul, solveu)
!
!---------------------------------------------------------------------
!
!     SUBROUTINE ASTER ORCHESTRANT LA METHODE DE SORENSEN: UN ARNOLDI
!     COMPLEXE AVEC REDEMARRAGE IMPLICITE VIA QR (VERSION ARPACK 2.8).
!     CAS QUADRATIQUE REEL
!---------------------------------------------------------------------
!     PARTANT DU PROBLEME QUADRATIQUE AUX VALEURS PROPRES
!      REDUCTION A UN PROBLEME GENERALISE
!
!         !K   0! !P!          !-C   -M! !P!
!         !     ! ! ! = LAMBDA !       ! ! ! <=> K.G*Z = LAMBDA*M.G*Z
!         !0  -M! !Q!          !-M    0! !Q!
!
!     AVEC
!      - LES MATRICES SYMETRIQUES (A) ET (B), CORRESPONDANT AUX
!     MATRICES DE RAIDEUR ET A CELLE DE MASSE (RESP. RAIDEUR GEOMETRIQUE
!     , EN FLAMBEMENT),
!      - LE COMPLEXE, LAMBDA, CORRESPONDANT AU CARRE DE LA PULSATION
!      - LE OU LES VECTEURS PROPRES COMPLEXES ASSOCIES, X,
!        (A- ET B- ORTHOGONAUX ENTRE EUX AINSI QU'AVEC CEUX DES AUTRES
!        VALEURS PROPRES).
!
!  ------------------------------------------------------------------
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
! IN  LMATRA : IS : DESCRIPTEUR MATRICE DE "RAIDEUR"-SHIFT"MASSE"
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
! OUT VECT   : R8 : VECT(1..NBEQ,1..NBVECT) MATRICE DES
!                   VECTEURS D'ARNOLDI.
! OUT RESID  : R8 : RESID(1..NBEQ) VECTEUR RESIDU.
! OUT WORKD  : R8 : WORKD(1..3*NBEQ) VECTEUR DE TRAVAIL PRIMAIRE IRAM
! OUT WORKL  : R8 : WORKL(1..LONWL) VECTEUR DE TRAVAIL SECONDAIRE IRAM
! OUT SELEC  : LS : SELEC(1..NBVECT) VECTEUR DE TRAVAIL POUR DNEUPD.
! OUT DSOR   : R8 : DSOR(1..NFREQ+1,1..2) MATRICE DES VALEURS PROPRES.
! OUT VAUX   : R8 : VAUX(1..NBEQ) VECTEUR DE TRAVAIL POUR WPSORN.
! OUT WORKV  : R8 : WORKV(1..3*NBVECT) VECTEUR DE TRAVAIL POUR DNEUPD
!                     ET VPGSKP.
! OUT NCONV  : IS : NOMBRE DE MODES CONVERGES.
! OUT FLAGE  : LO : FLAG PERMETTANT DE GERER LES IMPRESSIONS
! IN  SOLVEU : K19 : SD SOLVEUR POUR PARAMETRER LE SOLVEUR LINEAIRE
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
! aslint: disable=W1504
    implicit none
!
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dnaupd.h"
#include "asterfort/dneupd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/wp2ay1.h"
    character(len=1) :: appr
    integer(kind=8) :: lmasse, lmatra, nbeq, nbvect, nfreq, lonwl, ddlexc(*)
    integer(kind=8) :: ddllag(*), neqact, maxitr
    integer(kind=8) :: ifm, niv, priram(8), nconv, lamor
    real(kind=8) :: tolsor, resid(*), workd(*), workl(*), vaux(*), vaur(2*nbeq, *), vpr(*), vpi(*)
    real(kind=8) :: workv(*), alpha, dsor(nfreq+1, *)
    aster_logical :: selec(*), flage
    complex(kind=8) :: sigma, vect(nbeq, *), vauc(2*nbeq, *), vaul(2*nbeq, *)
    character(len=19) :: solveu
!--------------------------------------------------------------------
! DECLARATION VARIABLES LOCALES
!
! POUR LE FONCTIONNEMENT GLOBAL
    integer(kind=8) :: i, j
    integer(kind=8) :: vali(11)
    integer(kind=8) :: au1, au2, au3, av
!
! POUR ARPACK
    integer(kind=8) :: ido, info, ishfts, mode, iparam(11), ipntr(14)
    real(kind=8) :: sigmar, sigmai
    aster_logical :: rvec
    character(len=1) :: bmat
    character(len=2) :: which
!
    integer(kind=8) :: logfil, ndigit, mgetv0, mnaupd, mnaup2, mnaitr
    integer(kind=8) :: mneigh, mnapps, mngets, mneupd
    common/debug/&
     &  logfil, ndigit, mgetv0,&
     &  mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd
!------------------------------------------------------------------
    do i = 1, 11
        iparam(i) = 0
    end do
!
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
    if (appr .eq. 'R') then
        mode = 3
    else
        mode = 4
    end if
    sigmar = dble(sigma)
    sigmai = dimag(sigma)
    rvec = .true.
    bmat = 'G'
    which = 'LM'
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(4) = 1
    iparam(7) = mode
!------------------------------------------------------------------
! BOUCLE PRINCIPALE
!
    call jemarq()
!     ---- ALLOCATION DES ZONES DE TRAVAIL ---
    call wkvect('&&WPSORN.VECTEUR.AUX.U1R', 'V V R', nbeq, au1)
    call wkvect('&&WPSORN.VECTEUR.AUX.U2R', 'V V R', nbeq, au2)
    call wkvect('&&WPSORN.VECTEUR.AUX.U3R', 'V V R', nbeq, au3)
    call wkvect('&&WPSORN.VECTEUR.AUX.VC ', 'V V C', nbeq, av)
!******************************************************************
20  continue
!
! CALCUL DES VALEURS PROPRES DE (OP)
    call dnaupd(ido, bmat, 2*nbeq, which, nfreq, &
                tolsor, resid, nbvect, vaur, 2*nbeq, &
                iparam, ipntr, workd, workl, lonwl, &
                info, 2*neqact, alpha)
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
            vaux(j) = 0.d0*workd(ipntr(1)+j-1)*ddllag(j)
            vaux(j+nbeq) = workd(ipntr(1)+nbeq+j-1)*ddllag(j)
        end do
        call wp2ay1(appr, lmatra, lmasse, lamor, sigma, &
                    ddlexc, vaux(1), vaux(nbeq+1), workd(ipntr(1)), workd(ipntr(1)+nbeq), &
                    zr(au1), zr(au2), zr(au3), zc(av), nbeq, &
                    solveu)
        do j = 1, nbeq
            vaux(j) = workd(ipntr(1)+j-1)*ddlexc(j)
            vaux(j+nbeq) = workd(ipntr(1)+nbeq+j-1)*ddlexc(j)
        end do
        call wp2ay1(appr, lmatra, lmasse, lamor, sigma, &
                    ddlexc, vaux(1), vaux(nbeq+1), workd(ipntr(1)), workd(ipntr(1)+nbeq), &
                    zr(au1), zr(au2), zr(au3), zc(av), nbeq, &
                    solveu)
! RETOUR VERS DNAUPD
        do j = 1, nbeq
            workd(ipntr(2)+j-1) = workd(ipntr(1)+j-1)*ddlexc(j)
            workd(ipntr(2)+nbeq+j-1) = workd(ipntr(1)+nbeq+j-1)*ddlexc(j)
        end do
        goto 20
!
    else if (ido .eq. 1) then
! CALCUL DU Y = (OP)*X CONNAISSANT DEJA (B)*X (EN FAIT ON CONNAIT
! SEULEMENT (ID)*X VIA IDO= 2 CAR PRODUIT SCALAIRE= L2)
! X <- (OP)*X
        do j = 1, nbeq
            workd(ipntr(3)+j-1) = workd(ipntr(3)+j-1)*ddlexc(j)
            workd(ipntr(3)+nbeq+j-1) = workd(ipntr(3)+nbeq+j-1)*ddlexc(j)
        end do
        call wp2ay1(appr, lmatra, lmasse, lamor, sigma, &
                    ddlexc, workd(ipntr(3)), workd(ipntr(3)+nbeq), vaux(1), vaux(nbeq+1), &
                    zr(au1), zr(au2), zr(au3), zc(av), nbeq, &
                    solveu)
! RETOUR VERS DNAUPD
        do j = 1, nbeq
            workd(ipntr(2)+j-1) = vaux(j)*ddlexc(j)
            workd(ipntr(2)+nbeq+j-1) = vaux(j+nbeq)*ddlexc(j)
        end do
        goto 20
!
    else if (ido .eq. 2) then
! X <- X*DDL_BLOQUE  (PRODUIT SCALAIRE= L2)
        do j = 1, nbeq
            workd(ipntr(2)+j-1) = workd(ipntr(1)+j-1)*ddlexc(j)
            workd(ipntr(2)+nbeq+j-1) = workd(ipntr(1)+nbeq+j-1)*ddlexc( &
                                       j)
        end do
! RETOUR VERS DNAUPD
        goto 20
!
! GESTION DES MODES CONVERGES
    else if (ido .eq. 99) then
        if (nconv .lt. nfreq) then
            vali(1) = nconv
            vali(2) = nfreq
            call utmess('A', 'ALGELINE5_49', ni=2, vali=vali)
            flage = .false.
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
    call dneupd(rvec, 'A', selec, dsor, dsor(1, 2), &
                vaur, 2*nbeq, sigmar, sigmai, workv, &
                bmat, 2*nbeq, which, nfreq, tolsor, &
                resid, nbvect, vaur, 2*nbeq, iparam, &
                ipntr, workd, workl, lonwl, info)
!
! GESTION DES FLAGS D'ERREURS
    if (info .eq. 1) then
        call utmess('F', 'ALGELINE3_74')
    else if (info .eq. -7) then
        call utmess('F', 'ALGELINE3_73')
    else if (info .eq. -8) then
        call utmess('F', 'ALGELINE3_98')
    else if (info .eq. -9) then
        call utmess('F', 'ALGELINE3_99')
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
!       WRITE(IFM,*) '******** VALEUR DE RITZ N ********',J
!       WRITE(IFM,*) 'RE: LANDAJ/ FJ INIT',DSOR(J,1)
!       WRITE(IFM,*) 'IM: LANDAJ/ FJ INIT',DSOR(J,2)
!  59  CONTINUE
!
    do j = 1, nconv
        vpr(j) = dsor(j, 1)
        vpi(j) = dsor(j, 2)
    end do
!
!      REMPLISSAGE DE VAUC AVEC VAUR
!     DNEUPD SORT LES VECTEURS COMPLEXES PAR COLONNE
!     PARTIE REELLE PUIS PARTIE IMAGINAIRE SANS LES CONJUGUES
!
!
    do j = 2, nconv, 2
        do i = 1, 2*nbeq
            vaul(i, j/2) = dcmplx(vaur(i, j-1), vaur(i, j))
        end do
    end do
!
!
    do j = 2, nconv, 2
        do i = 1, 2*nbeq
            vauc(i, j-1) = vaul(i, j/2)
            vauc(i, j) = dconjg(vaul(i, j/2))
        end do
    end do
!*****************************************************************
!
    do j = 1, nconv
        do i = 1, nbeq
!     --- REMPLISSAGE DU VECT PAR LA PARTIE BASSE DE VECTA
            vect(i, j) = vauc(i+nbeq, j)
        end do
    end do
!
!
!     --- DESTRUCTION DES OJB TEMPORAIRES
    call jedetr('&&WPSORN.VECTEUR.AUX.U1R')
    call jedetr('&&WPSORN.VECTEUR.AUX.U2R')
    call jedetr('&&WPSORN.VECTEUR.AUX.U3R')
    call jedetr('&&WPSORN.VECTEUR.AUX.U5C')
    call jedetr('&&WPSORN.VECTEUR.AUX.VC ')
!
    call jedema()
!
! FIN DE WPSORN
!
end subroutine
