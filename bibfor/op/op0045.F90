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

subroutine op0045()
!        MODE_ITER_SIMULT
!        RECHERCHE DE MODES PAR ITERATION SIMULTANEE EN SOUS-ESPACE
!        (LANCZOS, JACOBI OU IRAM-ARPACK) OU METHODE DE TYPE QR (LAPACK)
!-----------------------------------------------------------------------
!        - POUR LE PROBLEME GENERALISE AUX VALEURS PROPRES :
!                         2
!                        L (M) Y  + (K) Y = 0
!
!          LES MATRICES (C) ET (M) SONT REELLES SYMETRIQUES
!          LA MATRICE (K) EST REELLE QCQ OU COMPLEXE SYMETRIQUE
!
!        - POUR LE PROBLEME QUADRATIQUE AUX VALEURS PROPRES :
!                         2
!                        L (M) Y  + L (C) Y + (K) Y = 0
!
!          LES MATRICES (C) ET (M) SONT REELLES SYMETRIQUES
!          LA MATRICE (K) EST REELLE OU COMPLEXE SYMETRIQUE
!
!          LES VALEURS PROPRES ET DES VECTEURS PROPRES SONT REELS OU
!          COMPLEXES CONJUGUEES OU NON
!-----------------------------------------------------------------------
! person_in_charge: olivier.boiteau at edf.fr
    implicit none
!
! --- INCLUDES DE MODE_ITER_SIMULT
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cresol.h"
#include "asterfort/infmaj.h"
#include "asterfort/vpcalj.h"
#include "asterfort/vpcalq.h"
#include "asterfort/vpcals.h"
#include "asterfort/vpcalt.h"
#include "asterfort/vpinis.h"
#include "asterfort/vpini0.h"
#include "asterfort/vpini1.h"
#include "asterfort/vpini2.h"
#include "asterfort/vpleci.h"
#include "asterfort/vpmpi.h"
#include "asterfort/vppost.h"
#include "asterfort/vpvers.h"
!
! --- DECLARATION DES VARIABLES LOCALES
!
    mpi_int :: mpicou, mpicow
    integer(kind=8) :: nbpari, nbparr, nbpark
    parameter(nbpari=8, nbparr=16, nbpark=3)
    integer(kind=8) :: iret, ibid, npivot, neqact, mxresf, nblagr, nstoc, nconv, nbvecg, nfreqg
    integer(kind=8) :: rangl, icom1, icom2
    real(kind=8) :: omemin, omemax, omeshi, vpinf, vpmax, rbid
    complex(kind=8) :: sigma
    character(len=4) :: mod45, mod45b
    character(len=8) :: modes, method
    character(len=16) :: typcon, compex
    character(len=19) :: matpsc, matopa, solveu, eigsol
    character(len=24) :: veclag, vecblo, vecrig, vecrer, vecrei, vecrek, vecvp, k24bid
    aster_logical :: flage, lcomod
!
!
! --------------------------------------------------------------------------------------------------
! --- ETAPE 0: INITS. ET LECTURE/VERIFICATION DES DONNEES
! --------------------------------------------------------------------------------------------------
    call infmaj()
    mod45 = 'OP45'
    mod45b = mod45
!
! --  ETAPE 0.0: INITIALISATIONS PROPRES A CET OPERATEUR
    call vpini0(compex, modes, typcon, solveu, eigsol, &
                matpsc, matopa, veclag, vecblo, vecrig, &
                vecrer, vecrei, vecrek, vecvp)
!
! --  ETAPE 0.1: PARALLELISME MULTI-NIVEAUX STEP 1
    call vpmpi(1, icom1_=icom1, icom2_=icom2, lcomod_=lcomod, &
               mpicou_=mpicou, mpicow_=mpicow, rangl_=rangl)
!
! --  ETAPE 0.2: LECTURE DES PARAMETRES SOLVEUR LINEAIRE ET CREATION DE LA SD SOLVEUR ASSOCIEE
    call cresol(solveu)
!
! --  ETAPE 0.3: LECTURE DES PARAMETRES SOLVEUR MODAL ET CREATION DE LA SD EIGENSOLVER ASSOCIEE
    call vpinis(eigsol)
!
! --  ETAPE 0.4: VERIFICATION DE LA COHERENCE DE LA SD EIGENSOLVER ET DES OBJETS SOUS-JACENTS
    call vpvers(eigsol, modes, .true._1)
!
! --------------------------------------------------------------------------------------------------
! --- ETAPES 1: PREPARATION DU CALCUL MODAL
! --------------------------------------------------------------------------------------------------
! --  ETAPE 1.1: TRAITEMENTS NUMERIQUES (SOLVEUR LINEAIRE, LAGRANGE, MODES RIGIDES,
! --             BORNES DE TRAVAIL EFFECTIVES, CALCUL DU NOMBRE DE MODES, FACTO. MATRICE SHIFTEE
! --             DETERMINATION TAILLE DE L'ESPACE DE PROJECTION)
    call vpini1(eigsol, modes, solveu, typcon, vecblo, &
                veclag, vecrig, matpsc, matopa, iret, &
                nblagr, neqact, npivot, nstoc, omemax, &
                omemin, omeshi, sigma, mod45)
    if (iret .ne. 0) goto 999
!
! --  ETAPE 1.2: PARALLELISME MULTI-NIVEAUX STEP 2
    call vpmpi(2, eigsol, lcomod_=lcomod, mpicou_=mpicou, mpicow_=mpicow, &
               nbvecg_=nbvecg, nfreqg_=nfreqg, rangl_=rangl)
!
! --  ETAPE 1.3: CREATION ET INITIALISATION DES SDS RESULTATS
    call vpini2(eigsol, lcomod, nbvecg, nfreqg, nbpark, &
                nbpari, nbparr, vecrer, vecrei, vecrek, &
                vecvp, mxresf)
!
! --------------------------------------------------------------------------------------------------
! --- ETAPE 2: CALCUL MODAL PROPREMENT DIT SUIVANT LA METHODE CHOISIE
! --------------------------------------------------------------------------------------------------
    call vpleci(eigsol, 'K', 6, k24bid, rbid, ibid)
    method = ''
    method = trim(k24bid)
    select case (method)
    case ('SORENSEN')
        call vpcals(eigsol, vecrer, vecrei, vecrek, vecvp, &
                    matopa, mxresf, neqact, nblagr, omemax, &
                    omemin, omeshi, solveu, vecblo, veclag, &
                    sigma, npivot, flage, nconv, vpinf, vpmax, mod45b, &
                    k24bid, k24bid, ibid, K24bid, ibid, rbid)
    case ('TRI_DIAG')
        call vpcalt(eigsol, vecrer, vecrei, vecrek, vecvp, &
                    matopa, matpsc, mxresf, nblagr, nstoc, &
                    omemax, omemin, omeshi, solveu, vecblo, &
                    veclag, vecrig, sigma, npivot, flage, &
                    nconv, vpinf, vpmax)
    case ('JACOBI')
        call vpcalj(eigsol, vecrer, vecrei, vecrek, vecvp, &
                    matopa, matpsc, mxresf, nblagr, omemax, &
                    omemin, omeshi, solveu, vecblo, npivot, &
                    flage, nconv, vpinf, vpmax)
    case ('QZ')
        call vpcalq(eigsol, vecrer, vecrei, vecrek, vecvp, &
                    mxresf, neqact, nblagr, omemax, omemin, &
                    omeshi, vecblo, sigma, npivot, flage, &
                    nconv, vpinf, vpmax)
    case default
        ASSERT(.false.)
    end select
!
! --------------------------------------------------------------------------------------------------
! --- ETAPE 3: POST-TRAITEMENTS/VERIFICATIONS GLOBAUX A TOUTES LES METHODES/SITUATIONS
! ---         + PARALLELISME MULTI-NIVEAUX STEP 3 ET 4
! ---         + NETTOYAGE EXPLICITE DES OBJETS JEVEUX GLOBAUX A L'OPERATEUR (BASE VOLATILE)
! --------------------------------------------------------------------------------------------------
    call vppost(vecrer, vecrei, vecrek, vecvp, nbpark, &
                nbpari, nbparr, mxresf, nconv, nblagr, &
                nfreqg, modes, typcon, compex, eigsol, &
                matopa, matpsc, solveu, vecblo, veclag, &
                flage, icom1, icom2, mpicou, mpicow, &
                omemax, omemin, vpinf, vpmax, lcomod, mod45)
999 continue
!
! --------------------------------------------------------------------------------------------------
! --- ETAPE 4: NETTOYAGE DES COMMUNICATEURS, PARALLELISME MULTI-NIVEAUX STEP 5
! --------------------------------------------------------------------------------------------------
    call vpmpi(5, lcomod_=lcomod, &
               mpicou_=mpicou, mpicow_=mpicow)
!
!
!     FIN DE OP0045
!
end subroutine
