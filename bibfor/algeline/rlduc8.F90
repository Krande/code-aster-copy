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
subroutine rlduc8(nommat, hcol, adia, ablo, neq, &
                  nbbloc, xsol, nbsol)
    implicit none
!
! aslint: disable=W1306
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: nommat
    integer(kind=8) :: neq
    integer(kind=8) :: hcol(*), adia(*), ablo(*)
    complex(kind=8) :: xsol(neq, *)
!
! BUT : RESOLUTION DU SYSTEME A COEFFICIENTS COMPLEXES :  A * X = B
!       ELLE CORRESPOND AU SOLVEUR 'LDLT' + COMPLEXE + NON SYM
!       (CETTE ROUTINE EST LA SOEUR JUMELLE DE RLDUR8)
! ======================================================================
!     RESOLUTION DU SYSTEME A COEFFICIENTS COMPLEXES :  A * X = B
!     LA MATRICE EST NON SYMETRIQUE ET A ETE FACTORISEE SOUS FORME L*D*U
!     LA RESOLUTION EST EN PLACE
!
!     ON PEUT RESOUDRE SUR UNE SOUS-MATRICE DE A :
!     ON PREND LES NEQ PREMIERES LIGNES ET COLONNES (NEQ PEUT ETRE
!     INFERIEUR A LA DIMENSION DE LA MATRICE).
!
!     ON PEUT RESOUDRE NBSOL SYSTEMES D'UN COUP A CONDITION
!     QUE LES VECETURS SOIENT CONSECUTIFS EN MEMOIRE :
!     XSOL EST UN VECTEUR DE NBSOL*NEQ REELS
!     ------------------------------------------------------------------
!
! IN  NOMMAT  :    : NOM UTILISATEUR DE LA MATRICE A FACTORISER
! IN  HCOL    : IS : HCOL DE LA MATRICE
!             HCOL(I) RENVOIE LA HAUTEUR DE LA I-EME COLONNE
! IN  ADIA    : IS : ADRESSE DU TERME DIAGONALE DANS SON BLOC
!             ADIA(I) RENVOIE L'ADRESSE DE LA I-EME LIGNE DANS SON BLOC
! IN  ABLO    :  :   POINTEUR DE BLOC
!             ABLO(I+1) RENVOIE LE NO DE LA DERNIERE LIGNE DU I-EME BLOC
!
! IN  NEQ     : IS : NOMBRE D'EQUATIONS PRISES EN COMPTE
!                    C'EST AUSSI LA DIMENSION DES VECTEURS XSOL.
! IN  NBBLOC  : IS : NOMBRE DE BLOC DE LA MATRICE
! VAR XSOL    : C8 : EN ENTREE LES SECONDS MEMBRES
!                    EN ENTREE LES SOLUTIONS
! IN  NBSOL   : IS : NOMBRE DE SOLUTIONS / SECONDS MEMBRES
!     ------------------------------------------------------------------
!  LES TABLEAUX RELATIFS AUX BLOCS ONT ETE CONSTRUITS A PARTIR
!  DES NBBLOCS PREMIERS BLOCS (I.E LES BLOCS SUP)
!  MAIS ILS RESTENT VALABLES POUR LES NBBLOCS BLOCS SUIVANTS
!   (I.E LES BLOCS INF) PUISQUE LES POFILS SUP ET INF SONT IDENTIQUES
!     ------------------------------------------------------------------
!
!     --- RAPPEL SOMMAIRE DE L'ALGORITHME ------------------------------
!
!     POUR I = 1,2, ... ,N
!     !  ACC = 0
!     !  POUR K = 1,2, ... ,I-1
!     !  !  ACC = ACC +  K(I,K)* X(K)
!     !  FIN_POUR
!     !  X(I) = X(I)-ACC/K(I,I)
!     FIN_POUR
!
!     POUR I = N,N-1,... ,1
!     !  ACC = 0
!     !  POUR K = I+1, ... ,N
!     !  !  ACC = ACC +  K(K,I)* X(K)
!     !  FIN_POUR
!     !  X(I) = X(I) - ACC
!     FIN_POUR
!     ------------------------------------------------------------------
!
!     REFERENCE (HISTORIQUE) :
!     (1) P.D. CROUT,
!         A SHORT METHOD FOR EVALUATING DETERMINANTS AND SOLVING SYSTEMS
!         OF LINEAR EQUATIONS WITH REAL OR COMPLEX COEFFICIENTS.
!         AIEE TRANSACTION VOL 60, PP 1235-1240  (1941)
!     ------------------------------------------------------------------
!
!
!
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
    complex(kind=8) :: c8val
    character(len=24) :: nomdia, ualf
!     ------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadia, ibloc, ide, iderbl, iequa, ilong
    integer(kind=8) :: isol, ixx, ldiag, lmat, nbbloc, nbsol
!
!-----------------------------------------------------------------------
    data ualf/'                   .UALF'/
    data nomdia/'                   .&VDI'/
!     ------------------------------------------------------------------
!
    call jemarq()
    ualf(1:19) = nommat
    nomdia(1:19) = nommat
!
!     --- CREATION/RAPPEL D'UN TABLEAU POUR STOCKER LA DIAGONALE -------
    call wkvect(nomdia, 'V V C', neq, ldiag)
!
!     ------------------------------------------------------------------
!     --- PREMIERE  PARTIE : RESOLUTION DESCENDANTE ---
!     --- ET REMPLI LE TABLEAU DIAGONAL POUR L'ETAPE SUIVANTE
!     --- DANS L'ETAPE DE RESOLUTION DESCENDANTE, SEULE EST
!     --- CONCERNEE LA PARTIE INFERIEURE DE LA MATRICE
!
    do ibloc = 1, nbbloc
        if (ablo(ibloc) .ge. neq) goto 101
!
!        -- IDERBL EST LE NUMERO DU DERNIER BLOC A TRAITER:
!
!-     ACCES EN LECTURE AU BLOC INF COURANT
!
        iderbl = ibloc+nbbloc
        call jeveuo(jexnum(ualf, ibloc+nbbloc), 'L', lmat)
        do iequa = ablo(ibloc)+1, ablo(ibloc+1)
            if (iequa .gt. neq) goto 111
            ilong = hcol(iequa)
!-  ADRESSE DANS LE BLOC INF DU TERME DIAGONAL CORRESPONDANT A
!-  L'EQUATION IEQUA
            iadia = lmat+adia(iequa)-1
!-  ADRESSE DANS LE BLOC INF DU PREMIER TERME NON NUL SUR LA
!-  LIGNE IEQUA
            ide = iadia-ilong+1
!-  INDICE DU PREMIER TERME NON NUL SUR LA LIGNE IEQUA
            ixx = iequa-ilong+1
            zc(ldiag+iequa-1) = zc(iadia)
            do isol = 1, nbsol
                c8val = (0.d0, 0.d0)
                do i = 0, ilong-2
                    c8val = c8val+xsol(ixx+i, isol)*zc(ide+i)
                end do
                xsol(iequa, isol) = xsol(iequa, isol)-c8val
                xsol(iequa, isol) = xsol(iequa, isol)/zc(ldiag+iequa-1)
            end do
        end do
111     continue
        call jelibe(jexnum(ualf, ibloc+nbbloc))
    end do
101 continue
!
!
!     --- DEUXIEME  PARTIE : RESOLUTION REMONTANTE ---
!     --- UTILISATION DES BLOCS SUP DE LA MATRICE   ---
    iderbl = nbbloc+nbbloc
    do ibloc = iderbl, nbbloc+1, -1
        call jeveuo(jexnum(ualf, ibloc-nbbloc), 'L', lmat)
        do iequa = ablo(ibloc-nbbloc+1), ablo(ibloc-nbbloc)+1, -1
            if (iequa .gt. neq) goto 310
            ilong = hcol(iequa)
!-  ADRESSE DANS LE BLOC SUP DU TERME DIAGONAL CORRESPONDANT A
!-  L'EQUATION IEQUA
            iadia = lmat+adia(iequa)-1
!-  ADRESSE DANS LE BLOC SUP DU PREMIER TERME NON NUL SUR LA
!-  COLONNE IEQUA
            ide = iadia-ilong+1
!-  INDICE DU PREMIER TERME NON NUL SUR LA COLONNE IEQUA
            ixx = iequa-ilong+1
            do isol = 1, nbsol
                c8val = -xsol(iequa, isol)
                if (c8val .ne. 0) then
                    do i = 0, ilong-2
                        xsol(ixx+i, isol) = xsol(ixx+i, isol)+c8val*zc( &
                                            ide+i)
!                    XSOL(IXX+I,ISOL)=XSOL(IXX+I,ISOL)+C8VAL*ZC(IDE+I)/
!     1              ZC(LDIAG+IXX+I-1)
                    end do
                end if
            end do
310         continue
        end do
        call jelibe(jexnum(ualf, ibloc-nbbloc))
    end do
!
    call jedetr(nomdia)
!
    call jedema()
end subroutine
