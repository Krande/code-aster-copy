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
subroutine tldur8(nommat, hcol, adia, ablo, npivot, &
                  neq, nbbloc, ildeb, ilfin, eps)
    implicit none
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
! BUT : DECOMPOSTION D'UNE MATRICE NON_SYMETRIQUE A COEFFICIENTS REELS
!       SOUS LA FORME DE CROUT LDU
!  ------------------------------------------------------------------
!
!  IN  NOMMAT  :    : NOM UTILISATEUR DE LA MATRICE A FACTORISER
!  IN  HCOL    : IS : HCOL DE LA MATRICE
!  HCOL(I) RENVOIE LA HAUTEUR DE LA I-EME COLONNE
!  IN  ADIA    : IS : ADRESSE DU TERME DIAGONALE DANS SON BLOC
!  ADIA(I) RENVOIE L'ADRESSE DE LA I-EME LIGNE DANS SON BLOC
!  IN  ABLO    :  :   POINTEUR DE BLOC
!  ABLO(I+1) RENVOIE LE NO DE LA DERNIERE LIGNE DU I-EME BLOC
!
!  LES TABLEAUX RELATIFS AUX BLOCS ONT ETE CONSTRUITS A PARTIR
!  DES NBBLOCS PREMIERS BLOCS (I.E LES BLOCS SUP)
!  MAIS ILS RESTENT VALABLES POUR LES NBBLOCS BLOCS SUIVANTS
!  (I.E LES BLOCS INF)
!  PUISQUE LES POFILS SUP ET INF SONT IDENTIQUES
!
!  VAR PIVOT   : IS :
!  : EN SORTIE : NPIVOT  = 0 ==> R.A.S.
!  :             NPIVOT  > 0 ==> MATRICE SINGULIERE
!                                POUR L'EQUATION DE NUMERO NPIVOT
!  :             NPIVOT  < 0 ==> -NPIVOT TERMES DIAGONAUX < 0
!
!  IN  NEQ     : IS : NOMBRE TOTAL D'EQUATION
!  IN  NBBLOC  : IS : NOMBRE DE BLOC DE LA MATRICE
!  IN  ILDEB   : IS : NUMERO DE LA LIGNE DE DEPART DE LA FACTORISATION
!  IN  ILFIN   : IS : NUMERO DE LA LIGNE DE FIN DE FACTORISITION
!  ------------------------------------------------------------------
!
!  CREATION DE DEUX OBJETS DE TRAVAIL (SUR LA VOLATILE)
!  1)  UN TABLEAU POUR LA DIAGONALE
!  2)  UN TABLEAU POUR LA COLONNE COURANTE
!
!
!  --- RAPPEL SOMMAIRE DE L'ALGORITHME ------------------------------
!
!  POUR S = 2,3, ... ,N
!  !  POUR I = 1,2, ... ,S-1
!  !  !  POUR M = 1,2, ... ,I-1
!  !  !  !  K(S,I) = K(S,I) - K(S,M)*K(M,I) % MODIFIE   LA LIGNE
!  !  !  !  K(I,S) = K(I,S) - K(I,M)*K(M,S) % MODIFIE   LA COLONNE
!  !  !  FIN_POUR
!  !  !  K(I,S) = K(I,S)/K(I,I)           % NORMALISATION DE LA COLONNE
!  !  FIN_POUR
!  !  POUR M = 1,2, ... ,S-1
!  !  !  K(S,S) = K(S,S) - K(S,M) * K(M,S)  % MODIFICATION DU PIVOT
!  !  FIN_POUR
!  FIN_POUR
!  ------------------------------------------------------------------
!  REFERENCE (HISTORIQUE) :
!  (1) P.D. CROUT,
!      A SHORT METHOD FOR EVALUATING DETERMINANTS AND SOLVING SYSTEMS
!      OF LINEAR EQUATIONS WITH REAL OR COMPLEX COEFFICIENTS.
!      AIEE TRANSACTION VOL 60, PP 1235-1240  (1941)
!  ------------------------------------------------------------------
!
!
!  ------------------------------------------------------------------
    integer(kind=8) :: hcol(*), adia(*), ablo(*)
    integer(kind=8) :: npivot, neq, nbbloc, ildeb, ilfin
    integer(kind=8) :: ldiag, ibloc, il1, il2, iaa, il
    integer(kind=8) :: iaas, iaai, kl1, imini, iequa, i, jblmin, jbloc, jl1, jl2
    integer(kind=8) :: iabs, iabi, ilong, iadiai, idei, iadias, ides, idl, jnmini
    integer(kind=8) :: jequa, jlong, jadias, jdes, jadiai, jdei, jdl, ibcl1
    integer(kind=8) :: lm, icai, icas, icbi, icbs, icd
!
    real(kind=8) :: eps
    real(kind=8) :: r8vali, r8vals
    character(len=24) :: nomdia, ualf
    character(len=19) :: noma19
    character(len=24), pointer :: refa(:) => null()
    real(kind=8), pointer :: digs(:) => null()
!  ------------------------------------------------------------------
    data ualf/'                   .UALF'/
    data nomdia/'                   .&VDI'/
!  ------------------------------------------------------------------
!
    call jemarq()
!
    noma19 = nommat
    ualf(1:19) = nommat
    nomdia(1:19) = nommat
!
!  --- CREATION D'UN TABLEAU POUR STOCKER LA DIAGONALE
    call wkvect(nomdia, 'V V R', neq, ldiag)
    call jeveuo(noma19//'.DIGS', 'E', vr=digs)
!
!  --- INITIALISATIONS ET ALLOCATION ---
!
    npivot = 0
!
!  ------- BOUCLE SUR LES BLOCS A TRANSFORMER -----------------------
    do ibloc = 1, nbbloc
!
        il1 = ablo(ibloc)+1
        il2 = ablo(ibloc+1)
!
!
        if (il2 .lt. ildeb) then
!        --- C'EST TROP TOT : MAIS ON REMPLIT LA DIAGONALE ---
            call jeveuo(jexnum(ualf, ibloc), 'L', iaa)
            do il = il1, il2
                zr(ldiag+il-1) = zr(iaa+adia(il)-1)
            end do
            call jelibe(jexnum(ualf, ibloc))
            goto 150
        else if (il1 .gt. ilfin) then
!        --- C'EST FINI ---
            goto 160
        else
!
!- RECUPERATION DES BLOCS SUP ET INF COURANTS
!
            call jeveuo(jexnum(ualf, ibloc), 'E', iaas)
            call jeveuo(jexnum(ualf, ibloc+nbbloc), 'E', iaai)
            if (il1 .lt. ildeb) then
                kl1 = ildeb
                do il = il1, kl1-1
                    zr(ldiag+il-1) = zr(iaai+adia(il)-1)
                end do
            else
                kl1 = il1
            end if
            if (il2 .gt. ilfin) il2 = ilfin
        end if
!
!
!     --- RECHERCHE DE LA PLUS PETITE EQUATION EN RELATION AVEC UNE
!     --- EQUATION DES LIGNES EFFECTIVES DU BLOC COURANT
        imini = il2-1
        do iequa = kl1, il2
            imini = min(iequa-hcol(iequa), imini)
        end do
        imini = imini+1
!
!     --- RECHERCHE DU BLOC D'APPARTENANCE DE L'EQUATION IMINI ---
        do i = 1, ibloc
            jblmin = i
            if (ablo(1+i) .ge. imini) goto 50
        end do
50      continue
!
!     --- BOUCLE  SUR  LES  BLOCS  DEJA  TRANSFORMES ---
        do jbloc = jblmin, ibloc-1
!
            jl1 = max(imini, ablo(jbloc)+1)
            jl2 = ablo(jbloc+1)
            call jeveuo(jexnum(ualf, jbloc), 'L', iabs)
            call jeveuo(jexnum(ualf, jbloc+nbbloc), 'L', iabi)
!
            do iequa = kl1, il2
!
!           --- RECUPERATION DE L'ADRESSE ET LA LONGUEUR DE LA LIGNE
                ilong = hcol(iequa)-1
!-   ADRESSE DU TERME DIAGONAL DE LA LIGNE IEQUA DANS LE BLOC INF
                iadiai = iaai+adia(iequa)-1
!-   ADRESSE DU DEBUT DE LA LIGNE IEQUA DANS LE BLOC INF
                idei = iadiai-ilong
!-   ADRESSE DU TERME DIAGONAL DE LA COLONNE IEQUA DANS LE BLOC SUP
                iadias = iaas+adia(iequa)-1
!-   ADRESSE DU DEBUT DE LA COLONNE IEQUA DANS LE BLOC SUP
                ides = iadias-ilong
!-   INDICE DU DEBUT DE LA LIGNE (COLONNE) IEQUA
                idl = iequa-ilong
!
!           --- UTILISATION DES LIGNES (IDL+1) A (JL2) ---
!
!- MODIFICATION DES LIGNES   IDLI+1 A JL2
!-          ET  DES COLONNES IDLS+1 A JL2
!- POUR LES BLOCS SUP ET INF COURANTS AUXQUELS ON A ACCEDE EN ECRITURE
!
                jnmini = max(idl, jl1)
                do jequa = jnmini, jl2
                    jlong = hcol(jequa)-1
!-    ADRESSE DU TERME DIAGONAL KII DANS LE BLOC SUP
                    jadias = iabs+adia(jequa)-1
!-    ADRESSE DANS LE BLOC SUP DU PREMIER TERME NON NUL DE LA COLONNE I
                    jdes = jadias-jlong
!-    ADRESSE DU TERME DIAGONAL KII DANS LE BLOC INF
                    jadiai = iabi+adia(jequa)-1
!-    ADRESSE DANS LE BLOC INF DU PREMIER TERME NON NUL DE LA LIGNE I
                    jdei = jadiai-jlong
!-    INDICE DU PREMIER TERME NON NUL DE LA COLONNE (LIGNE) I
                    jdl = jequa-jlong
!-    INDICE DU PREMIER TERME A PARTIR DUQUEL ON VA FAIRE LE PRODUIT
!-    SCALAIRE K(S,M)*K(M,I) DES TERMES LIGNE-COLONNE
                    ibcl1 = max(idl, jdl)
!-    LONGUEUR SUR LAQUELLE ON FAIT LE PRODUIT SCALAIRE LIGNE-COLONNE
                    lm = jequa-ibcl1
!-    ADRESSE DANS LE BLOC INF DU PREMIER TERME A MODIFIER SUR LA LIGNE
!-    IEQUA
                    icai = idei+ibcl1-idl
!-    ADRESSE DANS LE BLOC SUP DU PREMIER TERME A MODIFIER SUR LA
!-    COLONNE IEQUA
                    icas = ides+ibcl1-idl
!-    ADRESSE DANS LE BLOC INF DU PREMIER TERME DE LA LIGNE
!-    A PARTIR DUQUEL ON FAIT LE PRODUIT SCALAIRE POUR MODIFIER LA
!-    COLONNE = K(I,M)
                    icbi = jdei+ibcl1-jdl
!-    ADRESSE DANS LE BLOC SUP DU PREMIER TERME DE LA COLONNE
!-    A PARTIR DUQUEL ON FAIT LE PRODUIT SCALAIRE POUR MODIFIER LA
!-    LIGNE = K(M,I)
                    icbs = jdes+ibcl1-jdl
!-    TERME COURANT DE LA LIGNE S=IEQUA A MODIFIER ( = K(S,I))
                    r8vali = zr(icai+lm)
!-    TERME COURANT DE LA COLONNE S=IEQUA A MODIFIER ( = K(I,S))
                    r8vals = zr(icas+lm)
                    do i = 0, lm-1
!-                   K(S,I) = K(S,I) -  K(S,M)   * K(M,I)
                        r8vali = r8vali-zr(icai+i)*zr(icbs+i)
!-                   K(I,S) = K(I,S) -  K(M,S)   * K(I,M)
                        r8vals = r8vals-zr(icas+i)*zr(icbi+i)
                    end do
                    zr(icai+lm) = r8vali
                    zr(icas+lm) = r8vals
                    icd = ldiag+jequa-1
!                    ZR(ICAI+LM) = ZR(ICAI+LM) / ZR(ICD)
                    zr(icas+lm) = zr(icas+lm)/zr(icd)
!
                end do
            end do
!
!-  DEVEROUILLAGE DES BLOCS SUP ET INF DEJA TRANSFORMES
!
            call jelibe(jexnum(ualf, jbloc))
            call jelibe(jexnum(ualf, jbloc+nbbloc))
        end do
!
!     --- UTILISATION DU BLOC EN COURS DE TRANSFORMATION ---
        jl1 = max(imini, il1)
        do iequa = kl1, il2
!
!        --- RECUPERATION DE L ADRESSE ET LA LONGUEUR DE LA LIGNE
!            (COLONNE) ---
!           IADIAI : ADRESSE DU TERME DIAGONAL COURANT DS LE BLOC INF
!           IADIAS : ADRESSE DU TERME DIAGONAL COURANT DS LE BLOC SUP
!           IDEI   : ADRESSE DU DEBUT DE LA LIGNE COURANTE
!           IDES   : ADRESSE DU DEBUT DE LA COLONNE COURANTE
!           IDL   : 1-ER DDL A VALEUR NON NULLE DANS LA LIGNE (COLONNE)
            ilong = hcol(iequa)-1
!-  ADRESSE DE K(S,S) DANS LE BLOC INF
            iadiai = iaai+adia(iequa)-1
!-  ADRESSE DU PREMIER TERME NON NUL SUR LA LIGNE S (IEQUA) DANS
!-  LE BLOC INF
            idei = iadiai-ilong
!-  ADRESSE DE K(S,S) DANS LE BLOC SUP
            iadias = iaas+adia(iequa)-1
!-  ADRESSE DU PREMIER TERME NON NUL SUR LA COLONNE S (IEQUA)
!-  DANS LE BLOC SUP
            ides = iadias-ilong
!-  INDICE DU PREMIER TERME NON NUL SUR LA LIGNE (COLONNE) S (IEQUA)
            idl = iequa-ilong
!
!        --- UTILISATION DES LIGNES (IDL+1) A (IEQUA-1) ---
            jnmini = max(iequa-ilong, jl1)
            do jequa = jnmini, iequa-1
                jlong = hcol(jequa)-1
!-      ADRESSE DE K(I,I) DANS LE BLOC INF
                jadiai = iaai+adia(jequa)-1
!-      ADRESSE DU PREMIER TERME NON NUL SUR LA LIGNE I (JEQUA)
!-      DANS LE BLOC INF
                jdei = jadiai-jlong
!-      ADRESSE DE K(I,I) DANS LE BLOC SUP
                jadias = iaas+adia(jequa)-1
!-      ADRESSE DU PREMIER TERME NON NUL SUR LA COLONNE I (JEQUA)
!-      DANS LE BLOC SUP
                jdes = jadias-jlong
!-      INDICE DU PREMIER TERME NON NUL SUR LA LIGNE (COLONNE) I
!-      (JEQUA)
                jdl = jequa-jlong
!-      INDICE DU PREMIER TERME A PARTIR DUQUEL ON VA FAIRE LES
!-      PRODUITS SCALAIRES (LIGNE I X COLONNE S) ET
!-      (LIGNE S X COLONNE I)
                ibcl1 = max(idl, jdl)
!-      LONGUEUR SUR LAQUELLE ON FAIT LES PRODUITS SCALAIRES
!-      VECTEUR LIGNE X VECTEUR COLONNE
                lm = jequa-ibcl1
!-      ADRESSE DANS LE BLOC INF DU PREMIER TERME A PARTIR DUQUEL ON
!-      FAIT LE PRODUIT SCALAIRE POUR CALCULER LE TERME LIGNE K(S,I)
                icai = idei+ibcl1-idl
!-      ADRESSE DANS LE BLOC SUP DU PREMIER TERME A PARTIR DUQUEL ON
!-      FAIT LE PRODUIT SCALAIRE POUR CALCULER LE TERME LIGNE K(S,I)
                icbs = jdes+ibcl1-jdl
!-      ADRESSE DANS LE BLOC SUP DU PREMIER TERME A PARTIR DUQUEL ON
!-      FAIT LE PRODUIT SCALAIRE POUR CALCULER LE TERME COLONNE
                icas = ides+ibcl1-idl
!-      ADRESSE DANS LE BLOC INF DU PREMIER TERME A PARTIR DUQUEL ON
!-      FAIT LE PRODUIT SCALAIRE POUR CALCULER LE TERME COLONNE
                icbi = jdei+ibcl1-jdl
!-      TERME COURANT DE LA LIGNE S (IEQUA) = K(S,I)
                r8vali = zr(icai+lm)
!-      TERME COURANT DE LA COLONNE S (IEQUA) = K(I,S)
                r8vals = zr(icas+lm)
                do i = 0, lm-1
!-                K(S,I) = K(S,I) -  K(S,M)   * K(M,I)
                    r8vali = r8vali-zr(icai+i)*zr(icbs+i)
!-                K(I,S) = K(I,S) -  K(M,S)   * K(I,M)
                    r8vals = r8vals-zr(icas+i)*zr(icbi+i)
                end do
                zr(icai+lm) = r8vali
                zr(icas+lm) = r8vals
!
!        --- UTILISATION DE LA LIGNE IEQUA (CALCUL DU PIVOT) ---
                icd = ldiag+jequa-1
                zr(icas+lm) = zr(icas+lm)/zr(icd)
            end do
!
!
!        --- CALCUL DU TERME DIAGONAL ---
            lm = ilong-1
!-      INDICE DU PREMIER TERME NON NUL DE LA LIGNE COURANTE IEQUA
!-      DS LE BLOC
            icai = iadiai-ilong
            icas = iadias-ilong
!
!        --- SAUVEGARDE DE LA COLONNE ---
!        --- NORMALISATION DE LA LIGNE ---
!-      INDICE DU PREMIER TERME DIAGONAL DANS LE VECTEUR DES TERMES
!-      DIAGONAUX NON ENCORE DECOMPOSES CORRESPONDANT AU PREMIER
!-      TERME DE LA LIGNE COURANTE A NORMALISER
            icd = ldiag+iequa-ilong-1
            r8vali = zr(iadiai)
            do i = 0, lm
                r8vali = r8vali-zr(icai+i)*zr(icas+i)
            end do
            zr(iadiai) = r8vali
            digs(neq+iequa) = zr(iadiai)
            zr(ldiag+iequa-1) = r8vali
!
!
!           --- LE PIVOT EST-IL NUL ? ----------------------------------
            if (abs(r8vali) .le. eps) then
                npivot = iequa
                goto 999
            end if
!
!
!           --- ON COMPTE LES PIVOTS NEGATIFS --------------------------
            if (r8vali .lt. 0.d0) then
                npivot = npivot-1
            end if
        end do
        call jelibe(jexnum(ualf, ibloc))
        call jelibe(jexnum(ualf, ibloc+nbbloc))
150     continue
    end do
!
160 continue
999 continue
!
!
    call jeveuo(noma19//'.REFA', 'E', vk24=refa)
    if (ilfin .eq. neq) then
        refa(8) = 'DECT'
    else
        refa(8) = 'DECP'
    end if
!
    call jedetr(nomdia)
!
    call jedema()
end subroutine
