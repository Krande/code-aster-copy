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

subroutine bsthco(nomte, bsigth, indith)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btsig.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrc2.h"
#include "asterfort/moytem.h"
#include "asterfort/promat.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "asterfort/verifm.h"
#include "asterfort/vexpan.h"
    real(kind=8) :: bsigth(51)
    aster_logical :: indith
    character(len=16) :: nomte
!
!      CALCUL DU BSIGMA POUR LES CONTRAINTES THERMIQUES
!      (I.E. BT*D*ALPHA(T-TREF)) POUR LES ELEMENTS DE COQUE (COQUE_3D)
!     ------------------------------------------------------------------
!     IN  NOMTE  : NOM DU TYPE D'ELEMENT
!     OUT BSIGTH : BT*SIGMA POUR LES CONTRAINTES THERMIQUES
!     OUT INDITH : LOGICAL = .TRUE.  YA DES DEFORMATIONS THERMIQUES
!                          = .FALSE. SINON
!     ------------------------------------------------------------------
    integer(kind=8) :: icara, icompo, icou, imate, inte, intsn, intsr, jgeom, lzi
    integer(kind=8) :: lzr, nb1, nb2, nbcou, npge, npgsn, npgsr, kwgt, itab(8), iret
    parameter(npge=3)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: epsth(5), sigmth(5), bsigt1(42)
    real(kind=8) :: ksi3s2, kappa, matc(5, 5)
    real(kind=8) :: coef, deux, epais, eptot, quatre, trois, un, valpar, zero
    real(kind=8) :: zic, zmin, epsthe
!     ------------------------------------------------------------------
!
! --- INITIALISATIONS :
!     ---------------
    zero = 0.0d0
    un = 1.0d0
    deux = 2.0d0
    trois = 3.0d0
    quatre = 4.0d0
!
    indith = .false.
!
    epsth(1) = zero
    epsth(2) = zero
    epsth(3) = zero
    epsth(4) = zero
    epsth(5) = zero
!
    sigmth(1) = zero
    sigmth(2) = zero
    sigmth(3) = zero
    sigmth(4) = zero
    sigmth(5) = zero
!
    bsigt1(1:42) = zero
    bsigth(1:51) = zero
!
! --- RECUPERATION DE L'OBJET .DESI :
!     -----------------------------
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
!
! --- NOMBRE DE NOEUDS (NB1 : SERENDIP, NB2 : LAGRANGE) :
!     -------------------------------------------------
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
!
! --- NOMBRE DE POINTS D'INTEGRATION DANS LE PLAN MOYEN
! --- (INTEGRATION REDUITE) :
!     ---------------------
    npgsr = zi(lzi-1+3)
!
! --- NOMBRE DE POINTS D'INTEGRATION DANS LE PLAN MOYEN
! --- (INTEGRATION NORMALE) :
!     ---------------------
    npgsn = zi(lzi-1+4)
!
! --- RECUPERATION DE L'OBJET .DESR :
!     -----------------------------
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
! --- RECUPERATION DE LA CARTE DE COMPORTEMENT :
!     ----------------------------------------
    call tecach('NNO', 'PCOMPOR', 'L', iret, nval=8, &
                itab=itab)
    icompo = itab(1)
    if (icompo .eq. 0) then
        nbcou = 1
    else
!
! ------ NOMBRE DE COUCHES :
!        -----------------
        read (zk16(icompo+6-1), '(I3)') nbcou
!
        if (nbcou .le. 0) then
            call utmess('F', 'ELEMENTS_12')
        end if
!
    end if
!
!
! --- CARACTERISTIQUES DE COQUES :
!     --------------------------
    call jevech('PCACOQU', 'L', icara)
! ---   EPAISSEUR TOTALE :
    eptot = zr(icara)
! ---   COORDONNEE MINIMALE SUIVANT L'EPAISSEUR
    zmin = -eptot/deux
! ---   COEFFICIENT DE CORRECTION DU CISAILLEMENT
    kappa = zr(icara+3)
! ---   EPAISSEUR D'UNE COUCHE
    epais = eptot/nbcou
!
! --- RECUPERATION DES COORDONNEES DES NOEUDS DANS LA GEOMETRIE
! --- INITIALE :
!     --------
    call jevech('PGEOMER', 'L', jgeom)
!
! --- RECUPERATION DU MATERIAU :
!     ------------------------
    call jevech('PMATERC', 'L', imate)
!______________________________________________________________________
!
!---- RECUPERATION DE LA TEMPERATURE
!
!
    indith = .true.
!______________________________________________________________________
!
! --- DETERMINATION DES REPERES LOCAUX AUX NOEUDS DANS LA
! --- CONFIGURATION INITIALE
! --- VECTA DESIGNE LES VECTEURS COVARIANTS DANS LE PLAN MOYEN A
! ---       CHAQUE NOEUD
! --- VECTN DESIGNE LES VECTEURS NORMAUX AU PLAN MOYEN
! --- VECTPT DESIGNE LES REPERES LOCAUX ORTHORNORMES EN CHAQUE
! --- NOEUD DANS LA CONFIGURATION INITIALE :
!     ------------------------------------
    call vectan(nb1, nb2, zr(jgeom), zr(lzr), vecta, &
                vectn, vectpt)
!
! --- COMPTEUR SERVANT A L'INTEGRATION :
!     --------------------------------
    kwgt = 0
!
! --- BOUCLE SUR LES COUCHES :
!     ----------------------
    do icou = 1, nbcou
!
! ---   BOUCLE SUR LES POINTS D'INTEGRATION DANS L'EPAISSEUR :
!       ----------------------------------------------------
        do inte = 1, npge
!
! ---      POSITION DANS L'EPAISSEUR :
            if (inte .eq. 1) then
                zic = zmin+(icou-1)*epais
                coef = un/trois
            else if (inte .eq. 2) then
                zic = zmin+epais/deux+(icou-1)*epais
                coef = quatre/trois
            else if (inte .eq. 3) then
                zic = zmin+epais+(icou-1)*epais
                coef = un/trois
            end if
! ---      COORDONNEE ISOPARAMETRIQUE DANS L'EPAISSEUR DIVISEE PAR 2
            ksi3s2 = zic/epais
!
! ---      CALCUL POUR L'INTEGRATION REDUITE DES PARTIES MEMBRANE
! ---      BTDM ET CISAILLEMENT BTDS DE LA MATRICE B :
!          =========================================
!
! ---      BOUCLE SUR LES POINTS DE L'INTEGRATION REDUITE :
!          ----------------------------------------------
            do intsr = 1, npgsr
!
! ---       .D'UNE PART :
! ---        DETERMINATION DES REPERES LOCAUX AUX POINTS D'INTEGRATION
! ---        DANS LA CONFIGURATION INITIALE
! ---        VECTG DESIGNE LES VECTEURS COVARIANTS DANS LE PLAN MOYEN
! ---              EN CHAQUE POINT D'INTEGRATION
! ---        VECTT DESIGNE LES REPERES LOCAUX ORTHORNORMES EN CHAQUE
! ---        POINT D'INTEGRATION DANS LA CONFIGURATION INITIALE
! ---       .D'AUTRE-PART :
! ---        SOIT H LA MATRICE DE PASSAGE DU TENSEUR DE GREEN-LAGRANGE
! ---        DU REPERE LOCAL AU REPERE GLOBAL
! ---        SOIT S LA MATRICE CONSTANTE TELLE QUE [S]*(DU/DX)
! ---        REPRESENTE LA PARTIE LINEAIRE DU TENSEUR DE GREEN-LAGRANGE
! ---        ON CALCULE LES PRODUITS [HSFM] = [H]*[S] POUR LA PARTIE
! ---                                MEMBRANE-FLEXION
! ---                                [HSS]  = [H] * [S] POUR LA PARTIE
! ---                                CISAILLEMENT :
!           -------------------------------------
                call mahsms(0, nb1, zr(jgeom), ksi3s2, intsr, &
                            zr(lzr), epais, vectn, vectg, vectt, &
                            hsfm, hss)
!
! ---       MULTIPLICATION DES MATRICES [HSFM] ET [HSS] PAR L'INVERSE
! ---       DE LA MATRICE JACOBIENNE [JM1]:
! ---       [HSJ1M] = [HSFM]*[JM1] , [HSJ1S] = [HSS]*[JM1] :
!           ----------------------------------------------
                call hsj1ms(epais, vectg, vectt, hsfm, hss, &
                            hsj1m, hsj1s)
!
! ---       CALCUL POUR L'INTEGRATION REDUITE DES PARTIES MEMBRANE
! ---       BTDM ET CISAILLEMENT BTDS DE LA MATRICE B :
!           -----------------------------------------
                call btdmsr(nb1, nb2, ksi3s2, intsr, zr(lzr), &
                            epais, vectpt, hsj1m, hsj1s, btdm, &
                            btds)
            end do
!
! ---      CALCUL POUR L'INTEGRATION NORMALE DE LA PARTIE FLEXION
! ---      BTDFN DE LA MATRICE B :
!          =====================
!
! ---      BOUCLE SUR LES POINTS DE L'INTEGRATION NORMALE
! ---      DANS LE PLAN MOYEN :
!          ------------------
            do intsn = 1, npgsn
!
! ---       .D'UNE PART :
! ---        DETERMINATION DES REPERES LOCAUX AUX POINTS D'INTEGRATION
! ---        DANS LA CONFIGURATION INITIALE
! ---        VECTG DESIGNE LES VECTEURS COVARIANTS DANS LE PLAN MOYEN
! ---              EN CHAQUE POINT D'INTEGRATION
! ---        VECTT DESIGNE LES REPERES LOCAUX ORTHORNORMES EN CHAQUE
! ---        POINT D'INTEGRATION DANS LA CONFIGURATION INITIALE
! ---       .D'AUTRE-PART :
! ---        SOIT H LA MATRICE DE PASSAGE DU TENSEUR DE GREEN-LAGRANGE
! ---        DU REPERE LOCAL AU REPERE GLOBAL
! ---        SOIT S LA MATRICE CONSTANTE TELLE QUE [S]*(DU/DX)
! ---        REPRESENTE LA PARTIE LINEAIRE DU TENSEUR DE GREEN-LAGRANGE
! ---        ON CALCULE LE PRODUIT [HSF] = [H]*[S] POUR LA PARTIE
! ---                              FLEXION :
!           ------------------------------
                call mahsf(1, nb1, zr(jgeom), ksi3s2, intsn, &
                           zr(lzr), epais, vectn, vectg, vectt, &
                           hsf)
!
! ---       MULTIPLICATION DE LA MATRICE [HSF] PAR L'INVERSE
! ---       DE LA MATRICE JACOBIENNE [JM1]:
! ---       [HSJ1FX] = [HSF]*[JM1]  :
!           ----------------------
                call hsj1f(intsn, zr(lzr), epais, vectg, vectt, &
                           hsf, kwgt, hsj1fx, wgt)
!
! ---       PRODUIT DU POIDS DU POINT DE GAUSS DANS L'EPAISSEUR PAR WGT
! ---       QUI EST LE PRODUIT DU POIDS DU POINT DE GAUSS COURANT
! ---       DANS L'EPAISSEUR PAR LE JACOBIEN EN CE POINT :
!           --------------------------------------------
                wgt = coef*wgt
!
! ---       CALCUL POUR L'INTEGRATION NORMALE DE LA PARTIE FLEXION
! ---       BTDF DE LA MATRICE B :
!           --------------------
                call btdfn(1, nb1, nb2, ksi3s2, intsn, &
                           zr(lzr), epais, vectpt, hsj1fx, btdf)
!
! ---       CALCUL DE LA MATRICE B [BTILD] PAR INTEGRATION SELECTIVE
! ---       ET INSERTION DES PARTIES [BTDM] ET [BDTS] ET INSERTION
! ---       DE LA PARTIE [BTDF]  :
!           -------------------
                call btdmsn(1, nb1, intsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)
!
! ---       EVALUATION DES DEFORMATIONS THERMIQUES :
!           ======================================
                call verifm('RIGI', inte, 3, '+', zi(imate), &
                            epsthe)
                call moytem('RIGI', inte, 3, '+', valpar, &
                            iret)
!
                epsth(1) = epsthe
                epsth(2) = epsthe
!
! ---       CALCUL DE LA MATRICE DE COMPORTEMENT  MATC(5,5) :
!           ----------------------------------------------
                call matrc2(1, 'TEMP    ', [valpar], kappa, matc, &
                            vectt)
!
! ---       CALCUL DES CONTRAINTES THERMIQUES SIGMTH(5) :
!           -------------------------------------------
                call promat(matc, 5, 5, 5, epsth, &
                            5, 5, 1, sigmth)
!
! ---       CALCUL DES FORCES INTERNES DUES AUX CONTRAINTES THERMIQUES :
!           ----------------------------------------------------------
                call btsig(5*nb1+2, 5, wgt, btild, sigmth, &
                           bsigt1)
!
            end do
        end do
    end do
!
! --- EXPANSION DE BSIGT1 DANS BSIGTH :
!     -------------------------------
    call vexpan(nb1, bsigt1, bsigth)
    bsigth(6*nb1+1) = bsigt1(5*nb1+1)
    bsigth(6*nb1+2) = bsigt1(5*nb1+2)
!
end subroutine
