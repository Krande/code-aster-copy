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
subroutine dstci2(dci, carat3, hft2, dfc, dmc, &
                  bca, an, am)
    implicit none
#include "asterfort/dstbfa.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dxtbm.h"
#include "asterfort/mgauss.h"
    real(kind=8) :: dci(2, 2), carat3(*), hft2(2, 6), dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bca(2, 3), an(3, 9), am(3, 6)
!.======================================================================
!
!  DSTCI2 -- DETERMINATION DES MATRICES AN ET AM QUI SONT TELLES QUE
!            ALPHA = AN*UN + AM*UM
!            POUR OBTENIR CETTE EXPRESSION ON PASSE PAR LA RELATION :
!            AA*ALPHA = (AW + AB)*UN + AL*UM
!            UN DESIGNE LES DEPLACEMENTS DE FLEXION (W,BETAX,BETAY)
!            UM DESIGNE LES DEPLACEMENTS DE MEMBRANE (UX,UY)
!            FORMELLEMENT :
!
!            |L4 0  0|   |L4C4 L4S4|
!   AA = 2/3*|0 L5  0| - |L5C5 L5S5|*DCI*(BCA - DFC_T*BFA)
!            |0  0 L6|   |L6C6 L6S6|
!
!           |L4C4 L4S4|
!    AB = - |L5C5 L5S5|*DCI*DFC_T*BFB
!           |L6C6 L6S6|
!
!           |L4C4 L4S4|
!    AL = - |L5C5 L5S5|*DCI*DMC_T*BM
!           |L6C6 L6S6|
!
!              |-2  L4C4 L4S4   2  L4C4 L4S4   0     0      0|
!    AW = -1/2*| 0     0    0  -2  L5C5 L5S5   2  L5C5   L5S5|
!              | 2  L6C6 L6S6   0     0    0  -2  L6C6   L6S6|
!
!
!   ARGUMENT        E/S  TYPE         ROLE
!    DCI(2,2)       IN    R       INVERSE DE LA MATRICE DE CISAILLEMENT
!                                 DE HOOKE
!    R(*)           IN    R       TABLEAU DE CARACTERISTIQUES
!                                 GEOMETRIQUES DE L'ELEMENT :
!                                 COS ET SIN DES ANGLES, LONGUEUR
!                                 DES COTES ,...
!    HFT2(2,6)      IN    R       MATRICE DEFINIE PAR T = [HF.T2]
!                                 T = [HF.T2] BETA,QSI2
!                                 OU T EST LE VECTEUR DES EFFORTS
!                                 TRANCHANTS ET BETA,QSI2 DESIGNE LE
!                                 VECTEUR DES DERIVEES SECONDES DES
!                                 ROTATIONS BETA PAR RAPPORT AUX
!                                 COORDONNEES DE L'ELEMENT DE REFERENCE
!    DF(3,3)        IN    R       MATRICE DE FLEXION DE HOOKE
!    DFC(3,2)       IN    R       MATRICE DE COUPLAGE
!                                 FLEXION-CISAILLEMENT DE HOOKE
!    DMC(3,2)       IN    R       MATRICE DE COUPLAGE
!                                 MEMBRANE-CISAILLEMENT DE HOOKE
!    BCA(2,3)       OUT   R       MATRICE RELIANT LES EFFORTS TRANCHANTS
!                                 AUX INCONNUES ALPHA : T = [BCA]*ALPHA
!                                 PLUS PRECISEMMENT  :
!                                 T = [HFT.T2]*[TA]*ALPHA
!    AN(3,9)        OUT   R       MATRICE RELIANT LES ROTATIONS ALPHA
!                                 AUX INCONNUES DE FLEXION UN
!    AM(3,6)        OUT   R       MATRICE RELIANT LES ROTATIONS ALPHA
!                                 AUX INCONNUES DE MEMBRANE UM
!     ------------------------------------------------------------------
    integer(kind=8) :: i, j, k, ic, int, iret
    real(kind=8) :: l(3), c(3), s(3), x(3), y(3), qsi(2), det
    real(kind=8) :: ta(6, 3), db(2, 3), aa(3, 3), aai(3, 3), aw(3, 9)
    real(kind=8) :: ab(3, 9), al(3, 6), dfcbfb(2, 9), dfcbfa(2, 3), bfa(3, 3)
    real(kind=8) :: bm(3, 6), bfb(3, 9), dmctbm(3, 6), dcidmc(2, 6)
    real(kind=8) :: dcidfb(2, 9)
    real(kind=8) :: zero, undemi, un, deux, trois, quatre, huit, rac3, eta
!     ------------------------------------------------------------------
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     ---------------
    zero = 0.0d0
    undemi = 0.5d0
    un = 1.0d0
    deux = 2.0d0
    trois = 3.0d0
    quatre = 4.0d0
    huit = 8.0d0
    rac3 = sqrt(trois)
    qsi(1) = un/rac3
    qsi(2) = -qsi(1)
    eta = zero
!
    do i = 1, 3
        do j = 1, 6
            am(i, j) = zero
            dmctbm(i, j) = zero
        end do
    end do
!
    do i = 1, 3
        do j = 1, 9
            an(i, j) = zero
            aw(i, j) = zero
        end do
    end do
!
    do i = 1, 2
        do j = 1, 9
            dfcbfb(i, j) = zero
        end do
    end do
!
    do i = 1, 2
        do j = 1, 3
            bca(i, j) = zero
        end do
    end do
!
    do i = 1, 6
        do j = 1, 3
            ta(i, j) = zero
        end do
    end do
!
    c(1) = carat3(16)
    c(2) = carat3(17)
    c(3) = carat3(18)
    s(1) = carat3(19)
    s(2) = carat3(20)
    s(3) = carat3(21)
    l(1) = carat3(13)
    l(2) = carat3(14)
    l(3) = carat3(15)
    x(1) = carat3(1)
    x(2) = carat3(2)
    x(3) = carat3(3)
    y(1) = carat3(4)
    y(2) = carat3(5)
    y(3) = carat3(6)
!
    ta(1, 1) = -huit*c(1)
    ta(2, 3) = -huit*c(3)
    ta(3, 1) = -quatre*c(1)
    ta(3, 2) = quatre*c(2)
    ta(3, 3) = -quatre*c(3)
    ta(4, 1) = -huit*s(1)
    ta(5, 3) = -huit*s(3)
    ta(6, 1) = -quatre*s(1)
    ta(6, 2) = quatre*s(2)
    ta(6, 3) = -quatre*s(3)
!
! --- CALCUL DU PRODUIT HFT2.TA :
!     -------------------------
    do j = 1, 3
        do k = 1, 6
            bca(1, j) = bca(1, j)+hft2(1, k)*ta(k, j)
            bca(2, j) = bca(2, j)+hft2(2, k)*ta(k, j)
        end do
    end do
!
!================================================================
! --- DETERMINATION DE LA MATRICE AA QUI EST TELLE QUE          =
! --- AA*ALPHA = (AW + AB)*UN + AL*UM                           =
! --- UN DESIGNE LES DEPLACEMENTS DE FLEXION (W,BETAX,BETAY)    =
! --- UM DESIGNE LES DEPLACEMENTS DE MEMBRANE (UX,UY)           =
! ---  FORMELLEMENT                                             =
! ---            |L4 0  0|   |L4C4 L4S4|                        =
! ---   AA = 2/3*|0 L5  0| - |L5C5 L5S5|*DCI*(BCA - DFC_T*BFA)  =
! ---            |0  0 L6|   |L6C6 L6S6|                        =
!================================================================
!
! --- BOUCLE SUR LES COTES DU TRIANGLE :
!      --------------------------------
    do ic = 1, 3
!
! ---     INITIALISATION DE DFCBFA :
!         ------------------------
        do i = 1, 2
            do j = 1, 3
                dfcbfa(i, j) = zero
            end do
        end do
!
! ---   BOUCLE SUR LES POINTS D'INTEGRATION DU COTE COURANT :
!       ---------------------------------------------------
        do int = 1, 2
!
! ---     CALCUL DE LA MATRICE BFA AU POINT D'INTEGRATION COURANT
! ---     RELIANT LES COURBURES AUX INCONNUES ALPHA
! ---     I.E. X = BFB*UN + BFA*ALPHA :
!         ---------------------------
            call dstbfa(qsi(int), eta, carat3, bfa)
!
! ---     CALCUL DU PRODUIT DFC_T*BFA :
!         ---------------------------
            do j = 1, 3
                do k = 1, 3
                    dfcbfa(1, j) = dfcbfa(1, j)+dfc(k, 1)*bfa(k, j)
                    dfcbfa(2, j) = dfcbfa(2, j)+dfc(k, 2)*bfa(k, j)
                end do
            end do
!
            do i = 1, 2
                do j = 1, 3
                    dfcbfa(i, j) = undemi*dfcbfa(i, j)
                end do
            end do
!
        end do
!     -------------------------------------------------------------
! --  FIN DE LA BOUCLE SUR LES POINTS D'INTEGRATION DU COTE COURANT
!     -------------------------------------------------------------
!
! ---   CALCUL DU PRODUIT DCI*(BCA - DFC_T*BFA) :
!       --------------------------------------
        do j = 1, 3
            db(1, j) = dci(1, 1)*(bca(1, j)-dfcbfa(1, j))+dci(1, 2)*(bca(2, j)-dfcbfa(2, j))
            db(2, j) = dci(2, 1)*(bca(1, j)-dfcbfa(1, j))+dci(2, 2)*(bca(2, j)-dfcbfa(2, j))
        end do
!
! ---               |L4 0  0|   |L4C4 L4S4|
! --- CALCUL DE 2/3*|0 L5  0| - |L5C5 L5S5|*DCI*(BCA - DFC_T*BFA) :
! ---               |0  0 L6|   |L6C6 L6S6|
! --- LES LKCK SONT DANS X , LES LKSK SONT DANS Y :
!     -------------------------------------------
        do j = 1, 3
            aa(ic, j) = -(x(ic)*db(1, j)+y(ic)*db(2, j))
        end do
        aa(ic, ic) = aa(ic, ic)+deux/trois*l(ic)
!
    end do
!     -------------------------------------------
! --  FIN DE LA BOUCLE SUR LES COTES DE L'ELEMENT
!     -------------------------------------------
!
!================================================================
! --- DETERMINATION DE LA MATRICE AB QUI EST TELLE QUE          =
! --- AA*ALPHA = (AW + AB)*UN + AL*UM                           =
! --- UN DESIGNE LES DEPLACEMENTS DE FLEXION (W,BETAX,BETAY)    =
! --- UM DESIGNE LES DEPLACEMENTS DE MEMBRANE (UX,UY)           =
! ---  FORMELLEMENT                                             =
! ---          |L4C4 L4S4|                                      =
! ---   AB = - |L5C5 L5S5|*DCI*DFC_T*BFB                        =
! ---          |L6C6 L6S6|                                      =
!================================================================
!
! --- CALCUL DE LA MATRICE BFB RELIANT LES COURBURES AUX INCONNUES
! --- DE FLEXION UN (X = BFB*UN+BFA*ALPHA) :
!     ------------------------------------
    call dstbfb(carat3(9), bfb)
!
! --- CALCUL DU PRODUIT DFC_T*BFB :
!     ---------------------------
    do j = 1, 9
        do k = 1, 3
            dfcbfb(1, j) = dfcbfb(1, j)+dfc(k, 1)*bfb(k, j)
            dfcbfb(2, j) = dfcbfb(2, j)+dfc(k, 2)*bfb(k, j)
        end do
    end do
!
! --- CALCUL DU PRODUIT DCI*DFC_T*BFB :
!     -------------------------------
    do j = 1, 9
        dcidfb(1, j) = dci(1, 1)*dfcbfb(1, j)+dci(1, 2)*dfcbfb(2, j)
        dcidfb(2, j) = dci(2, 1)*dfcbfb(1, j)+dci(2, 2)*dfcbfb(2, j)
    end do
!
! ---  CALCUL DE :
! ---          |L4C4 L4S4|
! ---   AB = - |L5C5 L5S5|*DCI*DFC_T*BFB
! ---          |L6C6 L6S6|
! --- LES LKCK SONT DANS X , LES LKSK SONT DANS Y :
!     -------------------------------------------
    do i = 1, 3
        do j = 1, 9
            ab(i, j) = -(x(i)*dcidfb(1, j)+y(i)*dcidfb(2, j))
        end do
    end do
!
!================================================================
! --- DETERMINATION DE LA MATRICE AL QUI EST TELLE QUE          =
! --- AA*ALPHA = (AW + AB)*UN + AL*UM                           =
! --- UN DESIGNE LES DEPLACEMENTS DE FLEXION (W,BETAX,BETAY)    =
! --- UM DESIGNE LES DEPLACEMENTS DE MEMBRANE (UX,UY)           =
! ---  FORMELLEMENT                                             =
! ---          |L4C4 L4S4|                                      =
! ---   AL = - |L5C5 L5S5|*DCI*DMC_T*BM                         =
! ---          |L6C6 L6S6|                                      =
!================================================================
!
! --- CALCUL DE LA MATRICE BM RELIANT LES DEFORMATIONS MEMBRANAIRES
! --- AUX DEPLACEMENTS DE MEMBRANE UM (I.E. (UX,UY) ) :
!     -----------------------------------------------
    call dxtbm(carat3(9), bm)
!
! --- CALCUL DU PRODUIT DMC_T*BM :
!     --------------------------
    do j = 1, 6
        do k = 1, 3
            dmctbm(1, j) = dmctbm(1, j)+dmc(k, 1)*bm(k, j)
            dmctbm(2, j) = dmctbm(2, j)+dmc(k, 2)*bm(k, j)
        end do
    end do
!
! --- CALCUL DU PRODUIT DCI*DMC_T*BM :
!     ------------------------------
    do j = 1, 6
        dcidmc(1, j) = dci(1, 1)*dmctbm(1, j)+dci(1, 2)*dmctbm(2, j)
        dcidmc(2, j) = dci(2, 1)*dmctbm(1, j)+dci(2, 2)*dmctbm(2, j)
    end do
!
! ---  CALCUL DE :
! ---          |L4C4 L4S4|
! ---   AL = - |L5C5 L5S5|*DCI*DMC_T*BM
! ---          |L6C6 L6S6|
! --- LES LKCK SONT DANS X , LES LKSK SONT DANS Y :
!     -------------------------------------------
    do i = 1, 3
        do j = 1, 6
            al(i, j) = -(x(i)*dcidmc(1, j)+y(i)*dcidmc(2, j))
        end do
    end do
!
!=================================
! --- INVERSION DE LA MATRICE AA =
!=================================
!
    do i = 1, 3
        do j = 1, 3
            aai(i, j) = zero
        end do
    end do
    do i = 1, 3
        aai(i, i) = un
    end do
    call mgauss('NFVP', aa, aai, 3, 3, &
                3, det, iret)
!
!===================================================================
! --- DETERMINATION DE LA MATRICE AW QUI EST TELLE QUE             =
! --- AA*ALPHA = (AW + AB)*UN + AL*UM                              =
! --- UN DESIGNE LES DEPLACEMENTS DE FLEXION (W,BETAX,BETAY)       =
! --- UM DESIGNE LES DEPLACEMENTS DE MEMBRANE (UX,UY)              =
! ---  FORMELLEMENT                                                =
! ---             |-2  L4C4 L4S4   2  L4C4 L4S4   0     0      0|  =
! ---   AW = -1/2*| 0     0    0  -2  L5C5 L5S5   2  L5C5   L5S5|  =
! ---             | 2  L6C6 L6S6   0     0    0  -2  L6C6   L6S6|  =
!===================================================================
!
    aw(1, 1) = un
    aw(1, 2) = -undemi*x(1)
    aw(1, 3) = -undemi*y(1)
    aw(1, 4) = -un
    aw(1, 5) = -undemi*x(1)
    aw(1, 6) = -undemi*y(1)
    aw(2, 4) = un
    aw(2, 5) = -undemi*x(2)
    aw(2, 6) = -undemi*y(2)
    aw(2, 7) = -un
    aw(2, 8) = -undemi*x(2)
    aw(2, 9) = -undemi*y(2)
    aw(3, 1) = -un
    aw(3, 2) = -undemi*x(3)
    aw(3, 3) = -undemi*y(3)
    aw(3, 7) = un
    aw(3, 8) = -undemi*x(3)
    aw(3, 9) = -undemi*y(3)
!
!===================================================================
! --- DETERMINATION DE LA MATRICE AN QUI EST TELLE QUE             =
! --- ALPHA = AN*UN + AM*UM                                        =
! --- SOIT AN = AAI * (AW + AB)                                    =
! --- UN DESIGNE LES DEPLACEMENTS DE FLEXION (W,BETAX,BETAY)       =
! --- UM DESIGNE LES DEPLACEMENTS DE MEMBRANE (UX,UY)              =
!===================================================================
!
    do i = 1, 3
        do k = 1, 3
            do j = 1, 9
                an(i, j) = an(i, j)+aai(i, k)*(aw(k, j)+ab(k, j))
            end do
        end do
    end do
!
!===================================================================
! --- DETERMINATION DE LA MATRICE AM QUI EST TELLE QUE             =
! --- ALPHA = AN*UN + AM*UM                                        =
! --- SOIT AM = AAI*AL                                             =
!===================================================================
!
    do i = 1, 3
        do k = 1, 3
            do j = 1, 6
                am(i, j) = am(i, j)+aai(i, k)*al(k, j)
            end do
        end do
    end do
!
end subroutine
