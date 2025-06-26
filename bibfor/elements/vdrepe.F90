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
subroutine vdrepe(nomtez, matevn, matevg)
!.======================================================================
    implicit none
!
!      VDREPE   -- DETERMINATION DES MATRICES DE PASSAGE
!                  DES REPERES INTRINSEQUES AUX NOEUDS  DE L'ELEMENT
!                  AU REPERE UTILISATEUR (MATRICE MATEVN)
!                  ET DES REPERES INTRINSEQUES AUX POINTS D'INTEGRATION
!                  DE L'ELEMENT AU REPERE UTILISATEUR (MATRICE MATEVG)
!                  POUR LES ELEMENTS DE COQUE EPAISSE 3D .
!
!   ARGUMENT        E/S   TYPE         ROLE
!    NOMTE          IN     K*       NOM DU TYPE D'ELEMENT
!    MATEVN(2,2,10) OUT    R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX NOEUDS  DE
!                                   L'ELEMENT AU REPERE UTILISATEUR
!    MATEVG(2,2,10) OUT    R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX POINTS
!                                   D'INTEGRATION DE L'ELEMENT AU
!                                   REPERE UTILISATEUR
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "jeveux.h"
#include "asterc/r8dgrd.h"
#include "asterfort/coqrep.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
    character(len=*) :: nomtez
    real(kind=8) :: matevn(2, 2, 1), matevg(2, 2, 1)
! -----  VARIABLES LOCALES
    character(len=16) :: nomte
    real(kind=8) :: pgl(3, 3), r8bid4(4)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     ---------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idec, igau, ino, j, jcoqu, k
    integer(kind=8) :: lzi, lzr, nb2, npgsr
    real(kind=8) :: alpha, beta, c
    real(kind=8) :: s
!-----------------------------------------------------------------------
    nomte = nomtez
!
! --- RECUPERATION DES OBJETS DESCRIPTEURS DES ELEMENTS :
!     =================================================
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
!
! --- NOMBRE DE NOEUDS DE L'ELEMENT  :
!     -----------------------------
    nb2 = zi(lzi-1+2)
!
! --- NOMBRE DE POINTS D'INTEGRATION DE L'ELEMENT (SOUS-INTEGRE) :
!     ----------------------------------------------------------
    npgsr = zi(lzi-1+3)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
! --- RECUPERATION DES ANGLES DETERMINANT LE REPERE UTILISATEUR
! --- PAR RAPPORT AU REPERE GLOBAL :
!     ============================
    call jevech('PCACOQU', 'L', jcoqu)
!
    alpha = zr(jcoqu+1)*r8dgrd()
    beta = zr(jcoqu+2)*r8dgrd()
!
! --- DETERMINATION DES MATRICES DE PASSAGE DES REPERES INTRINSEQUES
! --- AUX NOEUDS DE L'ELEMENT AU REPERE UTILISATEUR :
!     =============================================
!
! --- ADRESSE DES MATRICES DE PASSAGE DU REPERE GLOBAL AUX REPERES
! --- INTRINSEQUES AUX NOEUDS DE L'ELEMENT DANS LE TABLEAU .DESR :
!     ----------------------------------------------------------
    idec = 1090
!
! --- BOUCLE SUR LES NOEUDS DE L'ELEMENT :
!     ----------------------------------
    do ino = 1, nb2
!
! ---   RECUPERATION DE LA MATRICE DE PASSAGE AU NOEUD COURANT :
!       ------------------------------------------------------
        k = 0
        do j = 1, 3
            do i = 1, 3
                k = k+1
                pgl(i, j) = zr(lzr+idec+(ino-1)*9+k-1)
            end do
        end do
!
! ---   DETERMINATION DE LA PROJECTION DU VECTEUR X DU REPERE
! ---   UTILISATEUR SUR LE FEUILLET TANGENT A LA COQUE AU NOEUD
! ---   COURANT :
!       -------
!
        call coqrep(pgl, alpha, beta, r8bid4, r8bid4, &
                    c, s)
!
!       -- (C,S) N'EST PAS TOUJOURS EXACTEMENT DE NORME=1:
        c = c/sqrt(c*c+s*s)
        s = s/sqrt(c*c+s*s)
!
        matevn(1, 1, ino) = c
        matevn(2, 1, ino) = s
        matevn(1, 2, ino) = -s
        matevn(2, 2, ino) = c
!
    end do
!
! --- DETERMINATION DES MATRICES DE PASSAGE DES REPERES INTRINSEQUES
! --- AUX POINTS D'INTEGRATION DE L'ELEMENT AU REPERE UTILISATEUR :
!     ===========================================================
!
! --- ADRESSE DES MATRICES DE PASSAGE DU REPERE GLOBAL AUX REPERES
! --- INTRINSEQUES AUX POINTS D'INTEGRATION DE L'ELEMENT
! --- DANS LE TABLEAU .DESR :
!     ---------------------
    idec = 2000
!
! --- BOUCLE SUR LES POINTS D'INTEGRATION DE L'ELEMENT (SOUS-INTEGRE) :
!     --------------------------------------------------------------
    do igau = 1, npgsr
!
! ---   RECUPERATION DE LA MATRICE DE PASSAGE AU POINT D'INTEGRATION
! ---   COURANT :
!       -------
        k = 0
        do j = 1, 3
            do i = 1, 3
                k = k+1
                pgl(i, j) = zr(lzr+idec+(igau-1)*9+k-1)
            end do
        end do
!
! ---   DETERMINATION DE LA PROJECTION DU VECTEUR X DU REPERE
! ---   UTILISATEUR SUR LE FEUILLET TANGENT A LA COQUE AU POINT
! ---   D'INTEGRATION COURANT :
!       ---------------------
        call coqrep(pgl, alpha, beta, r8bid4, r8bid4, &
                    c, s)
!
!       -- (C,S) N'EST PAS TOUJOURS EXACTEMENT DE NORME=1:
        c = c/sqrt(c*c+s*s)
        s = s/sqrt(c*c+s*s)
!
        matevg(1, 1, igau) = c
        matevg(2, 1, igau) = s
        matevg(1, 2, igau) = -s
        matevg(2, 2, igau) = c
!
    end do
!
end subroutine
