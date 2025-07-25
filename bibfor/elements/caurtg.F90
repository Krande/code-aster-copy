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
subroutine caurtg(nomte, ncmp, sigmau, sigrtg)
!.======================================================================
    implicit none
!
!      CAURTG  -- PASSAGE DES CONTRAINTES DE CAUCHY SIGMAU
!                 CALCULEES DANS LE REPERE UTILISATEUR
!                 VERS LE REPERE UTILISATEUR TOURNE DE
!                 LA ROTATION FAISANT PASSER DE L'ETAT
!                 INITIAL A L'ETAT DEFORME DANS LE CAS GROT_GDEP .
!                 SIGRTG DESIGNE LES CONTRAINTES DE CAUCHY DANS
!                 CE DERNIER REPERE .
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMTE          IN     K16      NOM DU TYPE D'ELEMENT
!    NCMP           IN     I        NOMBRE DE COMPOSANTES DU TENSEUR
!                                   DES CONTRAINTES
!    SIGMAU(NCMP,1) IN     R        VECTEUR DES CONTRAINTES
!                                   DE CAUCHY DANS LE REPERE UTILISATEUR
!    SIGRTG(NCMP,1) VAR    R        VECTEUR DES CONTRAINTES DE CAUCHY
!                                   TOURNE DU REPERE UTILISATEUR VERS
!                                   LE REPERE TRANSFORME DU REPERE
!                                   UTILISATEUR PAR LA GRANDE ROTATION
!                                   CALCULEE POUR COMP_ELAS - GROT_GDEP
!
!
!.========================= DEBUT DES DECLARATIONS ====================
! -----  ARGUMENTS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/marota.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "asterfort/vectan.h"
#include "asterfort/Behaviour_type.h"
    integer(kind=8) :: ncmp
    character(len=16) :: nomte
    real(kind=8) :: sigmau(ncmp, 1), sigrtg(ncmp, 1)
! -----  VARIABLES LOCALES
    real(kind=8) :: vecthe(9, 3), vecta(9, 2, 3)
    real(kind=8) :: vectpt(9, 2, 3), vectn(9, 3)
    real(kind=8) :: xab(3, 3), sigmad(3, 3), sigmat(3, 3)
    real(kind=8) :: drot(3, 3), tetag(3)
!
    aster_logical :: lgreen
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
! --- INITIALISATIONS :
!     ---------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idepl, igeom, ii, in, iret
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: lzi, lzr, nb1, nb2
!-----------------------------------------------------------------------
    lgreen = .false.
!
! --- RECUPERATION DE LA CARTE DE COMPORTEMENT :
!     ----------------------------------------
    call jevech('PCOMPOR', 'L', vk16=compor)
!
    if (compor(DEFO) .eq. 'GROT_GDEP') then
        lgreen = .true.
    end if
!
! --- RECUPERATION DU CHAMP DE DEPLACEMENT DANS LE CAS GROT_GDEP :
!     ---------------------------------------------------------
    if (lgreen) then
        call tecach('OOO', 'PDEPLAR', 'L', iret, iad=idepl)
    else
        goto 999
    end if
!
! --- RECUPERATION DES COORDONNEES DES NOEUDS DANS LA GEOMETRIE
! --- INITIALE :
!     --------
    call jevech('PGEOMER', 'L', igeom)
!
! --- RECUPERATION DES OBJETS INITIALISES :
!     -----------------------------------
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
!
! --- NOMBRE DE NOEUDS (NB1 : SERENDIP, NB2 : LAGRANGE) :
!     -------------------------------------------------
    nb1 = zi(lzi+1-1)
    nb2 = zi(lzi+2-1)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
! --- AFFECTATION DES VECTEURS DE TRANSLATION ET DE ROTATION :
!     ------------------------------------------------------
    do in = 1, nb1
        do ii = 1, 3
            vecthe(in, ii) = zr(idepl+6*(in-1)+ii+3-1)
        end do
    end do
!
    do ii = 1, 3
        vecthe(nb2, ii) = zr(idepl+6*nb1+ii-1)
    end do
!
! --- DETERMINATION DES REPERES LOCAUX AUX NOEUDS DANS LA
! --- CONFIGURATION INITIALE
! --- VECTA DESIGNE LES VECTEURS COVARIANTS DANS LE PLAN MOYEN A
! ---       CHAQUE NOEUD
! --- VECTN DESIGNE LES VECTEURS NORMAUX AU PLAN MOYEN
! --- VECTPT DESIGNE LES REPERES LOCAUX ORTHORNORMES EN CHAQUE
! --- NOEUD DANS LA CONFIGURATION INITIALE :
!     ------------------------------------
    call vectan(nb1, nb2, zr(igeom), zr(lzr), vecta, &
                vectn, vectpt)
!
! ---   MISE DU VECTEUR DES CONTRAINTES DANS LE REPERE UTILISATEUR
! ---   SOUS LA FORME D'UN TENSEUR 3X3 :
!       ------------------------------
    do i = 1, nb2
!
        tetag(1) = vecthe(i, 1)
        tetag(2) = vecthe(i, 2)
        tetag(3) = vecthe(i, 3)
        call marota(tetag, drot)
!
        sigmat(1, 1) = sigmau(1, i)
        sigmat(2, 2) = sigmau(2, i)
        sigmat(3, 3) = sigmau(3, i)
        sigmat(1, 2) = sigmau(4, i)
        sigmat(2, 1) = sigmat(1, 2)
        if (ncmp .eq. 6) then
            sigmat(1, 3) = sigmau(5, i)
            sigmat(2, 3) = sigmau(6, i)
            sigmat(3, 1) = sigmat(1, 3)
            sigmat(3, 2) = sigmat(2, 3)
        end if
!
! ---   ROTATION DU TENSEUR DES CONTRAINTES DE CAUCHY DE LA
! ---   ROTATION FAISANT PASSER DE L'ETAT INITAL A L'ETAT DEFORME :
!       ---------------------------------------------------------
        call utbtab('ZERO', 3, 3, sigmat, drot, &
                    xab, sigmad)
!
! ---   AFFECTATION DU VECTEUR EN SORTIE DES CONTRAINTES
! ---   DE CAUCHY DANS LE REPERE UTILISATEUR TOURNE :
!       -------------------------------------------
        sigrtg(1, i) = sigmad(1, 1)
        sigrtg(2, i) = sigmad(2, 2)
        sigrtg(3, i) = sigmad(3, 3)
        sigrtg(4, i) = sigmad(1, 2)
        if (ncmp .eq. 6) then
            sigrtg(5, i) = sigmad(1, 3)
            sigrtg(6, i) = sigmad(2, 3)
        end if
!
    end do
!
999 continue
!.============================ FIN DE LA ROUTINE ======================
end subroutine
