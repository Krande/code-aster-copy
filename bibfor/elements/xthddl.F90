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
subroutine xthddl(nfh, nddlno, nno, stano, option, &
                  nomte, mat, vect)
! person_in_charge: sam.cuvilliez at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/indent.h"
#include "asterfort/teattr.h"
    integer(kind=8), intent(in) :: nfh, nddlno, nno, stano(*)
    character(len=16), intent(in) :: option, nomte
    real(kind=8), optional, intent(inout) :: mat(*)
    real(kind=8), optional, intent(out) :: vect(*)
!
!     BUT: THERMIQUE + ELEMENTS X-FEM LINEAIRES, SUPPRIMER LES DDL "EN
!          TROP" (ATTENTION STOCKAGE SYMETRIQUE POUR LA MATRICE "MAT")
!          ROUTINE EQUIVALENTE EN MECA -> XTEDLL
!
! IN       NFH    : NOMBRE DE FONCTIONS HEAVYSIDE
! IN       NDDLNO : NOMBRE DE DDL PAR NOEUD
! IN       NNO    : NOMBRE DE NOEUDS
! IN       STANO  : STATUT DES NOEUDS
! IN       OPTION : OPTION DE CALCUL DU TE
! IN       NOMTE  : NOM DU TYPE ELEMENT
!
! IN/OUT   MAT    : MATRICE DE RIGIDITE OU DE MASSE
! OUT      VECT   : VECTEUR SECOND MEMBRE
!
!-----------------------------------------------------------------------
!---------------- DECLARATION DES VARIABLES LOCALES  -------------------
!
    integer(kind=8) :: ier, istatu, ino, i, j, ielim, in, ddlmax
    integer(kind=8) :: nddl
!      AU PLUS 8*3=24 DDL (MAX ATTEINT POUR L'HEXA8 XHT)
    parameter(ddlmax=24)
    integer(kind=8) :: posddl(ddlmax)
    character(len=8) :: tyenel
    aster_logical :: lelim, lmat, lvec
    real(kind=8) :: dmax, dmin, codia
!
!-------------------------------------------------------------
!
!
!-------------------------------------------------------------
!   NOMS DES OPTIONS AUTORISEES
!-------------------------------------------------------------
!
    lmat = .false.
    lvec = .false.
!
!   OPTIONS RELATIVES A UNE MATRICE
    if (option .eq. 'RIGI_THER' .or. option .eq. 'RIGI_THER_PARO_F' .or. option .eq. &
        'RIGI_THER_PARO_R' .or. option .eq. 'MASS_THER') then
        lmat = .true.
!   OPTIONS RELATIVES A UN VECTEUR
    elseif (option .eq. 'CHAR_THER_EVOL' &
            .or. option .eq. 'CHAR_THER_PARO_F' &
            .or. option .eq. 'CHAR_THER_PARO_R') then
        lvec = .true.
    else
        ASSERT(.false.)
    end if
!
!-------------------------------------------------------------
!   VERIFICATION DE LA COHERENCE OPTION / ARGUMENTS OPTIONNELS
!-------------------------------------------------------------
!
    if (present(mat) .and. .not. present(vect)) then
        ASSERT(lmat .and. .not. lvec)
    else if (.not. present(mat) .and. present(vect)) then
        ASSERT(.not. lmat .and. lvec)
!   EXACTEMENT UN DES 2 ARGUMENTS mat OU vect EST OBLIGATOIRE
    else
        ASSERT(.false.)
    end if
!
!-------------------------------------------------------------
!
! --- TYPE D'ENRICHISSEMENT DE L'ELEMENT ET TYPE D'ELIMINATION
!
    call teattr('S', 'XFEM', tyenel, ier, typel=nomte)
    if (tyenel(1:2) .eq. 'XH') ielim = 1
    if (tyenel(1:2) .eq. 'XT') ielim = 2
    if (tyenel(1:3) .eq. 'XHT') ielim = 3
!
!     REMPLISSAGE DU VECTEUR POS : POSITION DES DDLS A SUPPRIMER
!
    nddl = nddlno*nno
    ASSERT(nddl .le. ddlmax)
    do ino = 1, ddlmax
        posddl(ino) = 0
    end do
!
!     VRAI SI ON ELIMINE LES DDLS D'AU MOINS UN NOEUD
    lelim = .false.
!
    do ino = 1, nno
!
        call indent(ino, nddlno, 0, nno, in)
!
        if (ielim .eq. 1) then
!         1) CAS DES MAILLES 'ROND'
!         -------------------------
!         STATUT DES NOEUDS ENRICHIS
            istatu = abs(stano(ino))
            ASSERT(istatu .le. 1)
            if (istatu .eq. 0) then
!           ON SUPPRIME LES DDL H
                posddl(in+1+1) = 1
                lelim = .true.
            end if
!
        else if (ielim .eq. 2) then
!         2) CAS DES MAILLES 'CARRÉ'
!         --------------------------
!         STATUT DES NOEUDS ENRICHIS
            istatu = abs(stano(ino))
            ASSERT(istatu .le. 2 .and. istatu .ne. 1)
            if (istatu .eq. 0) then
!           ON SUPPRIME LES DDL E
                posddl(in+1+nfh+1) = 1
                lelim = .true.
            end if
!
        else if (ielim .eq. 3) then
!         3) CAS DES MAILLES 'ROND-CARRÉ'
!         ------------------------------
!         STATUT DES NOEUDS ENRICHIS
            istatu = abs(stano(ino))
            ASSERT(istatu .le. 3)
            if (istatu .eq. 2) then
!           ON SUPPRIME LES DDL H
                posddl(in+1+1) = 1
                lelim = .true.
            else if (istatu .eq. 1) then
!           ON SUPPRIME LES DDL E
                posddl(in+1+nfh+1) = 1
                lelim = .true.
            else if (istatu .eq. 0) then
!           ON SUPPRIME LES DDL H ET E
                posddl(in+1+1) = 1
                posddl(in+1+nfh+1) = 1
                lelim = .true.
            end if
!
        end if
!
    end do
!
    if (lelim) then
!
!     POUR LES OPTIONS CONCERNANT DES MATRICES :
!       CALCUL DU COEFFICIENT DIAGONAL POUR
!       L'ELIMINATION DES DDLS HEAVISIDE
        if (lmat) then
            dmin = r8maem()
            dmax = -r8maem()
            do i = 1, nddl
                codia = mat((i-1)*i/2+i)
                if (codia .gt. dmax) then
                    dmax = codia
                else if (codia .lt. dmin) then
                    dmin = codia
                end if
            end do
            codia = (dmax+dmin)/2.0d0
            if (codia .eq. 0) codia = 1
        end if
!
        do i = 1, nddl
            if (posddl(i) .eq. 0) goto 200
!         POUR LES OPTIONS CONCERNANT DES MATRICES :
!           MISE A ZERO DES TERMES HORS DIAGONAUX (I,J)
!           ET MISE A UN DES TERMES DIAGONAUX (I,I)
!           (ATTENTION AU STOCKAGE SYMETRIQUE)
            if (lmat) then
                do j = 1, nddl
                    if (j .lt. i) mat((i-1)*i/2+j) = 0.d0
                    if (j .eq. i) mat((i-1)*i/2+j) = codia
                    if (j .gt. i) mat((j-1)*j/2+i) = 0.d0
                end do
            end if
!         POUR LES OPTIONS CONCERNANT DES VECTEURS :
!           MISE A ZERO DES TERMES I
            if (lvec) vect(i) = 0.d0
200         continue
        end do
!
    end if
!
end subroutine
