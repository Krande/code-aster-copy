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

module elg_data_module
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
! person_in_charge: natacha.bereux at edf.fr
!
!
    use aster_petsc_module
    use elg_context_type
    use petsc_data_module
    use saddle_point_context_type
!
    implicit none
!
    private
#include "asterf.h"
#include "asterfort/assert.h"
!
!----------------------------------------------------------------------
!  Données globales utilisées par la fonctionnalité ELIM_LAGR='OUI'
!
!     on prévoit de pouvoir utiliser simultanément ELIM_LAGR='OUI'
!     avec nmax_ctxt matrice(s) Aster différente(s).
    integer(kind=8), parameter :: nmax_ctxt = 5
!     KE est l'indice à utiliser dans le tableau elg_ctxt(:)
!     C'est la variable "sensible" que l'on doit "positionner" avec
!     beaucoup de soins.
!     Aujourd'hui, KE est positionné au début de preres.f
#ifdef ASTER_HAVE_PETSC
    PetscInt, public :: ke
#else
    integer(kind=4), public :: ke
#endif
!     Tableau d'objets de type elim_lagr_ctxt.
!     Chaque objet contient toutes les données nécessaires pour
!     procéder à l'élimination des multiplicateurs de Lagrange.
    type(elim_lagr_ctxt), public, dimension(nmax_ctxt), target :: elg_context
!

! Routine de gestion des données
    public :: elg_gest_data
!
!----------------------------------------------------------------------
! Notations :
!-------------
! On peut résoudre le système dualisé suivant
! en "éliminant" les contraintes A*X=c.
!
!             ! B    A' !   (X)    (b)
!             !         ! *      =
!             ! A    0  !   (L)    (c)
!
!   1) on calcule T = noyau de A
!   2) on calcule X0 solution particulière de : A*X=c,
!       on cherche ensuite X sous la forme X=X0+TY
!   3) On résoud le sytème "réduit" :
!       [T'*B*T]*Y = T'*(b - B*X0)
!   4) On peut alors calculer X=X0 + T*Y
!   5) On peut alors calculer L = (A*A') \ A*(b - B*X)
!   6) La solution complète est : [X, L]
!
!   Dimensions des matrices et vecteurs :
!     n1 (= nphys) : nombre de ddls "physiques"
!     n2 (= nlag)  : nombre de contraintes de A
!                   (= nbre ddls "Lagrange 1" par exemple)
!         normalement :   n2 <= n1
!     n3 = n1-n2
!     B          : n1xn1
!     A          : n2xn1  (on suppose A de rang maximum)
!     T=ker(A)   : n1xn3  (si A est de rang maximum)
!     Kr=T'*B*T  : n3xn3
!     R          : n2xn2
!     X,b,X0     : n1
!     L,c        : n2
!     Fr=T'*(b - B*X0) : n3
!     Y          : n3
!----------------------------------------------------------------------
! Remarque très importante :
!  Certaines matrices (masse, amortissement, ...) ne contiennent pas A
!  Quand on veut les "réduire", il faut utiliser le A de la matrice de
!  rigidité associée.
!  Dans ce cas, on appelle elima1.F avec l'argument RIGI=MATRIG2
!  Dans apelim.F, on ne calcule pas les matrices Ctrans, Tfinal et RCt
!  mais on "pointe" sur celles de la matrice MATRIG2
!
!----------------------------------------------------------------------
! Correspondance entre les variables :
! ------------------------------------
! B      -> MatB
! A'     -> Ctrans
! T      -> Tfinal
! T'*B*T -> Kproj
! X0     -> VX0
! b      -> VecB
! c      -> VecC
!----------------------------------------------------------------------
contains
!
!--------------------------------------------------------------
! BUT :
!   * Gestion des variables globales du module elim_lagr_data_module
!   * Positionner l'indice KE correspondant à une matrice Aster
!
! IN  : ACTION :
!        / 'NOTE' : pour "déclarer" une nouvelle matrice
!                   (appelée par exemple dans preres.f)
!        / 'CHERCHE' : pour positionner KE
!                   (appelée par exemple dans resoud.f)
!        / 'EFFACE' : pour effacer une matrice
!                   (appelée par exemple dans detrsd.f)
! IN  : MAT1 :   / nom de la SD_MATR_ASSE complète
!                / ' ' si action='EFFACE'
! IN  : MAT2  : nom de la SD_MATR_ASSE réduite
! IN  : RIGI1 : / nom de la SD_MATR_ASSE complète qui contient
!                 réellement les relations linéaires
!               / ' '
!               Cet argument ne sert que pour action='NOTE'
!---------------------------------------------------------------
    subroutine elg_gest_data(action, mat1, mat2, rigi1)
!
!   Dummy arguments
        character(len=*), intent(in)  :: action, mat1, mat2, rigi1
#ifdef ASTER_HAVE_PETSC
!   Local variables
        integer(kind=8) :: k, ktrou, iprem, kpos
!
        save iprem
        data iprem/0/
!----------------------------------------------------------------
        iprem = iprem+1
!
        ASSERT(action .eq. 'NOTE' .or. action .eq. 'CHERCHE' .or. action .eq. 'EFFACE')
        if (action .ne. 'NOTE') then
            ASSERT(rigi1 .eq. ' ')
        end if
!
!
!     -- au 1er appel on initialise les données
!     -----------------------------------------
        if (iprem .eq. 1) then
            do k = 1, nmax_ctxt
                elg_context(k) = new_elg_context()
            end do
            ke = 0
        end if
!
        if (action .eq. 'NOTE') then
!       La matrice mat1 a-t-elle déjà été enregistrée ?
            kpos = 1
            do while ((trim(mat1) /= elg_context(kpos)%full_matas) .and. (kpos <= nmax_ctxt))
                kpos = kpos+1
            end do
            if (kpos <= nmax_ctxt) then
                ke = to_petsc_int(kpos)
!          On vérifie que la matrice de rigidité enregistrée a bien le même nom
!          que la matrice rigi1 fournie en entrée de la routine
                ASSERT(elg_context(ke)%k_matas == rigi1)
!          La matrice réduite n'est pas forcément la même
                elg_context(ke)%reduced_matas = mat2
            elseif (kpos > nmax_ctxt) then
!           la matrice n'a jamais été rencontrée
                ktrou = 0
!       -- on cherche une place libre
                do k = 1, nmax_ctxt
                    if (elg_context(k)%full_matas .eq. ' ') then
                        ktrou = k
                        goto 1
                    end if
                end do
1               continue
                ASSERT(ktrou .gt. 0)
                ke = to_petsc_int(ktrou)
                elg_context(ke)%full_matas = mat1
                elg_context(ke)%reduced_matas = mat2
                elg_context(ke)%k_matas = rigi1
            end if
        end if
!
!
        if (action .eq. 'CHERCHE') then
            ktrou = 0
            do k = 1, nmax_ctxt
                if (elg_context(k)%full_matas .eq. mat1) then
                    ktrou = k
                    goto 2
                end if
            end do
2           continue
            ASSERT(ktrou .gt. 0)
            if (mat2 .ne. ' ') then
                ASSERT(elg_context(ktrou)%reduced_matas .eq. mat2)
            end if
            ke = to_petsc_int(ktrou)
        end if
!
        if (action .eq. 'EFFACE') then
            ASSERT(mat1 .eq. ' ')
            ktrou = 0
            do k = 1, nmax_ctxt
                if (elg_context(k)%reduced_matas .eq. mat2) then
                    ktrou = k
                    goto 3
                end if
            end do
3           continue
            if (ktrou .eq. 0) goto 4
!
            call free_elg_context(elg_context(ktrou))
4           continue
        end if
!
#else
        character(len=1) :: kdummy
        kdummy = action//mat1//mat2//rigi1
#endif
    end subroutine elg_gest_data
!
!
end module elg_data_module
