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
subroutine codlet(entier, cadre, chaine, kstop)
    implicit none
#include "asterfort/assert.h"
    integer(kind=8), intent(in) :: entier
    character(len=*), intent(in) :: cadre
    character(len=*), intent(out) :: chaine
    character(len=*), optional :: kstop
!
!   ------------------------------------------------------------------
!   CODAGE D'UN ENTIER EN BASE 36 DANS UNE CHAINE DE CARACTERE
!   ------------------------------------------------------------------
! IN  ENTIER : IS    : ENTIER A CONVERTIR EN CHAINE
! IN  CADRE  : CH(*) : TYPE DE CADRAGE
!          D     CADRAGE A DROITE
!          D0    CADRAGE A DROITE ET ON COMPLETE A GAUCHE PAR DES ZERO
!          G     CADRAGE A GAUCHE
! OUT CHAINE : CH(*) : CHAINE RECEPTACLE, ON UTILISE TOUTE LA LONGUEUR
!                      DE LA CHAINE
!     ------------------------------------------------------------------
!     REMARQUES :
!      - EN CAS D'ERREUR (A LA TRANSCRIPTION OU DANS LE TYPE DE CADRAGE)
!        LA CHAINE EST REMPLIE D'ETOILE
!      - POUR LES ENTIERS NEGATIFS ==> VALEUR ABSOLUE
!     ------------------------------------------------------------------
!     ROUTINE(S) UTILISEE(S) :
!         -
!     ROUTINE(S) FORTRAN     :
!         LEN    MOD
!     ------------------------------------------------------------------
!
    integer(kind=8) :: lg, ent, ival, base, basmax, il1, il, ier, i
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    parameter(basmax=36, base=36)
    character(len=1) :: chiffr(0:basmax-1)
    data chiffr/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',&
     &                   'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',&
     &                   'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',&
     &                   'U', 'V', 'W', 'X', 'Y', 'Z'/
!
!
    if (present(kstop)) then
        ASSERT(kstop .eq. ' ' .or. kstop .eq. 'F')
    else
        kstop = 'F'
    end if
!
    ier = 0
    chaine = ' '
!
    ent = abs(entier)
    lg = len(chaine)
!
!     ON CADRE A DROITE A PRIORI   CADRAGE A DROITE
    il = lg+1
10  continue
    il = il-1
    if (il .le. 0) then
        ier = 1
        goto 99000
    else
        ival = mod(ent, base)
        chaine(il:il) = chiffr(ival)
        ent = ent/base
    end if
    if (ent .ne. 0) goto 10
!
!
    if (cadre(1:1) .eq. 'D') then
!        --- CADRAGE A DROITE ---
        if (len(cadre) .gt. 1) then
            if (cadre(2:2) .eq. '0') then
                do i = il-1, 1, -1
                    chaine(i:i) = '0'
                end do
            end if
        end if
!
    else if (cadre(1:1) .eq. 'G') then
!        --- CADRAGE A GAUCHE ---
        il1 = il-1
        do i = 1, lg-il1
            chaine(i:i) = chaine(i+il1:i+il1)
        end do
        chaine(lg-il1+1:) = ' '
    else
        ier = 1
    end if
!
!     SORTIE -----------------------------------------------------------
99000 continue
    if (ier .ne. 0) then
        if (kstop .eq. ' ') then
            do i = 1, lg
                chaine(i:i) = '*'
            end do
        else
            ASSERT(.false.)
        end if
    end if
!
end subroutine
