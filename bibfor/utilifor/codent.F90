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
subroutine codent(entier, cadre, chaine, kstop)
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
    integer(kind=8), intent(in) :: entier
    character(len=*), intent(in) :: cadre
    character(len=*), intent(out) :: chaine
    character(len=*), optional :: kstop
!
!     ------------------------------------------------------------------
!     CODAGE D'UN ENTIER SUR UNE CHAINE DE CARACTERE
!     ------------------------------------------------------------------
! IN  KSTOP  : CH(*) : ' '/'F' (VOIR CI-DESSOUS)
! IN  ENTIER : IS    : ENTIER A CONVERTIR EN CHAINE
! IN  CADRE  : CH(*) : TYPE DE CADRAGE
!          D     CADRAGE A DROITE
!          D0    CADRAGE A DROITE ET ON COMPLETE A GAUCHE PAR DES ZERO
!          G     CADRAGE A GAUCHE
! OUT CHAINE : CH(*) : CHAINE RECEPTACLE, ON UTILISE TOUTE LA LONGUEUR
!                      DE LA CHAINE
!     ------------------------------------------------------------------
!     REMARQUE :
!         SI KSTOP=' ' : EN CAS D'ERREUR (A LA TRANSCRIPTION OU DANS LE TYPE DE
!                        CADRAGE)   LA CHAINE EST REMPLIE D'ETOILE
!         SI KSTOP='F' : EN CAS D'ERREUR, ON ARRETE LE CODE
!     ROUTINE(S) UTILISEE(S) :
!         -
!     ROUTINE(S) FORTRAN     :
!         LEN    MOD
!     ------------------------------------------------------------------
! FIN CODENT
!     ------------------------------------------------------------------
!
!
    integer(kind=8) :: lg, ent, ival
    aster_logical :: neg
    character(len=1) :: chiffr(0:9)
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ier, il, il1
!-----------------------------------------------------------------------
    data chiffr/'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/
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
    ent = entier
    neg = ent .lt. 0
    if (neg) ent = -ent
    lg = len(chaine)
!
!     ON CADRE A DOITE A PRIORI   CADRAGE A DROITE
    il = lg+1
10  continue
    il = il-1
    if (il .le. 0) then
        ier = 1
        goto 99000
    else
        ival = mod(ent, 10)
        chaine(il:il) = chiffr(ival)
        ent = ent/10
    end if
    if (ent .ne. 0) goto 10
!
    if (neg) then
        il = il-1
        if (il .le. 0) then
            ier = 1
            goto 99000
        else
            chaine(il:il) = '-'
        end if
    end if
!
    if (cadre(1:1) .eq. 'D') then
!        --- CADRAGE A DROITE ---
        if (len(cadre) .gt. 1) then
            if (cadre(2:2) .eq. '0') then
                if (neg) chaine(il:il) = '0'
                do i = il-1, 1, -1
                    chaine(i:i) = '0'
                end do
                if (neg) chaine(1:1) = '-'
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
