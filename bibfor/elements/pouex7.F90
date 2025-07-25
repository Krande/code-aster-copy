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
subroutine pouex7(sk, ey, ez)
    implicit none
    real(kind=8) :: sk(105)
    real(kind=8) :: ey, ez
!    -------------------------------------------------------------------
!
!    * CE SOUS PROGRAMME FAIT LE CHANGEMENT DE VARIABLES :
!      POUR LES ELEMENTS DE POUTRES A 7 DDLS PAR NOEUD.
!     (VT,WT) --> (VG,WG) NECESSAIRE POUR LES POUTRES AVEC EXCENTRICITE
!     (CF BATOZ "MODELISATION DES STRUCTURES PAR ELEMENTS FINIS" TOME 2
!          EDITION HERMES 1990  P.181)
!
!    * REMARQUE :
!      LA MATRICE EST STOCKEE TRIANGULAIRE INFERIEURE DANS UN TABLEAU
!      UNICOLONNE
!    * ORDRE SUPPOSE DES DDLS :
!      DX,DY,DZ,DRX,DRY,DRZ,GRX,  DX,DY,...,GRX
!    -------------------------------------------------------------------
!  DONNEES NON MODIFIEES
!
! IN TYPE ! NOM    ! TABLEAU !             SIGNIFICATION
! IN -------------------------------------------------------------------
! IN R*8  ! EY     !     -   ! COMPOSANTE GT SUR Y PRINCIPAL
! IN R*8  ! EZ     !     -   ! COMPOSANTE GT SUR Z PRINCIPAL
!
! VAR TYPE ! NOM   ! TABLEAU !             SIGNIFICATION
! VAR ------------------------------------------------------------------
! VAR R*8 !   SK   !  105    ! MATRICE ELEMENTAIRE UNICOLONNE
!
!
! LOC TYPE !  NOM  ! TABLEAU !              SIGNIFICATION
! LOC ------------------------------------------------------------------
! LOC I   ! IP     !   14    ! POINTEUR SUR L'ELEMENT DIAGONAL PRECEDENT
! LOC R*8 ! SKP    !   105   ! MATRICE DE TRAVAIL
!     ------------------------------------------------------------------
    real(kind=8) :: skp(105)
    integer(kind=8) :: ip(14)
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i
!-----------------------------------------------------------------------
    data ip/0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 66, 78, 91/
! ---------------------------------------------------------------------
!
!
    if (ez .eq. 0.0d0 .and. ey .eq. 0.0d0) goto 999
!
    do i = 1, 105
        skp(i) = sk(i)
    end do
!
!
!     --LES INSTRUCTIONS SUIVANTES ONT ETE OBTENUES PAR MATHEMATICA
!       (ON NE SUPPOSE AUCUN TERME NUL DANS LA MATRICE SK(*))
!
    skp(01+ip(04)) = sk(01+ip(04))+sk(01+ip(03))*ey-sk(01+ip(02))*ez
    skp(01+ip(11)) = sk(01+ip(11))+sk(01+ip(10))*ey-sk(01+ip(09))*ez
    skp(02+ip(04)) = sk(02+ip(04))+sk(02+ip(03))*ey-sk(02+ip(02))*ez
    skp(02+ip(11)) = sk(02+ip(11))+sk(02+ip(10))*ey-sk(02+ip(09))*ez
    skp(03+ip(04)) = sk(03+ip(04))+sk(03+ip(03))*ey-sk(02+ip(03))*ez
    skp(03+ip(11)) = sk(03+ip(11))+sk(03+ip(10))*ey-sk(03+ip(09))*ez
    skp(04+ip(04)) = sk( &
                     04+ip(04))+sk(03+ip(04))*ey-sk(02+ip(04))*ez-ez*(sk(02+ip(04))+sk(02&
                     &+ip(03))*ey-sk(02+ip(02))*ez)+ey*(sk(03+ip(04))+sk(03+ip(03))*ey-sk(02+&
                     &ip(03))*ez &
                     )
    skp(04+ip(05)) = sk(04+ip(05))+sk(03+ip(05))*ey-sk(02+ip(05))*ez
    skp(04+ip(06)) = sk(04+ip(06))+sk(03+ip(06))*ey-sk(02+ip(06))*ez
    skp(04+ip(07)) = sk(04+ip(07))+sk(03+ip(07))*ey-sk(02+ip(07))*ez
    skp(04+ip(08)) = sk(04+ip(08))+sk(03+ip(08))*ey-sk(02+ip(08))*ez
    skp(04+ip(09)) = sk(04+ip(09))+sk(03+ip(09))*ey-sk(02+ip(09))*ez
    skp(04+ip(10)) = sk(04+ip(10))+sk(03+ip(10))*ey-sk(02+ip(10))*ez
    skp(04+ip(11)) = sk( &
                     04+ip(11))+sk(03+ip(11))*ey-sk(02+ip(11))*ez-ez*(sk(04+ip(09))+sk(03&
                     &+ip(09))*ey-sk(02+ip(09))*ez)+ey*(sk(04+ip(10))+sk(03+ip(10))*ey-sk(02+&
                     &ip(10))*ez &
                     )
    skp(04+ip(12)) = sk(04+ip(12))+sk(03+ip(12))*ey-sk(02+ip(12))*ez
    skp(04+ip(13)) = sk(04+ip(13))+sk(03+ip(13))*ey-sk(02+ip(13))*ez
    skp(04+ip(14)) = sk(04+ip(14))+sk(03+ip(14))*ey-sk(02+ip(14))*ez
    skp(05+ip(11)) = sk(05+ip(11))+sk(05+ip(10))*ey-sk(05+ip(09))*ez
    skp(06+ip(11)) = sk(06+ip(11))+sk(06+ip(10))*ey-sk(06+ip(09))*ez
    skp(07+ip(11)) = sk(07+ip(11))+sk(07+ip(10))*ey-sk(07+ip(09))*ez
    skp(08+ip(11)) = sk(08+ip(11))+sk(08+ip(10))*ey-sk(08+ip(09))*ez
    skp(09+ip(11)) = sk(09+ip(11))+sk(09+ip(10))*ey-sk(09+ip(09))*ez
    skp(10+ip(11)) = sk(10+ip(11))+sk(10+ip(10))*ey-sk(09+ip(10))*ez
    skp(11+ip(11)) = sk( &
                     11+ip(11))+sk(10+ip(11))*ey-sk(09+ip(11))*ez-ez*(sk(09+ip(11))+sk(09&
                     &+ip(10))*ey-sk(09+ip(09))*ez)+ey*(sk(10+ip(11))+sk(10+ip(10))*ey-sk(09+&
                     &ip(10))*ez &
                     )
    skp(11+ip(12)) = sk(11+ip(12))+sk(10+ip(12))*ey-sk(09+ip(12))*ez
    skp(11+ip(13)) = sk(11+ip(13))+sk(10+ip(13))*ey-sk(09+ip(13))*ez
    skp(11+ip(14)) = sk(11+ip(14))+sk(10+ip(14))*ey-sk(09+ip(14))*ez
!
    do i = 1, 105
        sk(i) = skp(i)
    end do
!
!
!
999 continue
end subroutine
