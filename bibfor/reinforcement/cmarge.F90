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

subroutine cmarge(N, M, Ned, Med, Ned_g, Med_g, margin, Nrd, Mrd, load_rate, c0_c, c0_crd)
    !_____________________________________________________________________
    !
    !      cmarge
    !
    !      CALCUL DE LA MARGE MECANIQUE (PAR FACETTE) POUR VERI_FERRAILLAGE
    !
    !      I N          Efforts N du diagramme d'interaction
    !      I M          Efforts M du diagramme d'interaction
    !      I Ned        Efforts N resultant du chargement etudié
    !      I Med        Efforts M resultant du chargement etudié
    !      I Ned_g      Efforts N resultant du chargement de référence
    !      I Med_g      Efforts M resultant du chargement de référence
    !      O margin     marge mecanique
    !      O Nrd,Mrd    Efforts N,M : projection de la droite des efforts sur le diagramme NM
    !      O load_rate  Taux d utilisation : 1-marge
    !      O c0_c       Distance C0C entre le point de ref. et le point d etude
    !      O c0_crd     Distance C0CRD entre le point de ref. et le point du diagramme

    implicit none
#include "asterfort/utmess.h"
#include "asterc/r8prem.h"
    real(kind=8), pointer, intent(in) :: N(:), M(:)
    real(kind=8), intent(in) :: Ned, Med, Ned_g, Med_g
    real(kind=8), intent(out) :: margin
    real(kind=8) :: Nrd, Mrd
    real(kind=8) :: load_rate, arg_origin_load, arg_origin_curvei
    real(kind=8) :: c0_c
    real(kind=8) :: c0_crd
    real(kind=8), allocatable :: absolute_val_array(:)
    real(kind=8) :: dot_product
    integer(kind=8) :: i, sn, smallest_product, j
    LOGICAL :: inside
    ! number of points of the diagram
    sn = size(N)
    ! Vérifier si le point de reference est à l'intérieur du diagramme
    inside = .FALSE.
    j = sn
    do i = 1, sn
        if ((abs(N(i)-Ned_g) .le. r8prem()*max(abs(N(i)), abs(Ned_g))) .and. &
            (abs(M(i)-Med_g) .le. r8prem()*max(abs(N(i)), abs(Ned_g)))) then
            inside = .FALSE.
            goto 888
        end if
        if ((M(i) > Med_g) .neqv. (M(j) > Med_g)) then
            if (Ned_g < (N(j)-N(i))*(Med_g-M(i))/(M(j)-M(i))+N(i)) then
                inside = .NOT. inside
            end if
        end if
        j = i
    end do
888 continue
    ! error : the reference point is outside the diagram
    if (inside .eqv. .FALSE.) then
        margin = 2.0
        load_rate = -1.0
        c0_c = sqrt((Ned-Ned_g)**2+(Med-Med_g)**2)
        Nrd = -1.d0
        Mrd = -1.d0
        c0_crd = -1.d0
        call utmess('A', 'VERIFERRAILLAGE_11')
        goto 998
    end if

    if ((abs(Ned-Ned_g) .le. r8prem()*max(abs(Ned), abs(Ned_g))) .and. &
        (abs(Med-Med_g) .le. r8prem()*max(abs(Med), abs(Med_g)))) then
        margin = 1.0
        load_rate = 0.d0
        c0_c = sqrt((Ned-Ned_g)**2+(Med-Med_g)**2)
        Nrd = -1.d0
        Mrd = -1.d0
        c0_crd = -1.d0
        call utmess('A', 'VERIFERRAILLAGE_18')
    else

        ! Slope origin_load
        if (abs(Ned-Ned_g) .gt. r8prem()*max(abs(Ned), abs(Ned_g))) then
            arg_origin_load = (Med-Med_g)/(Ned-Ned_g)
        else
            arg_origin_load = (Med-Med_g)/(r8prem()*max(abs(Med), abs(Med_g)))
        end if
        allocate (absolute_val_array(sn))

        !Initializing of Absolute_val_array's elements to a large value
        absolute_val_array = 1.0E20
        ! Loop through each value of N and M
        do i = 1, sn
            !to test the directions of the two vectors
            dot_product = (Ned-Ned_g)*(N(i)-Ned_g)+(Med-Med_g)*(M(i)-Med_g)
            ! If the vectors are in the same direction (positive dot product)
            if (dot_product > 0.0d0) then
                ! Calculation of the difference betwen the two slopes
                if (abs(N(i)-Ned_g) > r8prem()*max(abs(N(i)), abs(Ned_g))) then
                    arg_origin_curvei = (M(i)-Med_g)/(N(i)-Ned_g)
                    absolute_val_array(i) = abs(arg_origin_curvei-arg_origin_load)
                else
                    arg_origin_curvei = (M(i)-Med_g)/r8prem()*max(abs(M(i)), abs(Med_g))
                    absolute_val_array(i) = abs(arg_origin_curvei-arg_origin_load)
                end if
            end if
        end do

        smallest_product = minloc(absolute_val_array, dim=1)
        Nrd = N(smallest_product)
        Mrd = M(smallest_product)
        ! Calculation of margin using the norm (distance between points)
        c0_c = sqrt((Ned-Ned_g)**2+(Med-Med_g)**2)
        c0_crd = sqrt((Nrd-Ned_g)**2+(Mrd-Med_g)**2)
        load_rate = c0_c/c0_crd
        margin = 1-load_rate
        deallocate (absolute_val_array)

    end if
998 continue
end subroutine cmarge
