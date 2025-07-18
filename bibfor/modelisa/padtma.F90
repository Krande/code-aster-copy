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
subroutine padtma(coor1, coor2, nbnott, icoupl, dmin, &
                  no_err)
    implicit none
#include "jeveux.h"
#include "asterfort/padist.h"
#include "asterfort/utmess.h"
    real(kind=8) :: coor1(*), coor2(*), d
    integer(kind=8) :: icoupl(*), nbnott(3)
    integer(kind=8), intent(out), optional :: no_err
!     BUT: CALCULER LA DISTANCE ENTRE 2 MAILLES ET DONNER LE TABLEAU
!          D'INDIRECTION ENTRE LES NUM.DES NOEUDS QUI DESCRIT LE
!          VIS A VIS ENTRE LES 2 MAILLES REALISANT CETTE DISTANCE MINI.
!                  NBNOTT(1)
!  DMIN = MIN     ( +    DIST(N1(I),N2(PERM(I)) )
!         PERM E P I=1
!     P= ENSEMBLE DES PERMUTATIONS DE (1,2,...,NBNO) ENGENDREES PAR
!        LES COUPLAGES POSSIBLES(D'APRES LEUR ORIENTATION)DES 2 MAILLES
! ARGUMENTS
! IN   COOR1  R(*): COORDONNEES DES NBNO 1ERS NOEUDS DE LA MAILLE 1
!                   POUR INO = 1,NBNO
!                   COOR1(3*(INO-1)+1)= X1(NO(INO,MA1))
!                   COOR1(3*(INO-1)+2)= X2(NO(INO,MA1))
!                   COOR1(3*(INO-1)+3)= X3(NO(INO,MA1)) ( EN 2D 0)
! IN   COOR2  R(*): COORDONNEES DES NBNO 1ERS NOEUDS DE LA MAILLE 2
!                   POUR INO = 1,NBNO
!                   COOR2(3*(INO-1)+1)= X1(NO(INO,MA2))
!                   COOR2(3*(INO-1)+2)= X2(NO(INO,MA2))
!                   COOR2(3*(INO-1)+3)= X3(NO(INO,MA2)) ( EN 2D 0)
! IN   NBNOTT I   : (1) NOMBRE DE NOEUDS SOMMETS DE LA MAILLE A EXTRAIRE
!                   (2) NOMBRE DE NOEUDS SUR LES ARRETES
!                   (3) NOMBRE DE NOEUDS INTERIEURS
! OUT  ICOUPL I(*): PERMUTATION QUI REALISE LE MIN COMPLETEE JUSQU'AU
!                   NOMBRE TOTAL DE NOEUDS<=> VIS A VIS ENTRE LES 2
!                   MAILLES
! OUT  DMIN   R   : DISTANCE ENTRE LES 2 MAILLES
!
! OUT NO_ERR    :  present no error is thrown when encoutering perpendicular
!
    real(kind=8) :: x1(3), xn1(3), xn2(3), x2(3), x3(3), x4(3)
    integer(kind=8) :: iperm(106)
! --- DATA DES PERMUTATIONS REPRESENTANT LES VIS A VIS
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ideb, j, k, kdeb, kdeb0, kdeb1
    integer(kind=8) :: kdeb2, n, nbno, nbperm, nbsom, nno
    real(kind=8) :: dmin, s
!-----------------------------------------------------------------------
    data iperm/1, 2, 3, 2, 1, 3,&
     &            1, 2, 3, 4, 5, 6, 2, 3, 1, 5, 6, 4, 3, 1, 2, 6, 4, 5,&
     &            1, 3, 2, 6, 5, 4, 3, 2, 1, 5, 4, 6, 2, 1, 3, 4, 6, 5,&
     &            1, 2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 1, 6, 7, 8, 5,&
     &            3, 4, 1, 2, 7, 8, 5, 6, 4, 1, 2, 3, 8, 5, 6, 7,&
     &            2, 1, 4, 3, 5, 8, 7, 6, 3, 2, 1, 4, 6, 5, 8, 7,&
     &            4, 3, 2, 1, 7, 6, 5, 8, 1, 4, 3, 2, 8, 7, 6, 5/
! --- DEBUT
! --- ORIENTATION DES MAILLES
    no_err = 0
    nbsom = nbnott(1)
    nbno = nbnott(1)+nbnott(2)+nbnott(3)
    if (nbsom .eq. 2) then
        kdeb1 = 0
        kdeb2 = 3
        nno = 3
        do i = 1, 3
            xn1(i) = coor1(3+i)-coor1(i)
            xn2(i) = coor2(3+i)-coor2(i)
        end do
    else if (nbsom .le. 4) then
        if (nbsom .eq. 3) then
            kdeb1 = 6
            kdeb2 = 24
            nno = 6
        else
            kdeb1 = 42
            kdeb2 = 74
            nno = 8
        end if
        do i = 1, 3
            x1(i) = coor1(3+i)-coor1(i)
            x2(i) = coor1(6+i)-coor1(3+i)
            x3(i) = coor2(3+i)-coor2(i)
            x4(i) = coor2(6+i)-coor2(3+i)
        end do
        xn1(1) = x1(2)*x2(3)-x1(3)*x2(2)
        xn1(2) = x1(3)*x2(1)-x1(1)*x2(3)
        xn1(3) = x1(1)*x2(2)-x1(2)*x2(1)
        xn2(1) = x3(2)*x4(3)-x3(3)*x4(2)
        xn2(2) = x3(3)*x4(1)-x3(1)*x4(3)
        xn2(3) = x3(1)*x4(2)-x3(2)*x4(1)
    else
        call utmess('F', 'MODELISA6_7')
    end if
    s = 0.d0
    do i = 1, 3
        s = s+xn1(i)*xn2(i)
    end do
    if (s .gt. 0) then
        kdeb0 = kdeb1
    else if (s .lt. 0) then
        kdeb0 = kdeb2
    else
        if (present(no_err)) then
            no_err = 1
            GO TO 8
        else
            call utmess('F', 'MODELISA6_8')
        end if
    end if
    dmin = 99999999.d0
    nbperm = nbsom
    if (nbsom .eq. 2) nbperm = nbsom-1
    do n = 1, nbperm
        kdeb = kdeb0+(n-1)*nno
        d = 0.d0
        do j = 1, nbsom
            k = iperm(kdeb+j)
            d = d+padist(3, coor1(3*(j-1)+1), coor2(3*(k-1)+1))
        end do
        if (d .lt. dmin) then
            dmin = d
            ideb = kdeb
        end if
    end do
! --- VIS A VIS DES NOEUDS DES ARRETES ET INTERIEURS
    do i = 1, nbno
        icoupl(i) = iperm(ideb+i)
    end do
8   continue
end subroutine
