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
subroutine impre2(licoef, liddl, linoeu, libeta, indsur, &
                  ipntrl, nbterm, typcoe, typval, irela)
    implicit none
#include "jeveux.h"
#include "asterfort/iunifi.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: indsur, ipntrl, nbterm, irela
    character(len=24) :: licoef, liddl, linoeu, libeta
!
    integer(kind=8) :: idecal, idcoef, idnoeu, iddl, idbeta, ifm
    integer(kind=8) :: idcoe, idnoe, idl, i
    real(kind=8) :: dble, dimag
    character(len=4) :: typval, typcoe
!
    ifm = iunifi('MESSAGE')
!
! --- LISTE DES COEFFICIENTS :
!     ----------------------
    call jeveuo(licoef, 'L', idcoe)
!
! --- LISTE DES DDLS :
!     --------------
    call jeveuo(liddl, 'L', idl)
!
! --- LISTE DES NOMS DES NOEUDS :
!     -------------------------
    call jeveuo(linoeu, 'L', idnoe)
!
! --- NATURE ET VALEUR DU SECOND MEMBRE DE LA RELATION LINEAIRE :
!     ---------------------------------------------------------
    call jeveuo(libeta, 'L', idbeta)
!
    idecal = ipntrl-nbterm
    idcoef = idcoe+idecal
    idnoeu = idnoe+idecal
    iddl = idl+idecal
!
    if (indsur .eq. 1) then
        write (ifm, *) 'RELATION LINEAIRE REDONDANTE ET DONC SUPPRIMEE: '
!
! ---   IMPRESSION DE LA RELATION DANS LE CAS OU LES COEFFICIENTS
! ---   SONT REELS :
!       ----------
        if (typcoe .eq. 'REEL') then
            write (ifm, 10)
            do i = 1, nbterm-1
                write (ifm, 20) zr(idcoef+i-1), zk8(iddl+i-1), zk8(idnoeu+ &
                                                                   i-1)
            end do
            write (ifm, 90) zr(idcoef+nbterm-1), zk8(iddl+nbterm-1), &
                zk8(idnoeu+nbterm-1)
!
! ---   IMPRESSION DE LA RELATION DANS LE CAS OU LES COEFFICIENTS
! ---   SONT COMPLEXES :
!       --------------
        else if (typcoe .eq. 'COMP') then
            write (ifm, 30)
            do i = 1, nbterm-1
                write (ifm, 40) dble(zc(idcoef+i-1)), dimag(zc(idcoef+i- &
                                                               1)), zk8(iddl+i-1), zk8(idnoeu+i-1)
            end do
            write (ifm, 95) dble(zc(idcoef+nbterm-1)), dimag(zc(idcoef+ &
                                                nbterm-1)), zk8(iddl+nbterm-1), zk8(idnoeu+nbterm-1)
        end if
!
        if (typval .eq. 'REEL') then
            write (ifm, 50) zr(idbeta+irela-1)
        else if (typval .eq. 'FONC') then
            write (ifm, 60) zk24(idbeta+irela-1) (1:19)
        else if (typval .eq. 'COMP') then
            write (ifm, 70) dble(zc(idbeta+irela-1)), dimag(zc(idbeta+ &
                                                               irela-1))
        end if
    end if
!
    write (ifm, 80)
10  format(2x, '    COEF     ', '*', '   DDL  ', '(', ' NOEUD  ', ')')
30  format(13x, '    COEF     ', '*', '   DDL  ', '(', ' NOEUD  ', ')')
20  format(2x, 1pe12.5, ' * ', 2x, a6, '(', a8, ')', '+')
90  format(2x, 1pe12.5, ' * ', 2x, a6, '(', a8, ')')
40  format(2x, '(', 1pe12.5, ',', 1pe12.5, ')', ' * ',&
      &       2x, a6, '(', a8, ')', '+')
95  format(2x, '(', 1pe12.5, ',', 1pe12.5, ')', ' * ',&
      &       2x, a6, '(', a8, ')')
50  format(2x, '=', 1pe12.5)
60  format(2x, '=', a19)
70  format(2x, '=', '(', 1pe12.5, ',', 1pe12.5, ')')
80  format(2x, '______________________________________', //)
end subroutine
