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
subroutine lxinit()
! aslint: disable=
    implicit none
!     INITIALISATION DE L'ANALYSEUR LEXICAL
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     ROUTINE(S) UTILISEE(S) :
!         (CF ARGUMENT)
!     ROUTINE(S) FORTRAN     :
!         CHAR    ICHAR
!     ------------------------------------------------------------------
! FIN LXINIT
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
#include "asterfort/lxdeli.h"
    integer(kind=8) :: i, mxchar, mxclas, mxcols, mxdeli, nbdeli
!-----------------------------------------------------------------------
    parameter(mxclas=10, mxchar=255, mxdeli=15)
    integer(kind=8) :: clnum, cllet, clsig, clpnt, clexp, clquo, clbls, clbl, clill
    integer(kind=8) :: cleor
!
    common/lxcn01/clnum, cllet, clsig, clpnt, clexp, clquo,&
     &                  clbls, clbl, clill, cleor, nbdeli
!
    character(len=1) :: class(0:mxchar), cldeli(mxdeli)
    common/lxcc01/class, cldeli
!
!     ------------------------------------------------------------------
    parameter(mxcols=80)
    character(len=mxcols) :: chaine
    character(len=1) :: kclass
!     ------------------------------------------------------------------
!
!
!     ------------------------------------------------------------------
!                     DEFINITION DES CLASSES SIMPLES
!     CLNUM  =  1 : NUMERIQUES        CLLET  =  2 : LETTRES
!     CLSIG  =  3 : SIGNE + -         CLPNT  =  4 : POINT .
!     CLEXP  =  5 : EXPOSANT E D      CLQUO  =  6 : QUOTE '
!     CLBLS  =  7 : BLANC SOULIGNE _  CLBL   =  8 : BLANC
!     CLILL  =  9 : ILLEGAUX          CLEOR  = 10 : FIN D'ENREGISTREMENT
!     ------------------------------------------------------------------
!
    clnum = 1
    cllet = 2
    clsig = 3
    clpnt = 4
    clexp = 5
    clquo = 6
    clbls = 7
    clbl = 8
    clill = 9
    cleor = 10
!     ------------------------------------------------------------------
!
!
!     INITIALISATION DES CLASSES A ILLEGAL /* OPTION PAR DEFAUT */
    kclass = char(clill)
    do i = 0, mxchar
        class(i) = kclass
    end do
!
!     INITIALISATION DE LA CLASSE NUMERIQUE
    chaine = '0123456789'
    kclass = char(clnum)
    do i = 1, 10
        class(ichar(chaine(i:i))) = kclass
    end do
!
!     INITIALISATION DE LA CLASSE ALPHABETIQUE
    chaine = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'//'abcdefghijklmnopqrstuvwxyz'
    kclass = char(cllet)
    do i = 1, 52
        class(ichar(chaine(i:i))) = kclass
    end do
!
!     INITIALISATION DE LA CLASSE SIGNE
    class(ichar('+')) = char(clsig)
    class(ichar('-')) = char(clsig)
!
!     INITIALISATION DE LA CLASSE EXPOSANT
    class(ichar('E')) = char(clexp)
    class(ichar('e')) = char(clexp)
    class(ichar('D')) = char(clexp)
    class(ichar('d')) = char(clexp)
!
!     INITIALISATION DE LA CLASSE QUOTE BLANC BLANC_SOULIGNE ET POINT
    class(ichar('''')) = char(clquo)
    class(ichar(' ')) = char(clbl)
    class(ichar('_')) = char(clbls)
    class(ichar('.')) = char(clpnt)
!
!     TABULATION
    class(9) = char(clbl)
!
!     INITIALISATION DE LA CLASSE 'DELIMITEUR'
    nbdeli = mxdeli
    call lxdeli(cldeli, nbdeli)
    kclass = char(clill)
    do i = 1, nbdeli
        if (class(ichar(cldeli(i))) .eq. kclass) then
            class(ichar(cldeli(i))) = char(mxclas+i)
        end if
    end do
!
!     ------------------------------------------------------------------
!
end subroutine
