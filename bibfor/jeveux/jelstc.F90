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
subroutine jelstc(clas, souch, ipos, maxval, klst, &
                  nbval)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/utmess.h"
    character(len=*), intent(in) :: clas, souch
    character(len=*), intent(out) :: klst(*)
    integer(kind=8), intent(in) :: ipos, maxval
    integer(kind=8), intent(out) :: nbval
! ----------------------------------------------------------------------
!  BUT : RETROUVER LES NOMS DES OBJETS DONT LE NOM CONTIENT UNE CHAINE
!        DE CARATERES DONNEE, PRESENTS SUR UNE BASE JEVEUX.
!
!  IN  : CLAS : NOM DE LA BASE : 'G', 'V', ..( ' ' -> TOUTES LES BASES )
!  IN  : SOUCH: CHAINE DE CARACTERES A CHERCHER
!  IN  : IPOS : POSITION DU DEBUT DE LA CHAINE
!               SI IPOS=0 ON REND TOUS LES NOMS
!  IN  : MAXVAL: DIMENSION DU TABLEAU KLST
!  OUT : KLST  : TABLEAU DE K24 CONTENANT LES NOMS TROUVES
!  OUT : NBVAL : NOMBRE DE NOMS TROUVES (NBVAL = -NBVAL SI < MAXVAL)
!
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: jdocu, jgenr, jorig, jrnom, jtype, l, n
    integer(kind=8) :: nbl
!-----------------------------------------------------------------------
    parameter(n=5)
!
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    aster_logical :: trouve
    character(len=6) :: pgma
    common/kappje/pgma
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    integer(kind=8) :: nrhcod, nremax, nreuti
    common/icodje/nrhcod(n), nremax(n), nreuti(n)
!     ==================================================================
    integer(kind=8) :: ncla1, ncla2, ic, j
    character(len=32) :: crnom, k32val
    character(len=1) :: kclas
!     ==================================================================
    pgma = 'JELSTC'
    l = len(souch)
    if (ipos+l .gt. 25 .or. ipos .lt. 0 .or. l .eq. 0) then
        k32val = souch
        call utmess('F', 'JEVEUX1_11', sk=k32val)
    end if
    kclas = clas(1:min(1, len(clas)))
    if (kclas .eq. ' ') then
        ncla1 = 1
        ncla2 = index(classe, '$')-1
        if (ncla2 .lt. 0) ncla2 = n
    else
        ncla1 = index(classe, kclas)
        ncla2 = ncla1
    end if
    nbl = 0
    do ic = ncla1, ncla2
        do j = 1, nremax(ic)
            crnom = rnom(jrnom(ic)+j)
            if (crnom(1:1) .eq. '?' .or. crnom(25:32) .ne. '        ') goto 150
            if (ipos .eq. 0) then
                trouve = .true.
            else if (souch .eq. crnom(ipos:ipos+l-1)) then
                trouve = .true.
            else
                trouve = .false.
            end if
            if (trouve) then
                nbl = nbl+1
                if (nbl .le. maxval) then
                    klst(nbl) = crnom(1:24)
                end if
            end if
150         continue
        end do
    end do
    if (nbl .gt. maxval) then
        nbval = -nbl
    else
        nbval = nbl
    end if
!
end subroutine
