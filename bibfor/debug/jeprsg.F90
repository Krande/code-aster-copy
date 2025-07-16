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
subroutine jeprsg(cunit, tgr, info)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterfort/iunifi.h"
    character(len=*) :: cunit
    real(kind=8) :: tgr
    integer(kind=8) :: info
! ----------------------------------------------------------------------
! IMPRIME UNE IMAGE OU LA LISTE DES ESPACES DISPONIBLES EN MEMOIRE
!
! IN  CUNIT : NOM LOCAL DU FICHIER D'IMPRESSION
! IN  TGR   : TAILLE DE GRAIN EN MEGA MOTS (ENTIERS)
! IN  INFO  : =1 IMPRESSION D'UNE IMAGE
!             =2 IMPRESSION DE LA LISTE
!
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
    integer(kind=8) :: istat
    common/istaje/istat(4)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
    integer(kind=8) :: idinit, idxaxd, itrech, itiad, itcol, lmots, idfr
    common/ixadje/idinit(2), idxaxd(2), itrech, itiad, itcol, lmots, idfr
! ----------------------------------------------------------------------
    aster_logical :: lamov
    character(len=132) :: chaine, init, ente, diese
    integer(kind=8) :: k, id, ida, nbc, is, nn, nc
! DEB ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iz, julist, maplac, nc0
!-----------------------------------------------------------------------
    julist = iunifi(cunit)
    if (julist .eq. 0) goto 90
    do k = 1, 132
        init(k:k) = '.'
        diese(k:k) = '#'
    end do
    ente = ' '
    do k = 1, 132
        ente(k:k) = '-'
    end do
    do k = 1, 132, 10
        ente(k:k) = '+'
    end do
    do iz = 1, 2
        id = idinit(iz)
        if (id .eq. 0) goto 80
        ida = idinit(iz)
        lamov = .false.
        if (info .eq. 1) then
            write (julist, '(A)')&
     &      'REPRESENTATION CONDENSEE DE LA SEGMENTATION MEMOIRE'
            write (julist, '(A,2(1PE11.2,A))') '1 CARACTERE = ', tgr, &
                ' MEGA-MOTS (', tgr*lois, ' MEGA-OCTETS)'
            write (julist, '(A)') ' # -> ESPACE OCCUPE '
            write (julist, '(A,/)') ' . -> ESPACE DISPONIBLE   '
        else if (info .eq. 2) then
            write (julist, '(/,10X,A)') 'LISTE DES ESPACES DISPONIBLES : '
            write (julist, '(A)')&
     &      '----------------------------------------------------------'
            write (julist, '(A)')&
     &      '|ADRESSE DE DEPART |  ADDRESSE DE FIN  |       TAILLE    |'
            write (julist, '(A)')&
     &      '----------------------------------------------------------'
        end if
        chaine = init
        nbc = 0
        nc = 1
40      continue
        nc0 = nc
        is = iszon(jiszon+id)
        if (is .eq. 0) goto 60
        nn = nint((is*lois)/(tgr*1024*1024*lois))
        nc = mod(nn, 132)+1
        if (nn/132 .gt. nbc) then
            if (info .eq. 1) then
                if (iszon(jiszon+id+3) .eq. istat(2)) then
                    chaine(nc0:) = diese
                end if
                write (julist, '(A,/,A,/,A,/)') ente, chaine, ente
                if (iszon(jiszon+id+3) .eq. istat(2)) then
                    chaine = diese
                else
                    chaine = init
                end if
                do k = nbc+1, (nn/132)-1
                    write (julist, '(A,/,A,/,A,/)') ente, chaine, ente
                end do
                chaine = init
                if (iszon(jiszon+id+3) .eq. istat(2)) then
                    chaine(1:nc) = diese
                end if
            end if
            nbc = nn/132
        end if
        if (iszon(jiszon+id+3) .eq. istat(2)) then
            chaine(nc:nc) = '#'
        else
            if (.not. lamov) then
                lamov = .true.
                ida = id
            end if
            id = is
            goto 40
        end if
60      continue
        if (lamov) then
            lamov = .false.
            maplac = id-ida-8
            if (info .eq. 2) then
                write (julist, '(A1,4X,I10,4X,A1,6X,I10,3X,A1,6X,I10,1X,A1)')&
     &        '|', ida+4, '|', id-4, '|', maplac, '|'
            end if
        end if
        id = is
        if (is .ne. 0) goto 40
        nn = nint((liszon*lois)/(tgr*1024*1024*lois))
        nc = mod(nn, 132)
        if (nn/132 .gt. nbc) then
            if (info .eq. 1) then
                do k = nbc, (nn/132)-1
                    write (julist, '(A,/,A,/,A,/)') ente, chaine, ente
                    chaine = init
                end do
            end if
        end if
        chaine(nc:132) = ' '
        chaine(nc:nc) = '|'
        if (info .eq. 1) then
            write (julist, '(A,/,A,/,A,/)') ente, chaine, ente
        else if (info .eq. 2) then
            write (julist, '(A)')&
     &      '----------------------------------------------------------'
        end if
80      continue
    end do
90  continue
! FIN ------------------------------------------------------------------
end subroutine
