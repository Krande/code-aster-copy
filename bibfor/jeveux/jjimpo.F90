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
subroutine jjimpo(unit, iadmi, ideci, idatoc, genri, &
                  typei, lt, lonoi, mess)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=
    implicit none
#include "jeveux_private.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: unit
    integer(kind=8) :: iadmi, ideci, idatoc, lt, lonoi
    character(len=*) :: mess, genri, typei
! ----------------------------------------------------------------------
! ROUTINE UTILISATEUR : IMPRIME UN SEGMENT DE VALEURS
!
! IN  UNIT   : UNITE LOGIQUE D'IMPRESSION
! IN  IADMI  : ADRESSE DU PREMIER MOT DU SEGMENT DE VALEUR
! IN  IDECI  : DECALAGE PAR RAPPORT A IADMI (EN OCTETS)
! IN  IDATOC : IDENTIFICATEUR DE L'OBJET
! IN  GENRI  : GENRE DE L'OBJET
! IN  TYPEI  : TYPE DE L'OBJET
! IN  LT     : LONGUEUR DU TYPE
! IN  LONOI  : LONGEUR EN ENTIER DU SEGMENT
! IN  MESS   : MESSAGE D'INFORMATION
! ----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: icls, idco, idenom, idos, iep, ies, ips
    integer(kind=8) :: j, jdocu, jgenr, ji, jorig, jrnom, jtype
    integer(kind=8) :: k, kadm, l, ladm, n, nb, nd
    integer(kind=8) :: nl, nm, nu
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
! ----------------------------------------------------------------------
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
!
! ----------------------------------------------------------------------
    integer(kind=8) :: ideno, ilmax, iluti
    parameter(ideno=2, ilmax=4, iluti=5)
! ----------------------------------------------------------------------
    character(len=18) :: fmt
! DEB ------------------------------------------------------------------
!
    if (unit .eq. 0) goto 999
    fmt = ' '
    kadm = iadmi
    ladm = iszon(jiszon+kadm-3)
    ips = kadm-4
    iep = kadm-1
    ies = iszon(jiszon+ips)-4
    idos = iszon(jiszon+iep-1)
    icls = iszon(jiszon+ies+2)
    idco = iszon(jiszon+ies+1)
!
    if (idos .eq. 0) then
        call utmess('F', 'JEVEUX1_49')
    end if
    if (icls .lt. 0 .or. icls .gt. len(classe)) then
        call utmess('F', 'JEVEUX1_50', si=icls)
    else
        if (classe(icls:icls) .eq. ' ' .or. classe(icls:icls) .eq. '$') then
            call utmess('F', 'JEVEUX1_51', sk=classe(icls:icls))
        end if
    end if
    if (idatoc .eq. 0) then
        write (unit, "(/,' IMPRESSION SEGMENT DE VALEURS >',A,            '<')") &
            rnom(jrnom(icls)+idos) (1:32)
    else if (idatoc .eq. -1) then
        write (unit, "(/,' IMPRESSION COLLECTION ENTIERE >'            ,A,'<')") &
            rnom(jrnom(icls)+idos) (1:24)
    else if (idco .gt. 0) then
        write (unit, "(/,' IMPRESSION OBJET DE COLLECTION >'         ,A,'<  OC : ',I6)") &
            rnom(jrnom(icls)+idco) (1:24), idos
    else if (idco .eq. 0) then
        write (unit, "(/,' IMPRESSION OBJET DE COLLECTION CONTIGUE>'       ,A,'<  OC : ',I6)") &
            rnom(jrnom(icls)+idos) (1:24), idatoc
    end if
    write (unit, '(A,A)') ' >>>>> ', mess(1:min(50, len(mess)))
    if (genri .ne. 'N') then
        if (typei .eq. 'S') then
            ji = 1+((jiszon+kadm-1)*lois+ideci+ladm)/(lor8/2)
            nl = lonoi/(5*lor8/2)
            nd = mod(lonoi, (5*lor8/2))/(lor8/2)
            write (unit, '((I7,'' - '',5(I12,1X)))') (5*(l-1)+1, ( &
                                                      i4zon(ji+5*(l-1)+k-1), k=1, 5), l=1, nl)
            if (nd .ne. 0) then
                write (unit, '(I7,'' - '',5(I12,1X))') 5*nl+1, ( &
                    i4zon(ji+5*nl+k-1), k=1, nd)
            end if
        else if (typei .eq. 'I') then
            ji = jiszon+kadm+ideci/lois
            nl = lonoi/(5*lois)
            nd = mod(lonoi, (5*lois))/lois
            write (unit, '((I7,'' - '',5(I20,1X)))') (5*(l-1)+1, ( &
                                                      iszon(ji+5*(l-1)+k-1), k=1, 5), l=1, nl)
            if (nd .ne. 0) then
                write (unit, '(I7,'' - '',5(I20,1X))') 5*nl+1, ( &
                    iszon(ji+5*nl+k-1), k=1, nd)
            end if
        else if (typei .eq. 'R') then
            ji = 1+((jiszon+kadm-1)*lois+ideci+ladm)/lor8
            nl = lonoi/(5*lor8)
            nd = mod(lonoi, (5*lor8))/lor8
            write (unit, "((I7,' - ',5(1PD12.5,1X)))") (5*(l-1)+1, &
                                                        (r8zon(ji+5*(l-1)+k-1), k=1, 5), l=1, nl)
            if (nd .ne. 0) then
                write (unit, "(I7,' - ',5(1PD12.5,1X))") 5*nl+1, ( &
                    r8zon(ji+5*nl+k-1), k=1, nd)
            end if
        else if (typei .eq. 'C') then
            ji = 1+((jiszon+kadm-1)*lois+ideci+ladm)/lor8
            nl = lonoi/(2*loc8)
            nd = mod(lonoi, (2*loc8))/loc8
            write (unit, "((I7,' - ',1P,2(A1,D12.5,',',D12.5,A1)))")&
     &      (2*(l-1)+1, ('(', r8zon(ji+4*(l-1)+2*k),&
     &                      r8zon(ji+4*(l-1)+2*k+1), ')', k=0, 1), l=1, nl)
            if (nd .ne. 0) then
                write (unit, "((I7,' - ',1P,2(A1,D12.5,',',D12.5,A1)))")&
     &         2*nl+1, '(', r8zon(ji+4*(l-1)),&
     &                     r8zon(ji+4*(l-1)+1), ')'
            end if
        else if (typei .eq. 'L') then
            ji = 1+(jiszon+kadm-1)*lois+ideci+ladm
            nl = lonoi/(20*lols)
            nd = mod(lonoi, (20*lols))/lols
            write (unit, '((I7,'' - '',20(L1,1X)))') (20*(l-1)+1, ( &
                                                      lszon(ji+20*(l-1)+k-1), k=1, 20), l=1, nl)
            if (nd .ne. 0) then
                write (unit, '(I7,'' - '',20(L1,1X))') 20*nl+1, ( &
                    lszon(ji+20*nl+k-1), k=1, nd)
            end if
        else if (typei .eq. 'K') then
            ji = 1+(jiszon+kadm-1)*lois+ideci+ladm
            nb = max(65, min(81, lt+1))/(lt+1)
            nl = lonoi/(nb*lt)
            nd = (mod(lonoi, nb*lt))/lt
            write (fmt, '(I2,''(A1,'',I2,''A1,A1)'')') nb, lt
            do l = 1, nl
                write (unit, '(I7,'' - '','//fmt//')') nb*(l-1)+1, &
                    ('>', (k1zon(ji+lt*((k-1)+(l-1)*nb)+j-1), j=1, lt), '<', &
                     k=1, nb)
            end do
            if (nd .ne. 0) then
                write (unit, '(I7,'' - '','//fmt//')') nb*nl+1, &
                    ('>', (k1zon(ji+lt*((k-1)+nl*nb)+j-1), j=1, lt), '<', k=1, &
                     nd)
            end if
        else
            call utmess('F', 'JEVEUX1_43', sk=typei)
        end if
    else
        nm = iszon(jiszon+kadm-1+ilmax)
        nu = iszon(jiszon+kadm-1+iluti)
        if (nu .ne. 0) nm = nu
        idenom = iszon(jiszon+kadm-1+ideno)
        ji = 1+idenom+(jiszon+kadm-1)*lois
        nb = max(65, min(81, lt+1))/(lt+1)
        nl = nm/nb
        nd = mod(nm, nb)
        write (fmt, '(I2,''(A1,'',I2,''A1,A1)'')') nb, lt
        do l = 1, nl
            write (unit, '(I7,'' - '','//fmt//')') nb*(l-1)+1, &
                ('>', (k1zon(ji+lt*((k-1)+(l-1)*nb)+j-1), j=1, lt), '<', &
                 k=1, nb)
        end do
        if (nd .ne. 0) then
            write (unit, '(I7,'' - '','//fmt//')') nb*nl+1, ('>', ( &
                                             k1zon(ji+lt*((k-1)+nl*nb)+j-1), j=1, lt), '<', k=1, nd)
        end if
    end if
999 continue
! FIN ------------------------------------------------------------------
end subroutine
