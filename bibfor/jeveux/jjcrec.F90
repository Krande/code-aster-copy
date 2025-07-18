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

subroutine jjcrec(icl, ida, genri, typei, nb, &
                  iadmi)
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
#include "jeveux_private.h"
#include "asterfort/jjallt.h"
#include "asterfort/jjecrs.h"
#include "asterfort/jjprem.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: icl, ida, nb, iadmi
    character(len=*) :: genri, typei
! ----------------------------------------------------------------------
! CREATION D'UN OBJET SIMPLE ATTRIBUT COMPOSANT UNE COLLECTION
!
! IN  ICL   : CLASSE ASSOCIEE
! IN  IDA   : IDENTIFICATEUR
! IN  GENRI : GENRE DE L'OS
! IN  TYPEI : TYPE DE L'OS
! IN  NB    : LONGEUR EN MOTS DU SEGMENT DE VALEUR ASSOCIE
! OUT IADMI : ADRESSE DU SEGMENT DE VALEUR DANS ISZON
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iadyn, iv, jcara, jdate, jdocu, jgenr, jhcod
    integer(kind=8) :: jiadd, jiadm, jlong, jlono, jltyp, jluti
    integer(kind=8) :: jmarq, jorig, jrnom, jtype, l, lonoi, n
    integer(kind=8) :: nbl, nhc
!-----------------------------------------------------------------------
    parameter(n=5)
    common/jiatje/jltyp(n), jlong(n), jdate(n), jiadd(n), jiadm(n),&
     &                 jlono(n), jhcod(n), jcara(n), jluti(n), jmarq(n)
    common/jkatje/jgenr(n), jtype(n), jdocu(n), jorig(n), jrnom(n)
!     ------------------------------------------------------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
!     ------------------------------------------------------------------
    character(len=4) :: ifmt
    integer(kind=8) :: ilorep, ideno, ilnom, ilmax, iluti, idehc
    parameter(ilorep=1, ideno=2, ilnom=3, ilmax=4, iluti=5, idehc=6)
    integer(kind=8) :: irt
! DEB ------------------------------------------------------------------
    irt = 0
    genr(jgenr(icl)+ida) = genri(1:1)
    type(jtype(icl)+ida) = typei(1:1)
    if (genri .eq. 'N' .and. typei(1:1) .ne. 'K') then
        call utmess('F', 'JEVEUX1_38')
    end if
    if (typei(1:1) .eq. 'K') then
        l = len(typei)
        if (l .eq. 1) then
            call utmess('F', 'JEVEUX1_39')
        end if
        write (ifmt, '(''(I'',I1,'')'')') l-1
        read (typei(2:l), ifmt) iv
        if (iv .le. 0 .or. iv .gt. 512) then
            call utmess('F', 'JEVEUX1_40', si=iv)
        end if
        if (genri .eq. 'N') then
            if (mod(iv, lois) .ne. 0) then
                call utmess('F', 'JEVEUX1_41', si=iv)
            end if
            if (iv .gt. 24) then
                call utmess('F', 'JEVEUX1_42', si=iv)
            end if
        end if
    else if (typei(1:1) .eq. 'S') then
        iv = lor8/2
    else if (typei(1:1) .eq. 'I') then
        iv = lois
    else if (typei(1:1) .eq. 'R') then
        iv = lor8
    else if (typei(1:1) .eq. 'C') then
        iv = loc8
    else if (typei(1:1) .eq. 'L') then
        iv = lols
    else
        call utmess('F', 'JEVEUX1_43', sk=typei(1:1))
    end if
    ltyp(jltyp(icl)+ida) = iv
    iadm(jiadm(icl)+2*ida-1) = 0
    iadm(jiadm(icl)+2*ida) = 0
    if (nb .gt. 0) then
        if (genri .eq. 'N') then
            long(jlong(icl)+ida) = nb
            lonoi = (idehc+jjprem(nb, irt))*lois+(nb+1)*iv
            if (mod(lonoi, iv) .gt. 0) then
                lono(jlono(icl)+ida) = (lonoi/iv)+1
            else
                lono(jlono(icl)+ida) = (lonoi/iv)
            end if
        else if (typei(1:1) .eq. 'C') then
            long(jlong(icl)+ida) = 2*nb
            lono(jlono(icl)+ida) = 2*nb
        else
            long(jlong(icl)+ida) = nb
            lono(jlono(icl)+ida) = nb
        end if
        nbl = lono(jlono(icl)+ida)*iv
        call jjallt(nbl, icl, genri, typei, iv, &
                    'INIT', iadmi, iadyn)
        iadm(jiadm(icl)+2*ida-1) = iadmi
        iadm(jiadm(icl)+2*ida) = iadyn
        call jjecrs(iadmi, icl, ida, 0, 'E', &
                    imarq(jmarq(icl)+2*ida-1))
        if (genri .eq. 'N') then
            nhc = jjprem(nb, irt)
            iszon(jiszon+iadmi-1+ilorep) = nhc
            iszon(jiszon+iadmi-1+ideno) = (idehc+nhc)*lois
            iszon(jiszon+iadmi-1+ilnom) = iv
            iszon(jiszon+iadmi-1+ilmax) = nb
            iszon(jiszon+iadmi-1+iluti) = 0
            iszon(jiszon+iadmi-1+idehc) = idehc
        end if
    else if (genri .eq. 'E') then
        long(jlong(icl)+ida) = 1
    end if
! FIN ------------------------------------------------------------------
end subroutine
