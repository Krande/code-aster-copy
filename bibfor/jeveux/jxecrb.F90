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
subroutine jxecrb(ic, iaddi, iadmo, lso, idco, &
                  idos)
! person_in_charge: j-pierre.lefebvre at edf.fr
! aslint: disable=W1303
! for the path name
    implicit none
#include "asterf_types.h"
#include "jeveux_private.h"
#include "asterc/writdr.h"
#include "asterfort/get_jvbasename.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: ic, iaddi, iadmo, lso, idco, idos
! ----------------------------------------------------------------------
! ECRITURE DISQUE D'UN OU PLUSIEURS ENREGISTREMENTS SUR LE FICHIER
! D'ACCES DIRECT ASSOCIE A UNE BASE
!
! IN  IC    : NOM DE LA CLASSE
! IN  IADDI : ADRESSE DISQUE DU SEGMENT DE VALEURS
! IN  IADMO : ADRESSE MEMOIRE DU SEGMENT DE VALEURS EN OCTET
! IN  LSO   : LONGUEUR EN OCTET DU SEGMENT DE VALEURS
! IN  IDCO  : IDENTIFICATEUR DE COLLECTION
! IN  IDOS  : IDENTIFICATEUR D'OBJET SIMPLE OU D'OBJET DE COLLECTION
!                                              >0 GROS OBJET
!                                              =0 PETITS OBJETS
! ----------------------------------------------------------------------
    integer(kind=8) :: lk1zon, jk1zon, liszon, jiszon
    common/izonje/lk1zon, jk1zon, liszon, jiszon
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadloc, ib, ierr, jiacce
    integer(kind=8) :: jiecr, jusadi, n, nbacce, nblent, numext
!-----------------------------------------------------------------------
    parameter(n=5)
!     ------------------------------------------------------------------
    integer(kind=8) :: lbis, lois, lols, lor8, loc8
    common/ienvje/lbis, lois, lols, lor8, loc8
!     ------------------------------------------------------------------
    integer(kind=8) :: nblmax, nbluti, longbl, kitlec, kitecr, kiadm, iitlec, iitecr
    integer(kind=8) :: nitecr, kmarq
    common/ificje/nblmax(n), nbluti(n), longbl(n),&
     &                 kitlec(n), kitecr(n), kiadm(n),&
     &                 iitlec(n), iitecr(n), nitecr(n), kmarq(n)
!
    character(len=2) :: dn2
    character(len=5) :: classe
    character(len=8) :: nomfic, kstout, kstini
    common/kficje/classe, nomfic(n), kstout(n), kstini(n),&
     &                 dn2(n)
    character(len=8) :: nombas
    common/kbasje/nombas(n)
    character(len=128) :: repglo, repvol
    common/banvje/repglo, repvol
    integer(kind=8) :: lrepgl, lrepvo
    common/balvje/lrepgl, lrepvo
    integer(kind=8) :: idn, iext, nbenrg
    common/iextje/idn(n), iext(n), nbenrg(n)
    common/jiacce/jiacce(n), nbacce(2*n)
    common/jusadi/jusadi(n)
!     ------------------------------------------------------------------
    character(len=512) :: nom512
    aster_logical :: lrab
    integer(kind=8) :: lgbl, vali(3)
! DEB ------------------------------------------------------------------
    ib = 0
    ierr = 0
    lgbl = 1024*longbl(ic)*lois
    nblent = lso/lgbl
    lrab = (mod(lso, lgbl) .ne. 0)
!     ------------------------------------------------------------------
    do i = 1, nblent
        numext = (iaddi+i-2)/nbenrg(ic)
        iadloc = (iaddi+i-1)-(numext*nbenrg(ic))
        call get_jvbasename(nomfic(ic) (1:4), numext+1, nom512)
        jiecr = (jk1zon+iadmo-1+lgbl*(i-1))/lois+1
        call writdr(nom512, iszon(jiecr), lgbl, iadloc, ierr)
        if (ierr .ne. 0) then
            vali(1) = iaddi+i-1
            vali(2) = numext
            vali(3) = ierr
            call utmess('F', 'JEVEUX_40', sk=nombas(ic), ni=3, vali=vali)
        end if
        nbacce(2*ic) = nbacce(2*ic)+1
        iusadi(jusadi(ic)+3*(iaddi+i-1)-2) = idco
        iusadi(jusadi(ic)+3*(iaddi+i-1)-1) = idos
    end do
    iacce(jiacce(ic)+iaddi) = iacce(jiacce(ic)+iaddi)+1
    if (lrab) then
        numext = (iaddi+nblent-1)/nbenrg(ic)
        iadloc = (iaddi+nblent)-(numext*nbenrg(ic))
        call get_jvbasename(nomfic(ic) (1:4), numext+1, nom512)
        jiecr = (jk1zon+iadmo-1+lso-lgbl)/lois+1
        if (lso .lt. lgbl) then
            jiecr = (jk1zon+iadmo-1)/lois+1
        end if
        call writdr(nom512, iszon(jiecr), lgbl, iadloc, ierr)
        if (ierr .ne. 0) then
            vali(1) = iaddi+i-1
            vali(2) = numext
            vali(3) = ierr
            call utmess('F', 'JEVEUX_40', sk=nombas(ic), ni=3, vali=vali)
        end if
        nbacce(2*ic) = nbacce(2*ic)+1
        iusadi(jusadi(ic)+3*(iaddi+nblent)-2) = idco
        iusadi(jusadi(ic)+3*(iaddi+nblent)-1) = idos
    end if
! FIN ------------------------------------------------------------------
end subroutine
