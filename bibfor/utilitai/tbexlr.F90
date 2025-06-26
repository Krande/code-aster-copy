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
subroutine tbexlr(nomta, listr, basout)
    implicit none
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=*) :: nomta, listr, basout
!     TRANSFORMER UNE TABLE EN LISTR8 QUE SI CETTE TABLE EST
!                 "DIAGONALISABLE PAR BLOCS"
! ----------------------------------------------------------------------
! IN  : NOMTA  : NOM DE LA SD "TABLE".
! IN  : LISTR  : NOM DE LA SD "LISTR8" RESULTAT
! IN  : BASOUT : BASE DE CREATION DE "LISTR"
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
    integer(kind=8) :: iret, nbpara, nblign, nbpr, nblg, ipar, ndim
    integer(kind=8) :: i, j, jvale, jlogq, nbvale, k1, kcol, klig, ideb1, ideb2
    integer(kind=8) :: ifin1, ifin2, nbcl, ivide, ilig, ibloc, jpas, jnbp, jbor, k
    integer(kind=8) :: kcol1, kcol2
    character(len=1) :: base
    character(len=3) :: type
    character(len=19) :: nomtab, listr8
    character(len=24) :: nomjv, nomjvl
    integer(kind=8), pointer :: colonn(:) => null()
    integer(kind=8), pointer :: lignes(:) => null()
    integer(kind=8), pointer :: nume_lign(:) => null()
    integer(kind=8), pointer :: nume_para(:) => null()
    real(kind=8), pointer :: vale_r(:) => null()
    character(len=24), pointer :: tblp(:) => null()
    integer(kind=8), pointer :: tbnp(:) => null()
! ----------------------------------------------------------------------
!
    call jemarq()
!
    nomtab = ' '
    nomtab = nomta
    call jeexin(nomtab//'.TBBA', iret)
    if (iret .eq. 0) then
        call utmess('F', 'UTILITAI4_64')
    end if
    if (nomtab(18:19) .ne. '  ') then
        call utmess('F', 'UTILITAI4_68')
    end if
    base = basout(1:1)
!
    call jeveuo(nomtab//'.TBNP', 'E', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
    if (nbpara .eq. 0) then
        call utmess('F', 'UTILITAI4_65')
    end if
    if (nblign .eq. 0) then
        call utmess('F', 'UTILITAI4_76')
    end if
!
    call jeveuo(nomtab//'.TBLP', 'L', vk24=tblp)
!
!     --- ON NE RETIENT QUE LES PARAMETRES DE TYPE "I" ET "R" ---
!
    AS_ALLOCATE(vi=nume_para, size=nbpara)
    nbpr = 0
    do i = 1, nbpara
        type = tblp(1+4*(i-1)+1)
        if (type(1:1) .eq. 'I') then
            nbpr = nbpr+1
            nume_para(nbpr) = i
        else if (type(1:1) .eq. 'R') then
            nbpr = nbpr+1
            nume_para(nbpr) = i
        end if
    end do
    if (nbpr .eq. 0) then
        call utmess('F', 'UTILITAI4_81')
    end if
!
!     --- ON NE RETIENT QUE LES LIGNES NON VIDES ---
!
    AS_ALLOCATE(vi=nume_lign, size=nblign)
    nblg = 0
    do i = 1, nblign
        nbcl = 0
        do j = 1, nbpr
            ipar = nume_para(j)
            nomjvl = tblp(1+4*(ipar-1)+3)
            call jeveuo(nomjvl, 'L', jlogq)
            if (zi(jlogq+i-1) .eq. 1) nbcl = nbcl+1
        end do
        if (nbcl .ne. 0) then
            nblg = nblg+1
            nume_lign(nblg) = i
        end if
    end do
    if (nblg .eq. 0) then
        call utmess('F', 'UTILITAI4_82')
    end if
!
!     --- RECHERCHE DE BLOCS ---
!
    nbvale = nbpr*nblg
    AS_ALLOCATE(vr=vale_r, size=nbvale)
    AS_ALLOCATE(vi=colonn, size=nbpr)
    AS_ALLOCATE(vi=lignes, size=nblg)
!
    ibloc = 1
    k1 = 0
    do i = 1, nblg
        ilig = nume_lign(i)
        ideb1 = 0
        ifin1 = nbpr
        ivide = 0
        kcol1 = 0
        do j = 1, nbpr
            ipar = nume_para(j)
            type = tblp(1+4*(ipar-1)+1)
            nomjv = tblp(1+4*(ipar-1)+2)
            nomjvl = tblp(1+4*(ipar-1)+3)
            call jeveuo(nomjv, 'L', jvale)
            call jeveuo(nomjvl, 'L', jlogq)
            if (zi(jlogq+ilig-1) .eq. 1) then
                if (ideb1 .eq. 0) ideb1 = ipar
                kcol1 = kcol1+1
                if (ivide .eq. 1) then
                    call utmess('F', 'UTILITAI4_83')
                end if
                if (type(1:1) .eq. 'I') then
                    k1 = k1+1
                    vale_r(k1) = zi(jvale+ilig-1)
                else if (type(1:1) .eq. 'R') then
                    k1 = k1+1
                    vale_r(k1) = zr(jvale+ilig-1)
                end if
            else
                if (ideb1 .ne. 0) then
                    ivide = 1
                    if (ifin1 .eq. nbpr) ifin1 = nume_para(1+j-1-1)
                end if
!               IF ( IFIN1 .EQ. 0 ) IFIN1 = ZI(KPARA+J-1-1)
            end if
        end do
        if (i .eq. 1) then
            klig = 1
        else
            if (ideb1 .eq. ideb2 .and. ifin1 .eq. ifin2) then
                klig = klig+1
            else
!              --- NOUVEAU BLOC ---
                colonn(ibloc) = kcol2
                lignes(ibloc) = klig
                ibloc = ibloc+1
                klig = 1
            end if
        end if
        ideb2 = ideb1
        ifin2 = ifin1
        kcol2 = kcol1
    end do
    colonn(ibloc) = kcol2
    lignes(ibloc) = klig
!
!     --- ON STOCKE ---
!
    nbvale = 1+(2*ibloc)
    do i = 1, ibloc
        kcol = colonn(i)
        klig = lignes(i)
        nbvale = nbvale+(kcol*klig)
    end do
!
    listr8 = listr
    ndim = max(1, nbvale-1)
    call wkvect(listr8//'.LPAS', base//' V R', ndim, jpas)
    call wkvect(listr8//'.NBPA', base//' V I', ndim, jnbp)
    call wkvect(listr8//'.BINT', base//' V R', nbvale, jbor)
    call wkvect(listr8//'.VALE', base//' V R', nbvale, jvale)
!
    zr(jvale) = ibloc
    j = 1
    k1 = 0
    do i = 1, ibloc
        kcol = colonn(i)
        j = j+1
        zr(jvale+j-1) = kcol
        klig = lignes(i)
        j = j+1
        zr(jvale+j-1) = klig
        ndim = kcol*klig
        do k = 1, ndim
            k1 = k1+1
            j = j+1
            zr(jvale+j-1) = vale_r(k1)
        end do
    end do
!
    do i = 1, nbvale-1
        zr(jpas+i-1) = zr(jvale+i)-zr(jvale+i-1)
        zi(jnbp+i-1) = 1
        zr(jbor+i-1) = zr(jvale+i-1)
    end do
    zr(jbor+nbvale-1) = zr(jvale+nbvale-1)
!
    AS_DEALLOCATE(vi=nume_para)
    AS_DEALLOCATE(vi=nume_lign)
    AS_DEALLOCATE(vr=vale_r)
    AS_DEALLOCATE(vi=colonn)
    AS_DEALLOCATE(vi=lignes)
!
    call jedema()
end subroutine
