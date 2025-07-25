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

subroutine mtcmbi(typmat, lmat, coef, ccoef, lres)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtconl.h"
#include "asterfort/mtdsc2.h"
#include "asterfort/pteddl.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: lmat, lres
    character(len=*) :: typmat
    complex(kind=8) :: ccoef
!     DUPLIQUE LA MATRICE EN METTANT TOUTES LES TERMES A ZERO SAUF
!     LES "LAGRANGE" EN LEUR APPLIQUANT UN COEFFICIENT.
!     -----------------------------------------------------------------
! IN  K* TYPMAT = TYPE DE MATRICE   (R OU C)
! IN  I  LMAT   = POINTEUR DE MATRICE
! IN  I  LRES   = POINTEUR DE MATRICE RESULTAT
!     -----------------------------------------------------------------
!     NBBLIC = NOMBRE DE BLOCS POUR .VALI DE LA MATRICE
!     LGBLOC = LONGUEUR DES BLOCS
!     -----------------------------------------------------------------
    integer(kind=8) :: lgbloc
    real(kind=8) :: const(2)
    character(len=1) :: ch1, typcst
    character(len=8) :: nomddl
    character(len=14) :: nume
    character(len=19) :: matres, noma
    character(len=24) :: valm, valmr
    complex(kind=8) :: czero
    aster_logical :: matsym
!     -----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iatmai, iatmat, iatrei, iatres, ibid, icoef
    integer(kind=8) :: idebli, iequa, ifinli, ilig, ind, ival
    integer(kind=8) :: jsmdi, jsmhc, kin, lddl, neq
!
    real(kind=8) :: coef, zero
    character(len=24), pointer :: refa(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
!
    if (typmat(1:1) .ne. 'R' .and. typmat(1:1) .ne. 'C') then
        ch1 = typmat(1:1)
        call utmess('F', 'ALGELINE2_6', sk=ch1)
    end if
!
!     --- AFFE_CHAR_CINE ? ---
!
    if (zi(lmat+7) .ne. 0) then
        call utmess('F', 'ALGELINE2_7')
    end if
!
    zero = 0.d0
    czero = dcmplx(zero, zero)
    matsym = .true.
!
    if (zi(lmat+4) .ne. 1) matsym = .false.
    noma = zk24(zi(lmat+1)) (1:19)
    valm = noma//'.VALM'
!
    neq = zi(lres+2)
    call mtdsc2(zk24(zi(lres+1)), 'SMDI', 'L', jsmdi)
    lgbloc = zi(lres+14)
    matres = zk24(zi(lres+1)) (1:19)
    call jeveuo(matres//'.REFA', 'L', vk24=refa)
    call jeveuo(refa(2) (1:14)//'.SMOS.SMHC', 'L', jsmhc)
    call jeveuo(refa(2) (1:14)//'.SMOS.SMDI', 'L', ibid)
    ASSERT(ibid .eq. jsmdi)
!
    valmr = matres//'.VALM'
!
!     --- NOM DE LA NUMEROTATION ASSOCIEE A LA MATRICE ---
    call dismoi('NOM_NUME_DDL', noma, 'MATR_ASSE', repk=nume)
!
!
!     --- TOUTES COMPOSANTES A ZERO SAUF LES LAGRANGES ---
    nomddl = 'LAGR    '
    call wkvect('&&MTCMBI', 'V V I', neq, lddl)
    call pteddl('NUME_DDL', nume, 1, nomddl, neq, &
                list_equa=zi(lddl))
    do i = 0, neq-1
        zi(lddl+i) = 1-zi(lddl+i)
    end do
!
!
!
    call jeveuo(jexnum(valmr, 1), 'E', iatres)
    if (.not. matsym) then
        call jeveuo(jexnum(valmr, 2), 'E', iatrei)
    end if
!
    if (typmat(1:1) .eq. 'R') then
        do ival = iatres, iatres+lgbloc-1
            zr(ival) = zero
        end do
        if (.not. matsym) then
            do ival = iatrei, iatrei+lgbloc-1
                zr(ival) = zero
            end do
        end if
    else
        do ival = iatres, iatres+lgbloc-1
            zc(ival) = czero
        end do
    end if
!
    call jeveuo(jexnum(valm, 1), 'L', iatmat)
    if (.not. matsym) then
        call jeveuo(jexnum(valm, 2), 'E', iatmai)
    end if
!
!
    if (typmat(1:1) .eq. 'R') then
        kin = 0
        idebli = 1
        do iequa = 1, neq
            ifinli = zi(jsmdi+iequa-1)
            do ind = idebli, ifinli
                kin = kin+1
                ilig = zi4(jsmhc-1+kin)
                icoef = min((2-zi(lddl+ilig-1)-zi(lddl+iequa-1)), 1)
                zr(iatres+kin-1) = zr(iatres+kin-1)+zr(iatmat+kin-1)*icoef*coef
            end do
            idebli = zi(jsmdi+iequa-1)+1
        end do
!
!
    else if (typmat(1:1) .eq. 'C') then
        kin = 0
        idebli = 1
        do iequa = 1, neq
            ifinli = zi(jsmdi+iequa-1)
            do ind = idebli, ifinli
                kin = kin+1
                ilig = zi4(jsmhc-1+kin)
                icoef = min((2-zi(lddl+ilig-1)-zi(lddl+iequa-1)), 1)
                zc(iatres+kin-1) = zc(iatres+kin-1)+zc(iatmat+kin-1)*icoef*ccoef
            end do
            idebli = zi(jsmdi+iequa-1)+1
        end do
    end if
!
!
    call jelibe(jexnum(valm, 1))
    if (.not. matsym) then
        call jelibe(jexnum(valm, 2))
    end if
    call jelibe(jexnum(valmr, 1))
    if (.not. matsym) then
        call jelibe(jexnum(valmr, 2))
    end if
!
!
!     --- ACTUALISATION DU .CONL ----
    if (typmat(1:1) .eq. 'R') then
        typcst = 'R'
        const(1) = 1.d0
    else
        typcst = 'C'
        const(1) = 1.d0
        const(2) = 1.d0
    end if
    call mtconl(1, typcst, const, [lmat], typmat, &
                lres)
!
    call jedetr('&&MTCMBI')
!
!
    call jedema()
end subroutine
