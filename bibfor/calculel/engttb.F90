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
subroutine engttb(ific, nomsd, typtes, preci, formr)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
#include "asterfort/tbexip.h"
!
    integer(kind=8) :: ific
    character(len=8) :: typtes
    character(len=10) :: preci, formr
    character(len=19) :: nomsd
!     COMMANDE:  ENGENDRE_TEST
!                TRAITEMENT DES SD TABLE
!
! IN  : IFIC   : NUMERO D'UNITE IMPRESSION
! IN  : NOMSD : NOM D'UNE SD RESULTAT
! IN  : TYPTES : TYPE DU TEST = SOMM_ABS, SOMM
! IN  : PRECI  : PRECISION POUR LE TEST_TABLE
! IN  : FORMR  : FORMAT D'IMPRESSION DU CHAMP VALE REEL
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbpara, nblign, vali, ipar, lg, lg1, lg2, i, jvale, jvall
    real(kind=8) :: valr
    aster_logical :: exist
    character(len=3) :: type
    character(len=16) :: nomsym
    character(len=90) :: form1, form2, form3
    integer(kind=8), pointer :: tbnp(:) => null()
    character(len=24), pointer :: tblp(:) => null()
!     ------------------------------------------------------------------
!
    call jemarq()
!
    lg1 = lxlgut(formr)
    lg2 = lxlgut(typtes)
    form1 = '(&
&            '' TYPE_TEST= '''''//typtes(1:lg2)//''''', VALE_CALC= '', '//formr(1:lg1)// &
            ', '' )''  )'
    form2 = '( '' TYPE_TEST= '''''//typtes(1:lg2)//''''', VALE_CALC_I = '', I9, '' )'' )'
!
    call jeveuo(nomsd//'.TBLP', 'L', vk24=tblp)
    call jeveuo(nomsd//'.TBNP', 'L', vi=tbnp)
    nbpara = tbnp(1)
    nblign = tbnp(2)
!
    do ipar = 1, nbpara
!
        nomsym = tblp(1+4*(ipar-1))
        call tbexip(nomsd, nomsym, exist, type)
        if (.not. exist) goto 400
!
        lg = lxlgut(nomsym)
        call jeveuo(tblp(1+4*(ipar-1)+2), 'L', jvale)
        call jeveuo(tblp(1+4*(ipar-1)+3), 'L', jvall)
!
        form3 = '( ''TEST_TABLE(TABLE= '',A8,'', NOM_PARA= '''''//nomsym(1:lg)//''''', '' )'
!
        if (type .eq. 'I') then
!             -------------
            write (ific, form3) nomsd(1:8)
            write (ific, 402) preci
!
            if (typtes .eq. 'SOMM_ABS') then
                vali = 0
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) vali = vali+abs(zi(jvale+i-1))
                end do
            else if (typtes .eq. 'SOMM') then
                vali = 0
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) vali = vali+zi(jvale+i-1)
                end do
            else if (typtes .eq. 'MAX') then
                vali = -ismaem()
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) vali = max(vali, zi(jvale+i-1))
                end do
            else if (typtes .eq. 'MIN') then
                vali = ismaem()
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) vali = min(vali, zi(jvale+i-1))
                end do
            end if
            if (vali .eq. 0) write (ific, 401)
            write (ific, form2) vali
!
        else if (type .eq. 'R') then
!                 -------------
            write (ific, form3) nomsd(1:8)
            write (ific, 402) preci
!
            if (typtes .eq. 'SOMM_ABS') then
                valr = 0.d0
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) valr = valr+abs(zr(jvale+i-1))
                end do
            else if (typtes .eq. 'SOMM') then
                valr = 0.d0
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) valr = valr+zr(jvale+i-1)
                end do
            else if (typtes .eq. 'MAX') then
                valr = -r8maem()
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) valr = max(valr, zr(jvale+i-1))
                end do
            else if (typtes .eq. 'MIN') then
                valr = r8maem()
                do i = 1, nblign
                    if (zi(jvall+i-1) .eq. 1) valr = min(valr, zr(jvale+i-1))
                end do
            end if
            if (abs(valr) .le. r8prem()) write (ific, 401)
            write (ific, form1) valr
        end if
400     continue
    end do
!
    call jedema()
!
401 format('            CRITERE= ''ABSOLU'', ')
!
402 format('            TOLE_MACHINE= ', a10, ',')
!
end subroutine
