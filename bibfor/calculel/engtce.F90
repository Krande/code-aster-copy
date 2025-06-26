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
subroutine engtce(ific, chamel, typtes, preci, formr)
    implicit none
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/lxlgut.h"
    integer(kind=8) :: ific
    character(len=8) :: typtes
    character(len=10) :: preci, formr
    character(len=19) :: chamel
!     COMMANDE:  ENGENDRE_TEST
!                TRAITEMENT DES SD CHAM_ELEM
!
! IN  : IFIC   : NUMERO D'UNITE IMPRESSION
! IN  : NOMSTR : NOM D'UNE SD RESULTAT
! IN  : TYPTES : TYPE DU TEST = SOMM_ABS, SOMM
! IN  : PRECI  : PRECISION POUR LE TEST_RESU
! IN  : FORMR  : FORMAT D'IMPRESSION DU CHAMP VALE REEL
! ----------------------------------------------------------------------
!
    integer(kind=8) :: vali, i, jvale, long, lg1, lg2
    real(kind=8) :: valr
    character(len=3) :: type
    character(len=80) :: form1, form2
!     ------------------------------------------------------------------
!
    call jemarq()
!
    lg1 = lxlgut(formr)
    lg2 = lxlgut(typtes)
    form1 = '(&
&            '' TYPE_TEST= '''''//typtes(1:lg2)//''''', VALE_CALC= '', '//formr(1:lg1)// &
            ', '' ), '' )'
    form2 = '( '' TYPE_TEST= '''''//typtes(1:lg2)//''''', VALE_CALC_I = '', I9, '' ), '' )'
!
    write (ific, 100)
!
    call jeveuo(chamel//'.CELV', 'L', jvale)
    call jelira(chamel//'.CELV', 'LONMAX', long)
    call jelira(chamel//'.CELV', 'TYPE', cval=type)
!
    write (ific, 101) chamel(1:8)
    write (ific, 102) preci
!
    if (type .eq. 'I') then
        if (typtes .eq. 'SOMM_ABS') then
            vali = 0
            do i = 1, long
                vali = vali+abs(zi(jvale+i-1))
            end do
        else if (typtes .eq. 'SOMM') then
            vali = 0
            do i = 1, long
                vali = vali+zi(jvale+i-1)
            end do
        else if (typtes .eq. 'MAX') then
            vali = -ismaem()
            do i = 1, long
                vali = max(vali, zi(jvale+i-1))
            end do
        else if (typtes .eq. 'MIN') then
            vali = ismaem()
            do i = 1, long
                vali = min(vali, zi(jvale+i-1))
            end do
        end if
        if (vali .eq. 0) write (ific, 112)
        write (ific, form2) vali
!
    else if (type .eq. 'R') then
        if (typtes .eq. 'SOMM_ABS') then
            valr = 0.d0
            do i = 1, long
                valr = valr+abs(zr(jvale+i-1))
            end do
        else if (typtes .eq. 'SOMM') then
            valr = 0.d0
            do i = 1, long
                valr = valr+zr(jvale+i-1)
            end do
        else if (typtes .eq. 'MAX') then
            valr = -r8maem()
            do i = 1, long
                valr = max(valr, zr(jvale+i-1))
            end do
        else if (typtes .eq. 'MIN') then
            valr = r8maem()
            do i = 1, long
                valr = min(valr, zr(jvale+i-1))
            end do
        end if
        if (abs(valr) .le. r8prem()) write (ific, 112)
        write (ific, form1) valr
    end if
!
    write (ific, 103)
!
    call jedema()
!
100 format('TEST_RESU(CHAM_ELEM= ')
!
101 format('          _F( CHAM_GD= ', a8, ', ')
!
112 format('              CRITERE= ''ABSOLU'', ')
102 format('              TOLE_MACHINE= ', a10, ',')
!
103 format('          )')
!
end subroutine
