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

subroutine diagav(noma19, neq, ilfin, typvar, eps)
    implicit none
!     BUT : AJOUTER L'OBJET .DIGS A UNE MATR_ASSE
!           ET CALCULER UN EPSILON NUMERIQUE POUR LA FACTORISATION
!     IN  : NOMA19 : MATR_ASSE QUE L'ON COMPLETERA PAR L'OBJET .DIGS
!     IN  : NEQ    : NOMBRE D'EQUATIONS
!     IN  : ILFIN  : NUMERO DE LA LIGNE DE FIN DE FACTORISITION
!     IN  : TYPVAR : REEL/COMPLEXE
!     OUT : EPS    : 'EPSILON' TEL QU'UN TERME DIAGONAL APRES
!                    FACTORISATION SERA CONSIDERE COMME NUL
!     ------------------------------------------------------------------

#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=19) :: noma19
    character(len=14) :: nu
    character(len=1) :: base
    character(len=4) :: kmpic
    real(kind=8) :: eps, diamax, diamin, vabs
    integer(kind=8) :: neq, ilfin, typvar, ifm, niv, iret, iadigs
    integer(kind=8) :: jsxdi, nbbloc, ibloc, iavale, idern, iprem, i
    integer(kind=8), pointer :: scbl(:) => null()
    integer(kind=8), pointer :: scib(:) => null()
    character(len=24), pointer :: refa(:) => null()
!     ------------------------------------------------------------------
!
!
    call jemarq()
    call infdbg('FACTOR', ifm, niv)
!
    call dismoi('MPI_COMPLET', noma19, 'MATR_ASSE', repk=kmpic)
    if (kmpic .ne. 'OUI') then
        call utmess('F', 'CALCULEL6_54')
    end if
    call jeveuo(noma19//'.REFA', 'L', vk24=refa)
    call jelira(noma19//'.REFA', 'CLAS', cval=base)
    ASSERT(refa(3) .ne. 'ELIML')
!
!
!     -- ALLOCATION ET CALCUL DE L'OBJET .DIGS :
!        CET OBJET CONTIENDRA LES TERMES DIAGONAUX
!        AVANT ET APRES FACTORISATION :
!        (1->NEQ : AVANT , NEQ+1 ->2*NEQ : APRES )
!     -----------------------------------------
    call jedetr(noma19//'.DIGS')
    if (typvar .eq. 1) then
        call wkvect(noma19//'.DIGS', base//' V R', 2*neq, iadigs)
    else
        call wkvect(noma19//'.DIGS', base//' V C', 2*neq, iadigs)
    end if
    call dismoi('NOM_NUME_DDL', noma19, 'MATR_ASSE', repk=nu)
!
!
!     CAS STOCKAGE MORSE DISPONIBLE (OBJET .VALM):
!     ---------------------------------------------
    call jeexin(noma19//'.VALM', iret)
    if (iret .gt. 0) then
        call jeveuo(nu//'.SMOS.SMDI', 'L', jsxdi)
        call jeveuo(jexnum(noma19//'.VALM', 1), 'L', iavale)
        if (typvar .eq. 1) then
            do i = 1, neq
                zr(iadigs-1+i) = zr(iavale-1+zi(jsxdi+i-1))
            end do
        else if (typvar .eq. 2) then
            do i = 1, neq
                zc(iadigs-1+i) = zc(iavale-1+zi(jsxdi+i-1))
            end do
        else
            ASSERT(.false.)
        end if
        goto 9998
    end if
!
!
!     CAS STOCKAGE MORSE INDISPONIBLE (OBJET .VALM):
!     ---------------------------------------------
    ASSERT((noma19 .eq. '&&OP0070.RESOC.MATC') .or. (noma19 .eq. '&&OP0070.RESUC.MATC'))
    call jeveuo(nu//'.SLCS.SCDI', 'L', jsxdi)
    call jeveuo(nu//'.SLCS.SCBL', 'L', vi=scbl)
    call jeveuo(nu//'.SLCS.SCIB', 'L', vi=scib)
    nbbloc = scib(ilfin)
    do ibloc = 1, nbbloc
        call jeveuo(jexnum(noma19//'.UALF', ibloc), 'L', iavale)
        idern = scbl(ibloc+1)
        ASSERT(idern .le. neq)
        iprem = scbl(ibloc)+1
        if (typvar .eq. 1) then
            do i = iprem, idern
                zr(iadigs-1+i) = zr(iavale-1+zi(jsxdi+i-1))
            end do
        else if (typvar .eq. 2) then
            do i = iprem, idern
                zc(iadigs-1+i) = zc(iavale-1+zi(jsxdi+i-1))
            end do
        else
            ASSERT(.false.)
        end if
        call jelibe(jexnum(noma19//'.UALF', ibloc))
    end do
!
!
!
9998 continue
!     -- CALCUL DE EPS :
!     ------------------
!     ON AVAIT PENSE CALCULER EPS COMME:
!     1.D-15 FOIS LE TERME DIAGONAL MIN (/=0)
!     MAIS IL ARRIVE QU'AVEC MULT_FRONT ON PASSE EN
!     DESSOUS SANS QUE CELA FASSE D'OVERFLOW
!     DONC ON PREND UNE VALEUR ARBITRAIRE :
    eps = 1.d0/r8gaem()
!
    if (niv .gt. 1) then
        diamax = 0.d0
        diamin = r8maem()
        do i = 1, neq
            if (typvar .eq. 1) then
                vabs = abs(zr(iadigs-1+i))
            else
                vabs = abs(zc(iadigs-1+i))
            end if
            diamax = max(diamax, vabs)
            if (vabs .ne. 0.d0) diamin = min(diamin, vabs)
!
!
        end do
        write (ifm, *) '<FACTOR> AVANT FACTORISATION :'
        write (ifm, *) '<FACTOR>   NB EQUATIONS : ', neq
        write (ifm, *) '<FACTOR>   TERME DIAGONAL MAXIMUM :  ', diamax
        write (ifm, *) '<FACTOR>   TERME DIAGONAL (NON NUL) MINIMUM : ',&
     &                 diamin
        write (ifm, *) '<FACTOR>   EPSILON CHOISI  : ', eps
    end if
!
    call jedema()
end subroutine
