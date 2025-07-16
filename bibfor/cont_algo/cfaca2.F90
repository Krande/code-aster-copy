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

subroutine cfaca2(nbliac, spliai, &
                  indfac, sdcont_solv, lmat, &
                  xjvmax)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/caladu.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer(kind=8) :: nbliac
    integer(kind=8) :: spliai, indfac
    integer(kind=8) :: lmat
    real(kind=8) :: xjvmax
    character(len=24) :: sdcont_solv
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Solve
!
! Discrete methods - Compute A.C-1.AT
!
! --------------------------------------------------------------------------------------------------
!
! IN  RESOCO : SD DE RESOLUTION DU CONTACT
! IN  LMAT   : DESCRIPTEUR DE LA MATR_ASSE DU SYSTEME MECANIQUE
! IN  NBLIAC : NOMBRE DE LIAISONS ACTIVES
! I/O XJVMAX : VALEUR DU PIVOT MAX
! I/O SPLIAI : INDICE DANS LA LISTE DES LIAISONS ACTIVES DE LA DERNIERE
!              LIAISON AYANT ETE CALCULEE POUR LE VECTEUR CM1A
! I/O INDFAC : INDICE DE DEBUT DE LA FACTORISATION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jdecal
    integer(kind=8) :: nbddl, jva, jvale, neq
    integer(kind=8) :: iliac, jj, lliac, lljac, ii, dercol, bloc
    integer(kind=8) :: nbbloc
    real(kind=8) :: val
    character(len=19) :: liac, cm1a
    integer(kind=8) :: jliac, jcm1a
    character(len=19) :: stoc
    character(len=19) :: ouvert, macont
    integer(kind=8) :: jouv
    character(len=24) :: appoin, apddl, apcoef
    integer(kind=8) :: japptr, japddl, japcoe
    integer(kind=8), pointer :: scbl(:) => null()
    integer(kind=8), pointer :: scde(:) => null()
    integer(kind=8), pointer :: scib(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! --- RECUPERATION D'OBJETS JEVEUX
!
    cm1a = sdcont_solv(1:14)//'.CM1A'
    appoin = sdcont_solv(1:14)//'.APPOIN'
    apddl = sdcont_solv(1:14)//'.APDDL'
    liac = sdcont_solv(1:14)//'.LIAC'
    apcoef = sdcont_solv(1:14)//'.APCOEF'
    macont = sdcont_solv(1:14)//'.MATC'
    stoc = sdcont_solv(1:14)//'.SLCS'
    call jeveuo(appoin, 'L', japptr)
    call jeveuo(apddl, 'L', japddl)
    call jeveuo(liac, 'L', jliac)
    call jeveuo(apcoef, 'L', japcoe)
    call jeveuo(stoc//'.SCIB', 'L', vi=scib)
    call jeveuo(stoc//'.SCBL', 'L', vi=scbl)
    call jeveuo(stoc//'.SCDE', 'L', vi=scde)
!
! --- INITIALISATIONS
!
    neq = zi(lmat+2)
    nbbloc = scde(3)
    ouvert = '&CFACA2.TRAV'
!
    call wkvect(ouvert, 'V V L', nbbloc, jouv)

! ======================================================================
! --- CALCUL DE -A.C-1.AT (REDUITE AUX LIAISONS ACTIVES) ---------------
! --- (STOCKAGE DE LA MOITIE PAR SYMETRIE) -----------------------------
! ======================================================================
    indfac = min(indfac, spliai+1)
    do iliac = spliai+1, nbliac
        lliac = zi(jliac-1+iliac)
        call jeveuo(jexnum(cm1a, lliac), 'L', jcm1a)
        ii = scib(iliac)
        dercol = scbl(ii)
        bloc = dercol*(dercol+1)/2
        if (.not. zl(jouv-1+ii)) then
            if (ii .gt. 1) then
                call jelibe(jexnum(macont//'.UALF', (ii-1)))
                zl(jouv-2+ii) = .false.
            end if
            call jeveuo(jexnum(macont//'.UALF', ii), 'E', jvale)
            zl(jouv-1+ii) = .true.
        end if
        jva = jvale-1+(iliac-1)*(iliac)/2-bloc
        do jj = 1, iliac
            lljac = zi(jliac-1+jj)
            jdecal = zi(japptr+lljac-1)
            nbddl = zi(japptr+lljac)-zi(japptr+lljac-1)
            jva = jva+1
            zr(jva) = 0.0d0
            call caladu(neq, nbddl, zr(japcoe+jdecal), zi(japddl+jdecal), zr(jcm1a), &
                        val)
            zr(jva) = zr(jva)-val
            if (abs(zr(jva)) .gt. xjvmax) xjvmax = abs(zr(jva))
        end do
        call jelibe(jexnum(cm1a, lliac))
    end do
!
! ======================================================================
    spliai = nbliac
    call jedetr(ouvert)
    call jedema()
!
end subroutine
