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

subroutine cpdepl(melflu, base, nuor, nbm)
    implicit none
!  RECOPIE LES CHAMPS DE DEPLACEMENTS PRIS DANS UN CONCEPT MODE_MECA
!  LES DDL DE LAGRANGE SONT ELIMINES
!  APPELANT : FLUST1 , FLUST2
!-----------------------------------------------------------------------
!  IN : MELFLU : NOM DU CONCEPT DE TYPE MELASFLU PRODUIT
!  IN : BASE   : NOM DU CONCEPT DE TYPE MODE_MECA DEFINISSANT LA BASE
!                MODALE DU SYSTEME AVANT PRISE EN COMPTE DU COUPLAGE
!  IN : NUOR   : LISTE DES NUMEROS D'ORDRE DES MODES SELECTIONNES POUR
!                LE COUPLAGE (SUR LESQUELS PORTE L'EXTRACTION)
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/extmod_sorted.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbm, nuor(nbm)
    character(len=8) :: base
    character(len=19) :: melflu
!
    integer(kind=8) :: iddl(6)
    character(len=8) :: mailla
    character(len=14) :: numddl
    character(len=24) :: nomcha, matria, nomnoe
!-----------------------------------------------------------------------
    integer(kind=8) :: icham, im, imod, lnoe
    integer(kind=8) :: neq
!-----------------------------------------------------------------------
    data iddl/1, 2, 3, 4, 5, 6/
!
!-----------------------------------------------------------------------
    call jemarq()
!
    nomcha(1:13) = melflu(1:8)//'.C01.'
    nomcha(17:24) = '001.VALE'
!
    call wkvect('&&CPDEPL.TEMP.NUOR', 'V V I', 1, imod)
!
!
    call dismoi('REF_RIGI_PREM', base, 'RESU_DYNA', repk=matria)
!
    call dismoi('NOM_NUME_DDL', matria, 'MATR_ASSE', repk=numddl)
    call dismoi('NB_EQUA', matria, 'MATR_ASSE', repi=neq)
    call dismoi('NOM_MAILLA', matria, 'MATR_ASSE', repk=mailla)
    nomnoe = mailla//'.COORDO    .VALE'
    call jelira(nomnoe, 'LONMAX', lnoe)
    lnoe = lnoe/3
!
    do im = 1, nbm
        write (nomcha(14:16), '(I3.3)') nuor(im)
        call jeveuo(nomcha, 'E', icham)
        zi(imod) = nuor(im)
        call extmod_sorted(base, numddl, zi(imod), 1, zr(icham), &
                           neq, lnoe, iddl, 6)
        call jelibe(nomcha)
    end do
!
!     MENAGE
    call jedetr('&&CPDEPL.TEMP.NUOR')
!
    call jedema()
!
end subroutine
