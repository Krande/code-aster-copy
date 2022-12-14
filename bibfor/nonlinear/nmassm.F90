! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmassm(lischa, numedd, numfix, typmat, optasz,&
                  meelem, matass)
!
use NonLin_Datastructure_type
!
implicit none
!
#include "jeveux.h"
#include "asterfort/asmaam.h"
#include "asterfort/asmama.h"
#include "asterfort/asmatr.h"
#include "asterfort/assert.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mtdscr.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/utmess.h"
!
character(len=19) :: lischa
character(len=24) :: numedd, numfix
character(len=6) :: typmat
character(len=*) :: optasz
character(len=19) :: meelem(8)
character(len=19) :: matass
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL)
!
! ASSEMBLAGE DES MATRICES ELEMENTAIRES
!
! --------------------------------------------------------------------------------------------------
!
! IN  LISCHA : LISTE DES CHARGEMENTS
! IN  OPTASS : OPTION D'ASSEMBLAGE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! IN  MEELEM : ARIABLE CHAPEAU POUR NOM DES MATR_ELEM
! OUT MATASS : MATR_ASSE CALCULEE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19) :: mediri, memass, meamor, messtr
    integer :: ifm, niv
    character(len=16) :: optass
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
!
! --- INITIALISATIONS
!
    optass = optasz
!
! --- DECOMPACTION DES VARIABLES CHAPEAUX
!
    if (meelem(1)(1:1) .ne. ' ') then
        call nmchex(meelem, 'MEELEM', 'MEDIRI', mediri)
        call nmchex(meelem, 'MEELEM', 'MEMASS', memass)
        call nmchex(meelem, 'MEELEM', 'MEAMOR', meamor)
        call nmchex(meelem, 'MEELEM', 'MESSTR', messtr)
    endif
!
! --- ASSEMBLAGE MATRICES ELEMENTAIRES
!
    if (typmat.eq.'MEAMOR') then
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_71')
        endif
        call asmaam(meamor, numedd, lischa, matass)
        call mtdscr(matass)
    else if (typmat.eq.'MEMASS') then
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_72')
        endif
        if (optass .eq. ' ') then
            call asmama(memass, ' ', numfix, lischa,&
                        matass)
        else if (optass.eq.'AVEC_DIRICHLET') then
            call asmama(memass, mediri, numedd, lischa,&
                        matass)
        endif
    else if (typmat.eq.'MESSTR') then
        if (niv .ge. 2) then
            call utmess('I', 'MECANONLINE13_73')
        endif
        call asmatr(1, messtr, ' ', numfix, &
                    lischa, 'ZERO', 'V', 1, matass)
        call mtdscr(matass)
    else if (typmat.eq.'MERIGI') then
! ----- Direct with asmari
        ASSERT(ASTER_FALSE)
    else
        ASSERT(ASTER_FALSE)
    endif
!
! - DEBUG
!
    if (niv .eq. 2) then
        call nmdebg('MATA', matass, ifm)
    endif
!
    call jedema()
!
end subroutine
