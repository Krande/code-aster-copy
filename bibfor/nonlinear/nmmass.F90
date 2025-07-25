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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine nmmass(lischa, sddyna, numedd, &
                  numfix, meelem, masse)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mtdscr.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/asmama.h"
!
    character(len=19) :: lischa, sddyna
    character(len=24) :: numedd, numfix
    character(len=19) :: meelem(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (CALCUL)
!
! PRE-CALCUL DES MATRICES ASSEMBLEES CONSTANTES AU COURS DU CALCUL
!
! ----------------------------------------------------------------------
!
!
! IN  FONACT : FONCTIONNALITES ACTIVEES (VOIR NMFONC)
! IN  LISCHA : LISTE DES CHARGEMENTS
! IN  SDDYNA : SD DYNAMIQUE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! In  ds_algopara      : datastructure for algorithm parameters
! IN  NUMFIX : NUME_DDL (FIXE AU COURS DU CALCUL)
! IN  MEELEM : MATRICES ELEMENTAIRES
! OUT MASSE  : MATRICE MASSE POUR LE CALCUL DES ENERGIES
!
! ----------------------------------------------------------------------
!
    aster_logical :: ldyna, lexpl, limpl
    character(len=16) :: optass
    character(len=19) :: masse, memass, mediri
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    ldyna = ndynlo(sddyna, 'DYNAMIQUE')
    lexpl = ndynlo(sddyna, 'EXPLICITE')
    limpl = ndynlo(sddyna, 'IMPLICITE')
    call nmchex(meelem, 'MEELEM', 'MEMASS', memass)
    call nmchex(meelem, 'MEELEM', 'MEDIRI', mediri)
!
! --- ASSEMBLAGE DE LA MATRICE MASSE
!
    if (ldyna) then
        if (limpl) then
            optass = ' '
        else if (lexpl) then
            optass = 'AVEC_DIRICHLET'
        else
            ASSERT(.false.)
        end if
        masse = '&&NMMASS.MASSENER'
        if (optass .eq. ' ') then
            call asmama(memass, ' ', numfix, lischa, masse)
        else if (optass .eq. 'AVEC_DIRICHLET') then
            call asmama(memass, mediri, numedd, lischa, masse)
        end if
        call mtdscr(masse)
    end if
!
    call jedema()
end subroutine
