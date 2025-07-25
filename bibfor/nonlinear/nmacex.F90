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

subroutine nmacex(sddisc, iterat, lextra, valext)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmdcrg.h"
#include "asterfort/nmlere.h"
#include "asterfort/nmlerr.h"
#include "asterfort/utmess.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
    character(len=19) :: sddisc
    integer(kind=8) :: iterat
    aster_logical :: lextra
    real(kind=8) :: valext(4)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - GESTION DES EVENEMENTS)
!
! EXTRAPOLATION LINEAIRE DES RESIDUS
!
! ----------------------------------------------------------------------
!
!
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  ITERAT : NUMERO D'ITERATION DE NEWTON
! OUT LEXTRA : .TRUE. SI EXTRAPOLATION OK
! OUT VALEXT : VALEURS DE L'EXTRAPOLATION (XA0 + ITER*XA1) / XDET
!               VALEXT(1): XA0
!               VALEXT(2): XA1
!               VALEXT(3): XDET
!               VALEXT(4): CRESI (RESIDU CIBLE)
!
!
!
!
    integer(kind=8) :: ibid, regres, depart
    real(kind=8) :: cresi, crela, cmaxi
    real(kind=8) :: vrela(1), vmaxi(1)
    real(kind=8) :: r8bid
    real(kind=8) :: xa0, xa1, xdet
    integer(kind=8) :: nbiter, mniter, mxiter
    integer(kind=8) :: nbigno
    real(kind=8), pointer :: erreurs(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- INITIALISATIONS
!
    regres = 0
    lextra = .false.
    valext(1) = 0.d0
    valext(2) = 0.d0
    valext(3) = 0.d0
    valext(4) = 0.d0
!
! --- AFFICHAGE
!
    call utmess('I', 'EXTRAPOLATION_1')
!
! --- LECTURE DES INFOS SUR LES CONVERGENCES
!
    call nmlerr(sddisc, 'L', 'MXITER', r8bid, mxiter)
    call nmlerr(sddisc, 'L', 'MNITER', r8bid, mniter)
    call nmlerr(sddisc, 'L', 'NBITER', r8bid, nbiter)
    call nmlerr(sddisc, 'L', 'RESI_GLOB_RELA', crela, ibid)
    call nmlerr(sddisc, 'L', 'RESI_GLOB_MAXI', cmaxi, ibid)
!
! --- REGRESSION SUR GLOB_RELA OU GLOB_MAXI ?
!
    call nmlerr(sddisc, 'L', 'TYPE_RESI', r8bid, regres)
    call nmlere(sddisc, 'L', 'VRELA', iterat, vrela)
    call nmlere(sddisc, 'L', 'VMAXI', iterat, vmaxi)
!
! --- SI REGRES=3 ON DOIT FAIRE LA REGRESSION SUR LES 2, MAIS ON
! --- COMMENCE PAR LA FAIRE SUR GLOB_RELA
! --- SI VRELA > RGRELA ON MET REGRES=1 POUR FAIRE LA REGRESSION
! ---  SUR GLOB_RELA
! --- SI VMAXI > RGMAXI ON MET REGRES=2 POUR FAIRE LA REGRESSION
! ---  SUR GLOB_MAXI
!
    if (regres .eq. 3) then
        regres = 0
        if (vrela(1) .gt. crela) then
            regres = 1
        else if (vmaxi(1) .gt. cmaxi) then
            regres = 2
        end if
    end if
!
! --- LES CRITERES D'ERREUR SONT OK, MAIS PAS DETECTE DANS LE
! --- STAT_NON_LINE -> NI GLOB_RELA, NI GLOB_MAXI !
!
    if (regres .eq. 0) then
        call utmess('A', 'EXTRAPOLATION_2')
        lextra = .false.
        goto 999
    end if
!
! --- PARAMETRES DE LA METHODE D'EXTRAPOLATION
!
    nbigno = 3
!
! --- ASSEZ D'ITERATIONS POUR FAIRE L'EXTRAPOLATION ?
!
    if ((nbigno+3) .le. iterat) then
        depart = nbigno
    else
        lextra = .false.
        call utmess('I', 'EXTRAPOLATION_3')
        goto 999
    end if
!
! --- TOUTES LES RESIDUS AU COURS DES ITERATIONS [0,ITERAT]
!
    AS_ALLOCATE(vr=erreurs, size=iterat+1)
    if (regres .eq. 1) then
        cresi = crela
        call nmlere(sddisc, 'L', 'VRELA_TOUS', iterat, erreurs)
    else if (regres .eq. 2) then
        cresi = cmaxi
        call nmlere(sddisc, 'L', 'VMAXI_TOUS', iterat, erreurs)
    else
        ASSERT(.false.)
    end if
!
! --- CALCUL DE L'EXTRAPOLATION LINEAIRE
!
    call nmdcrg(depart, iterat, erreurs, xa0, xa1, &
                xdet)
    AS_DEALLOCATE(vr=erreurs)
!
! --- EXTRAPOLATION REUSSIE ?
!
    if (xdet .le. r8prem()) then
        call utmess('I', 'EXTRAPOLATION_10')
        lextra = .false.
    else
        valext(1) = xa0
        valext(2) = xa1
        valext(3) = xdet
        valext(4) = cresi
        lextra = .true.
    end if
!
999 continue
!
    call jedema()
end subroutine
