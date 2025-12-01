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
subroutine nmacex(sddisc, iterNewt, lExtrapol, extrapolVale)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/nmdcrg.h"
#include "asterfort/nmlere.h"
#include "asterfort/nmlerr.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: iterNewt
    aster_logical, intent(out) :: lExtrapol
    real(kind=8), intent(out) :: extrapolVale(4)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - GESTION DES EVENEMENTS)
!
! EXTRAPOLATION LINEAIRE DES RESIDUS
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc          : datastructure for time discretization
! In  iterNewt        : index of current Newton iteration
! OUT LEXTRA : .TRUE. SI EXTRAPOLATION OK
! OUT VALEXT : VALEURS DE L'EXTRAPOLATION (XA0 + ITER*XA1) / XDET
!               VALEXT(1): XA0
!               VALEXT(2): XA1
!               VALEXT(3): XDET
!               VALEXT(4): CRESI (RESIDU CIBLE)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: regres, depart
    real(kind=8) :: cresi, crela, cmaxi
    real(kind=8) :: vrela(1), vmaxi(1)
    real(kind=8) :: xa0, xa1, xdet
    integer(kind=8) :: nbIter, minIter, maxIter
    integer(kind=8) :: nbigno
    real(kind=8), pointer :: erreurs(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    regres = 0
    lExtrapol = ASTER_FALSE
    extrapolVale = 0.d0

! - AFFICHAGE
    call utmess('I', 'EXTRAPOLATION_1')

! - LECTURE DES INFOS SUR LES CONVERGENCES
    call nmlerr(sddisc, 'MXITER', paraValeI_=maxIter)
    call nmlerr(sddisc, 'MNITER', paraValeI_=minIter)
    call nmlerr(sddisc, 'NBITER', paraValeI_=nbIter)
    call nmlerr(sddisc, 'RESI_GLOB_RELA', paraValeR_=crela)
    call nmlerr(sddisc, 'RESI_GLOB_MAXI', paraValeR_=cmaxi)

! - REGRESSION SUR GLOB_RELA OU GLOB_MAXI ?
    call nmlerr(sddisc, 'TYPE_RESI', paraValeI_=regres)
    call nmlere(sddisc, 'L', 'VRELA', iterNewt, vrela)
    call nmlere(sddisc, 'L', 'VMAXI', iterNewt, vmaxi)
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
        lExtrapol = ASTER_FALSE
        goto 999
    end if
!
! --- PARAMETRES DE LA METHODE D'EXTRAPOLATION
!
    nbigno = 3
!
! --- ASSEZ D'ITERATIONS POUR FAIRE L'EXTRAPOLATION ?
!
    if ((nbigno+3) .le. iterNewt) then
        depart = nbigno
    else
        lExtrapol = ASTER_FALSE
        call utmess('I', 'EXTRAPOLATION_3')
        goto 999
    end if
!
! --- TOUTES LES RESIDUS AU COURS DES ITERATIONS [0,ITERAT]
!
    AS_ALLOCATE(vr=erreurs, size=iterNewt+1)
    if (regres .eq. 1) then
        cresi = crela
        call nmlere(sddisc, 'L', 'VRELA_TOUS', iterNewt, erreurs)
    else if (regres .eq. 2) then
        cresi = cmaxi
        call nmlere(sddisc, 'L', 'VMAXI_TOUS', iterNewt, erreurs)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! --- CALCUL DE L'EXTRAPOLATION LINEAIRE
!
    call nmdcrg(depart, iterNewt, erreurs, xa0, xa1, &
                xdet)
    AS_DEALLOCATE(vr=erreurs)
!
! --- EXTRAPOLATION REUSSIE ?
!
    if (xdet .le. r8prem()) then
        call utmess('I', 'EXTRAPOLATION_10')
        lExtrapol = ASTER_FALSE
    else
        extrapolVale(1) = xa0
        extrapolVale(2) = xa1
        extrapolVale(3) = xdet
        extrapolVale(4) = cresi
        lExtrapol = ASTER_TRUE
    end if
!
999 continue
!
    call jedema()
end subroutine
