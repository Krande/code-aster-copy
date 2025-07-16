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
subroutine te0543(option, nomte)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/pipepe.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES COEFFICIENTS A0 ET A1
!                          POUR LE PILOTAGE PAR CRITERE ELASTIQUE
!                          OU PAR INCREMENT DE DEFORMATION POUR LES
!                          ELEMENTS A VARIABLES LOCALES
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    character(len=8) :: typmod(2)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: rela_comp, pilo
    character(len=4), parameter :: fami = 'RIGI'
    integer(kind=8) :: jgano, ndim, nno, nnos, npg, lgpg, jtab(7), itype
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate, icarcr
    integer(kind=8) :: icontm, ivarim, icopil, iborne, ictau
    integer(kind=8) :: ideplm, iddepl, idepl0, idepl1, iret
    real(kind=8) :: instam, instap
    type(Behaviour_Integ) :: BEHinteg

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)
!
! - TYPE DE MODELISATION
!
    if (lteatt('DIM_TOPO_MODELI', '3')) then
        typmod(1) = '3D'
    else if (lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (lteatt('C_PLAN', 'OUI')) then
        typmod(1) = 'C_PLAN'
    else if (lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    end if
!
    typmod(2) = 'DEPLA'
!
! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    ASSERT(nno .le. 27)
    ASSERT(npg .le. 27)

! - PARAMETRES EN ENTREE
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDDEPLR', 'L', iddepl)
    call jevech('PDEPL0R', 'L', idepl0)
    call jevech('PDEPL1R', 'L', idepl1)
    call jevech('PTYPEPI', 'L', itype)
    call jevech('PCARCRI', 'L', icarcr)
!
    pilo = zk16(itype)
    rela_comp = compor(RELA_NAME)
    if (pilo .eq. 'PRED_ELAS') then
        call jevech('PCDTAU', 'L', ictau)
        call jevech('PBORNPI', 'L', iborne)
    end if
!
! -- NOMBRE DE VARIABLES INTERNES
!
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Continuation method: no time !
    instam = r8vide()
    instap = r8vide()

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndim, typmod, option, &
                              compor, zr(icarcr), &
                              instam, instap, &
                              fami, zi(imate), &
                              BEHinteg)

! ----- Prepare external state variables (geometry)
    if (rela_comp .eq. 'BETON_DOUBLE_DP') then
        call behaviourPrepESVAGeom(nno, npg, ndim, &
                                   ipoids, ivf, idfde, &
                                   zr(igeom), BEHinteg)

    end if
!
! PARAMETRES EN SORTIE
!
    call jevech('PCOPILO', 'E', icopil)
!
    call pipepe(BEHinteg, &
                pilo, ndim, nno, npg, ipoids, &
                ivf, idfde, zr(igeom), typmod, zi(imate), &
                compor, lgpg, zr(ideplm), zr(icontm), zr(ivarim), &
                zr(iddepl), zr(idepl0), zr(idepl1), zr(icopil), &
                iborne, ictau)
!
end subroutine
