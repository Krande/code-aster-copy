! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
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
    use MaterialPara_module
    use MaterialPara_type
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
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DES COEFFICIENTS A0 ET A1 POUR LE PILOTAGE PAR CRITERE ELASTIQUE
!  OU PAR INCREMENT DE DEFORMATION POUR LES ELEMENTS A VARIABLES LOCALES
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    character(len=8) :: typmod(2)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relaComp, pilo
    integer(kind=8) :: ndim, nno, npg, lgpg, jtab(7), itype
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom, jvMaterc, jvCarcri
    integer(kind=8) :: icontm, ivarim, icopil, iborne, ictau
    integer(kind=8) :: ideplm, iddepl, idepl0, idepl1, iret
    real(kind=8) :: instam, instap
    type(Material_Para) :: materPara
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!

! - TYPE DE MODELISATION
    typmod = " "
    if (lteatt('DIM_TOPO_MODELI', '3')) then
        typmod(1) = '3D'
    else if (lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (lteatt('C_PLAN', 'OUI')) then
        typmod(1) = 'C_PLAN'
    else if (lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    end if
    typmod(2) = 'DEPLA'

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, &
                     ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(nno .le. MT_NNOMAX)
    ASSERT(npg .le. 27)

! - PARAMETRES EN ENTREE
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PCONTMR', 'L', icontm)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PDDEPLR', 'L', iddepl)
    call jevech('PDEPL0R', 'L', idepl0)
    call jevech('PDEPL1R', 'L', idepl1)

! - Continuation method: no time !
    instam = r8vide()
    instap = r8vide()

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Get fields for non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Properties of behaviour
    relaComp = compor(RELA_NAME)

! - Continuation method: no time !
    instam = r8vide()
    instap = r8vide()

! - Type of continuation
    call jevech('PTYPEPI', 'L', itype)
    pilo = zk16(itype)
    if (pilo .eq. 'PRED_ELAS') then
        call jevech('PCDTAU', 'L', ictau)
        call jevech('PBORNPI', 'L', iborne)
    end if
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              instam, instap, &
                              materPara, BEHInteg)

! - Prepare external state variables (geometry)
    if (relaComp .eq. 'BETON_DOUBLE_DP') then
        call behaviourPrepESVAGeom(nno, npg, ndim, &
                                   ipoids, ivf, idfde, &
                                   zr(jvGeom), BEHInteg)
    end if

! - Output field
    call jevech('PCOPILO', 'E', icopil)

! - Main subroutine to compute coefficients
    call pipepe(BEHInteg, &
                typmod, compor, &
                pilo, ndim, nno, npg, &
                ipoids, ivf, idfde, zr(jvGeom), &
                lgpg, zr(ideplm), zr(icontm), zr(ivarim), &
                zr(iddepl), zr(idepl0), zr(idepl1), zr(icopil), &
                iborne, ictau)
!
end subroutine
