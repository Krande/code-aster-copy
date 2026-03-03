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
subroutine te0414(option, nomte)
!
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/cosiro.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matpgl.h"
#include "asterfort/tranlg.h"
#include "asterfort/utmess.h"
#include "asterfort/vdgnlr.h"
#include "asterfort/vdpnlr.h"
#include "asterfort/vdxnlr.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_3D
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "MASS"
    character(len=8), parameter :: typmod(2) = (/"C_PLAN  ", "        "/)
    integer(kind=8) :: nb1, jcret, codret
    real(kind=8) :: matloc(51, 51), plg(9, 3, 3)
    integer(kind=8) :: ibid, ideplm, ideplp, jvMaterc, jvCarcri
    integer(kind=8) :: jvGeom, jmatr, lzr, nb2, nddlet, lzi
    integer(kind=8) :: jvInstmr, jvInstpr
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: defoComp, relaComp
    aster_logical :: lVect, lMatr, lVari, lSigm
    type(Behaviour_Integ) :: BEHInteg
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb2 = zi(lzi-1+2)

! - Get input fields
    call cosiro(nomte, 'PCONTMR', 'L', 'UI', 'G', ibid, 'S')
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PDEPLMR', 'L', ideplm)
    call jevech('PDEPLPR', 'L', ideplp)
    call jevech('PINSTMR', 'L', jvInstmr)
    call jevech('PINSTPR', 'L', jvInstpr)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - No definition of local coordinate system
    call initLCSNone(materPara)

! - Get fields for behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', jvCarcri)

! - Get parameters for behaviour
    relaComp = compor(RELA_NAME)
    defoComp = compor(DEFO)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              zr(jvInstmr), zr(jvInstpr), &
                              materPara, BEHInteg)

! - Get output fields
    if (lMatr) then
        call jevech('PMATUUR', 'E', jmatr)
    end if
    if (lSigm) then
        call jevech('PCODRET', 'E', jcret)
    end if

! - Some checks
    if (relaComp(1:5) .eq. 'ELAS_') then
        call utmess('F', 'PLATE1_12', sk=relaComp)
    end if

! - Compute
    if (defoComp .eq. 'GROT_GDEP') then
        if (relaComp .eq. 'ELAS ') then
            call vdgnlr(materPara, &
                        lMatr, lVect, lSigm, lVari, relaComp, nomte)
            codret = 0
        else
            call vdpnlr(BEHInteg, option, nomte, codret)
        end if
    else if (defoComp(1:5) .eq. 'PETIT') then
        call vdxnlr(BEHInteg, &
                    option, nomte, zr(jvGeom), matloc, nb1, &
                    codret)
        if (lMatr) then
! -----    MATRICE DE PASSAGE REPERE GLOBAL REPERE LOCAL
            call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
            call matpgl(nb2, zr(lzr), plg)
! -----    OPERATION DE TRANFORMATION DE MATLOC DANS LE REPERE GLOBAL ET STOCKAGE DANS ZR
            nddlet = 6*nb1+3
            call tranlg(nb1, 51, nddlet, plg, matloc, zr(jmatr))
        end if
    else
        call utmess('F', 'PLATE1_14', sk=defoComp)
    end if
!
    if (lSigm) then
        zi(jcret) = codret
    end if
!
    if (lSigm) then
        call cosiro(nomte, 'PCONTPR', 'E', 'IU', 'G', ibid, 'R')
    end if
!
end subroutine
