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
subroutine te0033(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    use plateMaterial_module
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/cosiro.h"
#include "asterfort/dkqedg.h"
#include "asterfort/dkqsie.h"
#include "asterfort/dktedg.h"
#include "asterfort/dktsie.h"
#include "asterfort/dsqedg.h"
#include "asterfort/dsqsie.h"
#include "asterfort/dstedg.h"
#include "asterfort/dstsie.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxsiro.h"
#include "asterfort/dxsit2.h"
#include "asterfort/dxsit3.h"
#include "asterfort/dxsith.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/q4gedg.h"
#include "asterfort/q4gsie.h"
#include "asterfort/t3gedg.h"
#include "asterfort/t3gsie.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DE CONTRAINTES, DEFORMATIONS, EFFORTS ET DEFORMATIONS
!     GENERALISES POUR LES ELEMENTS DKT, DKTG, DST, DKQ, DSQ ET Q4G
!     POUR UN MATERIAU ISOTROPE OU MULTICOUCHE
!         OPTIONS TRAITEES  ==>  SIEF_ELGA
!                                EPSI_ELGA
!                                DEGE_ELGA
!                                DEGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: ipoids, ivf, idfdx
    integer(kind=8) :: jvCacoqu, jvDisp, jeffg, jvGeom, jvMaterc, jsigm
    integer(kind=8) :: np, multic, nbLayer
    real(kind=8) :: alpha, beta
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), r8bid
    real(kind=8) :: depl(24)
    real(kind=8) :: effgt(32), effpg(32)
    real(kind=8) :: t2iu(4), t2ui(4), c, s
    aster_logical :: lDKTG
    character(len=8) :: fami
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    r8bid = 0.d0
!
    if (option(6:9) .eq. 'ELNO') then
        fami = 'NOEU'
    else
        fami = 'RIGI'
    end if
    call elrefe_info(fami=fami, &
                     ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx)
!
    if (option .ne. 'SIEF_ELGA' .and. option .ne. 'EPSI_ELGA' .and. option .ne. 'DEGE_ELNO' &
        .and. option .ne. 'DEGE_ELGA') then
        ASSERT(ASTER_FALSE)
    end if
!
    lDKTG = ASTER_FALSE
    if ((nomte .eq. 'MEDKTG3') .or. (nomte .eq. 'MEDKQG4')) then
        lDKTG = ASTER_TRUE
    end if
!
    effgt = 0.d0

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    if (option .eq. 'SIEF_ELGA' .or. option .eq. 'EPSI_ELGA') then
! ----- Get material parameters
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), materPara)

! ----- Set local coordinate system from user
        call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)
    end if

! - For multi-layers
    nbLayer = 0
    if (option .eq. 'SIEF_ELGA' .or. option .eq. 'EPSI_ELGA') then
        call getMultiLayerNbLayer(materPara, lDKTG, nbLayer)
    end if
!
    if (option(8:9) .eq. 'GA') then
        np = npg
    else if (option(8:9) .eq. 'NO') then
        np = nno
    end if

! - Management of local coordinate system (intrinsec, for plate)
    if (nno .eq. 3) then
        call dxtpgl(zr(jvGeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(jvGeom), pgl)
    end if
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)
    call jevech('PCACOQU', 'L', jvCacoqu)
    alpha = zr(jvCacoqu+1)*r8dgrd()
    beta = zr(jvCacoqu+2)*r8dgrd()
    call coqrep(pgl, alpha, beta, t2iu, t2ui, c, s)
!
    call jevech('PDEPLAR', 'L', jvDisp)
    call utpvgl(nno, 6, pgl, zr(jvDisp), depl)
!
    if (option(1:9) .eq. 'SIEF_ELGA') then
        call jevech('PCONTRR', 'E', jsigm)
        ASSERT(nbLayer .ge. 1)
        if (nomte .eq. 'MEDKTR3') then
            call dktsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEDSTR3') then
            call dstsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEDKQU4') then
            call dkqsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEDSQU4') then
            call dsqsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MET3TR3') then
            call t3gsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (materPara%elasID .eq. ELAS_ISOT .or. &
            materPara%elasID .eq. ELAS_ORTH .or. &
            materPara%elasID .eq. ELAS_ISTR) then
            call dxsith(nomte, materPara, zr(jsigm))

        else if (materPara%elasID .eq. ELAS_COMPOSITE) then
            call dxsit2(nomte, pgl, zr(jsigm))

        elseif (materPara%elasID .eq. ELAS_SHELL) then
            call dxsit3(nomte, zi(jvMaterc), pgl, zr(jsigm))

        else
            call utmess('F', 'PLATE1_1', nk=2, valk=[option, materPara%elasKeyword])
        end if

    else if (option(1:9) .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jsigm)
        ASSERT(nbLayer .ge. 1)
        if (nomte .eq. 'MEDKTR3') then
            call dktsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEDSTR3') then
            call dstsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEDKQU4') then
            call dkqsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEDSQU4') then
            call dsqsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MEQ4QU4') then
            call q4gsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        else if (nomte .eq. 'MET3TR3') then
            call t3gsie(option, fami, xyzl, pgl, depl, &
                        nbLayer, zr(jsigm))
        end if
        call dxsiro(np*nbLayer*3, t2iu, zr(jsigm), zr(jsigm))

    else if (option(1:9) .eq. 'DEGE_ELNO') then
        call jevech('PDEFOGR', 'E', jeffg)
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(xyzl, option, pgl, depl, effgt, multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(xyzl, option, pgl, depl, effgt)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(xyzl, option, pgl, depl, effgt)
        end if
        call dxefro(np, t2iu, effgt, zr(jeffg))

    else if (option(1:9) .eq. 'DEGE_ELGA') then
        call jevech('PDEFOPG', 'E', jeffg)
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(xyzl, option, pgl, depl, effpg, multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(xyzl, option, pgl, depl, effpg)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(xyzl, option, pgl, depl, effpg)
        end if
        call dxefro(np, t2iu, effpg, zr(jeffg))

    end if
!
    if (option .eq. 'SIEF_ELGA') then
        call cosiro(nomte, 'PCONTRR', 'E', 'IU', 'G', &
                    jsigm, 'S')
    end if
!
end subroutine
