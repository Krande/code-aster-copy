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
!
subroutine te0121(option, nomte)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/pmfmats.h"
#include "asterfort/lteatt.h"
!
character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: all except DIS_*
!
! Option: AMOR_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer, parameter :: nbPara = 3
    character(len=8), parameter :: paraName(nbPara) = (/'X', 'Y', 'Z'/)
    real(kind=8) :: paraVale(nbPara)
    integer, parameter :: nbResu = 2
    character(len=16), parameter :: resuName(nbResu) = (/'AMOR_ALPHA',&
                                                         'AMOR_BETA '/)
    real(kind=8) :: resuVale(nbResu)
    integer :: icodre(nbResu)
    integer :: iret, nbddl, geomDime, nbParaEff
    integer :: i, j, iGeomDime, kns, ks
    integer :: jvMate, jvMateCod, jvGeom
    integer :: nbNode, iNode
    real(kind=8) :: alpha, beta, vxyz
    integer :: matrRigiSize, matrResuSize, matrMassSize
    aster_logical :: lAbso
    aster_logical :: lRigiSyme, lMatrRigi, lMatrMass, lMatrResuSyme
    character(len=8) ::  materPMF
    character(len=32) :: elasKeyword
    integer, parameter :: tecachNbVal = 5
    integer :: itab(tecachNbVal)
    integer :: jvMass, jvRigi, jvResu
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'AMOR_MECA')
!
    call elrefe_info(fami='RIGI', nno=nbNode)
    lAbso = lteatt('ABSO', 'OUI')

! - Get RIGI_MECA (symmetric)
    call tecach('NNO', 'PRIGIEL', 'L', iret, nval=tecachNbVal, itab=itab)
    lMatrRigi = iret .eq. 0
    lRigiSyme = ASTER_FALSE
    jvRigi = 0
    matrRigiSize = 0
    if (lMatrRigi) then
        call tecach('OOO', 'PRIGIEL', 'L', iret, nval=tecachNbVal, itab=itab)
        lRigiSyme = ASTER_TRUE
        jvRigi = itab(1)
        matrRigiSize = itab(2)
    endif

! - Get MASS_MECA
    call tecach('ONO', 'PMASSEL', 'L', iret, nval=tecachNbVal, itab=itab)
    lMatrMass = iret .eq. 0
    jvMass = 0
    matrMassSize = 0
    if (lMatrMass) then
        call tecach('OOO', 'PMASSEL', 'L', iret, nval=tecachNbVal, itab=itab)
        jvMass = itab(1)
        matrMassSize = itab(2)
    endif

! - Select output matrix (symmetric or not)
    lMatrResuSyme = ASTER_TRUE
    if (lRigiSyme) then
        call tecach('ONO', 'PMATUUR', 'E', iret, nval=tecachNbVal, itab=itab)
        lMatrResuSyme = ASTER_TRUE
    else
        call tecach('NNO', 'PMATUNS', 'E', iret, nval=tecachNbVal, itab=itab)
        lMatrResuSyme = ASTER_FALSE
        if (iret .ne. 0) then
            lMatrResuSyme = ASTER_TRUE
            call tecach('ONO', 'PMATUUR', 'E', iret,&
                        nval=tecachNbVal, itab=itab)
        endif
    endif
    jvResu = itab(1)
    matrResuSize = itab(2)

    if (.not.lMatrRigi .and. .not. lMatrResuSyme) then
        call tecach('OOO', 'PRIGINS', 'L', iret, nval=tecachNbVal, itab=itab)
        lMatrRigi = ASTER_TRUE
        jvRigi = itab(1)
        matrRigiSize = itab(2)
    endif

! - Access to geometry
    call tecach('ONO', 'PGEOMER', 'L', iret, tecachNbVal, itab=itab)
    jvGeom = itab(1)
    geomDime = itab(2)/nbNode
    ASSERT(geomDime .eq. 2 .or. geomDime .eq. 3)

! - Get (average) coordinates of cell
    nbParaEff = geomDime
    do iGeomDime = 1, geomDime
        vxyz = 0.d0
        do iNode = 1, nbNode
            vxyz = vxyz+zr(jvGeom + geomDime*(iNode-1)+iGeomDime -1)
        enddo
        paraVale(iGeomDime) = vxyz/nbNode
    enddo

! - Get material parameters
    call jevech('PMATERC', 'L', jvMate)
    jvMateCod = zi(jvMate)
    call rccoma(jvMateCod, 'ELAS', 1, elasKeyword, icodre(1))
    if(.not.(elasKeyword .eq. 'ELAS'       .or. elasKeyword .eq. 'ELAS_COQMU' .or. &
             elasKeyword .eq. 'ELAS_GLRC'  .or. elasKeyword .eq. 'ELAS_DHRC'  .or. &
             elasKeyword .eq. 'ELAS_ORTH')) then
        call utmess('F', 'PLATE1_1', nk=2, valk=[option, elasKeyword])
    endif

!   si l'élément est multifibre, il faut prendre le materiau "section"
!   pour récupérer les coefficients de dilatation :
    call pmfmats(jvMateCod, materPMF)
    resuVale = 0.d0
    call rcvalb('RIGI', 1, 1, '+', jvMateCod, materPMF, elasKeyword,&
                nbParaEff, paraName, paraVale,&
                nbResu, resuName, resuVale,&
                icodre, 0, nan='NON')
    alpha = resuVale(1)
    beta = resuVale(2)

! - Some checks
    if (lMatrResuSyme) then
        if (lMatrMass) then
            ASSERT(matrMassSize .eq. matrResuSize)
        endif
    else
        nbddl = int(sqrt(dble(matrResuSize)))
        if (lMatrMass) then
            ASSERT(matrMassSize .eq. nbddl*(nbddl+1)/2)
        endif
    endif
    if (lMatrRigi) then
        ASSERT(matrRigiSize .eq. matrResuSize)
    endif
    ASSERT(jvRigi .ne. 0 .or. jvMass .ne. 0)
    if (jvMass .eq. 0 .and. .not. lAbso) then
        call utmess('F', "MATRICE0_12")
    endif

! - Compute
    if (lMatrResuSyme) then
        do i = 1, matrResuSize
            if (jvRigi .eq. 0) then
                zr(jvResu-1+i) = beta*zr(jvMass-1+i)
            else
                if (jvMass .eq. 0) then
                    zr(jvResu-1+i) = alpha*zr(jvRigi-1+i)
                else
                    zr(jvResu-1+i) = alpha*zr(jvRigi-1+i) + beta*zr(jvMass-1+i)
                endif
            endif
        enddo
    else
        ASSERT(jvRigi .ne. 0)
        ASSERT(nbddl .gt. 0)
        do i = 1, nbddl
            kns = (i-1)*nbddl
            do j = 1, nbddl
                if (j .le. i) then
                    ks = (i-1)*i/2+j
                else
                    ks = (j-1)*j/2+i
                endif
                if (jvMass .eq. 0) then
                    zr(jvResu-1+kns+j) = alpha*zr(jvRigi-1+kns+j)
                else
                    zr(jvResu-1+kns+j) = alpha*zr(jvRigi-1+kns+j) + beta*zr(jvMass-1+ks)
                endif
            enddo
        enddo
    endif
!
end subroutine
