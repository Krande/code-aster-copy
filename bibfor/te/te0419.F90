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
subroutine te0419(option, nomte)
!
    implicit none
!
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/btldth.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/matrth.h"
#include "asterfort/moytem.h"
#include "asterfort/rcvarc.h"
#include "asterfort/trnflg.h"
#include "asterfort/vdxtemp.h"
#include "asterfort/vectan.h"
#include "asterfort/vexpan.h"
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
! Options: CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npge = 2
    real(kind=8), parameter :: epsval(npge) = (/-0.577350269189626d0, 0.577350269189626d0/)
    integer(kind=8) :: nb1, nb2, npgsn, npgsr, nddle
    integer(kind=8) :: i, j, ib, iret
    integer(kind=8) :: jvCacoqu, jvNbsp, jvMater, jvGeom, jvVect
    integer(kind=8) :: lzr, lzi
    integer(kind=8) :: kInf, kMoy, kSup
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3), vecpt(9, 3, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: forthi(42), forcth(42), vecl(51)
    real(kind=8) :: young, nu, alpha
    integer(kind=8) :: kpge, kpgsn, kpgsr, kwgt, nbLayer
    real(kind=8) :: epais, tempKpg, tempRefe, tempMoy
    real(kind=8) :: ksi3s2, ksi3
    integer(kind=8) :: elasID
    character(len=16) :: elasKeyword
    aster_logical :: hasTempRefe, hasTemp
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', jvGeom)
    call jevech('PVECTUR', 'E', jvVect)

! - Get objects
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    nddle = 5*nb1+2

! - Get reference temperature
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tempRefe, iret)
    hasTempRefe = iret .eq. 0

! - Get properties of shell
    call jevech('PNBSP_I', 'L', jvNbsp)
    nbLayer = zi(jvNbsp)
    kInf = 1
    kMoy = (3*nbLayer+1)/2
    kSup = 3*nbLayer

! - Get thickness
    call jevech('PCACOQU', 'L', jvCacoqu)
    epais = zr(jvCacoqu)

! - Compute local basis
    call vectan(nb1, nb2, zr(jvGeom), zr(lzr), vecta, &
                vectn, vectpt)

! - Get elasticity
    call jevech('PMATERC', 'L', jvMater)
    call get_elas_id(zi(jvMater), elasID, elasKeyword)
!
    kwgt = 0
    forcth = 0.d0
    do kpge = 1, npge
        ksi3 = epsval(kpge)
        ksi3s2 = epsval(kpge)/2.d0

! ----- MEMBRANE ET CISAILLEMENT
        do kpgsr = 1, npgsr
            call mahsms(0, nb1, zr(jvGeom), ksi3s2, kpgsr, &
                        zr(lzr), epais, vectn, vectg, vectt, &
                        hsfm, hss)
            call hsj1ms(epais, vectg, vectt, hsfm, hss, &
                        hsj1m, hsj1s)
            call btdmsr(nb1, nb2, ksi3s2, kpgsr, zr(lzr), &
                        epais, vectpt, hsj1m, hsj1s, btdm, &
                        btds)
        end do

        do kpgsn = 1, npgsn
! --------- FLEXION, CISAILLEMENT, NORMAL
            call mahsf(1, nb1, zr(jvGeom), ksi3s2, kpgsn, &
                       zr(lzr), epais, vectn, vectg, vectt, &
                       hsf)
            call hsj1f(kpgsn, zr(lzr), epais, vectg, vectt, &
                       hsf, kwgt, hsj1fx, wgt)
            call btdfn(1, nb1, nb2, ksi3s2, kpgsn, &
                       zr(lzr), epais, vectpt, hsj1fx, btdf)

! --------- Final btild
            call btdmsn(1, nb1, kpgsn, npgsr, zr(lzr), &
                        btdm, btdf, btds, btild)

! --------- Compute medium temperature
            call moytem('MASS', npgsn, 3*nbLayer, '+', tempMoy, iret)

! --------- Get temperature
            call vdxtemp(kInf, kMoy, kSup, &
                         kpgsn, ksi3, &
                         hasTempRefe, tempRefe, &
                         hasTemp, tempKpg)

! --------- Get elasticity parameters
            call matrth('MASS', &
                        elasID, elasKeyword, zi(jvMater), &
                        hasTemp, tempMoy, alpha, &
                        young, nu)

! --------- Compute "thermal" force
            call btldth(nb1, btild, wgt, &
                        hasTemp, tempKpg, &
                        young, nu, alpha, &
                        forthi)
            do i = 1, nddle
                forcth(i) = forcth(i)+forthi(i)
            end do
        end do
    end do
!
    call vexpan(nb1, forcth, vecl)
    do i = 1, 3
        vecl(6*nb1+i) = 0.d0
    end do
!
    do ib = 1, nb2
        do i = 1, 2
            do j = 1, 3
                vecpt(ib, i, j) = vectpt(ib, i, j)
            end do
        end do
        vecpt(ib, 3, 1) = vectn(ib, 1)
        vecpt(ib, 3, 2) = vectn(ib, 2)
        vecpt(ib, 3, 3) = vectn(ib, 3)
    end do
!
    call trnflg(nb2, vecpt, vecl, zr(jvVect))
!
end subroutine
