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
subroutine vdxsig(nomte, nodeCoor, &
                  nbLayer, siefElga)
!
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
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
#include "asterfort/trndgl.h"
#include "asterfort/utmess.h"
#include "asterfort/vdefge.h"
#include "asterfort/vdesga.h"
#include "asterfort/vdxtemp.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: nodeCoor(3, 9)
    integer(kind=8), intent(in) :: nbLayer
    real(kind=8), intent(out) :: siefElga(6*27*nbLayer)
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, deux = 2.d0
! NOMBRE DE POINTS DE GAUSS DANS LA TRANCHE
    ! (POUR RESTER COHERENT AVEC SIEF_ELGA EN PLASTICITE )
    integer(kind=8), parameter :: npgLayer = 3
    real(kind=8), parameter :: epsval(3) = (/-1.d0, 0.d0, +1.d0/)
    integer(kind=8) :: nb1, nb2, npgsn, npgsr
    integer(kind=8) :: i, j, k
    integer(kind=8) :: jvCacoqu, jvDisp, jvNbsp, jvMater
    integer(kind=8) :: lzi, lzr
    integer(kind=8) :: iret
    integer(kind=8) :: iLayer, kpgLayer, kpgsn, kpgsr, kpgs, kwgt
    integer(kind=8) :: kInf, kMoy, kSup
    real(kind=8) :: siefKpg(6*27*nbLayer)
    real(kind=8) :: vecta(9, 2, 3), vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: disp(42), rotf(9)
    real(kind=8) :: alpha, epais
    real(kind=8) :: ksi3, ksi3s2, hic, zmin
    real(kind=8) :: tempRefe, tempKpg, tempMoy
    aster_logical :: hasTemp, hasTempRefe
    integer(kind=8) :: elasID
    character(len=16) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    siefElga = 0.d0

! - Get objects
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get thickness
    call jevech('PCACOQU', 'L', jvCacoqu)
    epais = zr(jvCacoqu)

! - Get properties of shell
    call jevech('PNBSP_I', 'L', jvNbsp)
    hic = un/nbLayer
    zmin = -0.5d0

! - Get reference temperature
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tempRefe, iret)
    hasTempRefe = iret .eq. 0

! - Get displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Compute local basis
    call vectan(nb1, nb2, nodeCoor, zr(lzr), vecta, &
                vectn, vectpt)
    call trndgl(nb2, vectn, vectpt, zr(jvDisp), disp, &
                rotf)

! - Get elasticity
    call jevech('PMATERC', 'L', jvMater)
    call get_elas_id(zi(jvMater), elasID, elasKeyword)
!
    kwgt = 0
    kpgs = 0
!
    do iLayer = 1, nbLayer
        kInf = 3*iLayer-2
        kMoy = 3*iLayer-1
        kSup = 3*iLayer
        do kpgLayer = 1, npgLayer
! --------- COTE DES POINTS D'INTEGRATION
            ksi3 = epsval(kpgLayer)
            if (kpgLayer .eq. 1) then
                ksi3s2 = zmin+(iLayer-1)*hic
            else if (kpgLayer .eq. 2) then
                ksi3s2 = zmin+hic/deux+(iLayer-1)*hic
            else
                ksi3s2 = zmin+hic+(iLayer-1)*hic
            end if

! --------- MEMBRANE ET CISAILLEMENT
            do kpgsr = 1, npgsr
                call mahsms(0, nb1, nodeCoor, ksi3s2, kpgsr, &
                            zr(lzr), epais, vectn, vectg, vectt, &
                            hsfm, hss)
                call hsj1ms(epais, vectg, vectt, hsfm, hss, &
                            hsj1m, hsj1s)
                call btdmsr(nb1, nb2, ksi3s2, kpgsr, zr(lzr), &
                            epais, vectpt, hsj1m, hsj1s, btdm, &
                            btds)
            end do

            do kpgsn = 1, npgsn
! ------------- FLEXION, CISAILLEMENT, NORMAL
                call mahsf(1, nb1, nodeCoor, ksi3s2, kpgsn, &
                           zr(lzr), epais, vectn, vectg, vectt, &
                           hsf)
                call hsj1f(kpgsn, zr(lzr), epais, vectg, vectt, &
                           hsf, kwgt, hsj1fx, wgt)
                call btdfn(1, nb1, nb2, ksi3s2, kpgsn, &
                           zr(lzr), epais, vectpt, hsj1fx, btdf)

! ------------- Final btild
                call btdmsn(1, nb1, kpgsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)

! ------------- Compute medium temperature
                call moytem('MASS', npgsn, 3*nbLayer, '+', tempMoy, iret)

! ------------- Get temperature
                call vdxtemp(kInf, kMoy, kSup, &
                             kpgsn, ksi3, &
                             hasTempRefe, tempRefe, &
                             hasTemp, tempKpg)

! ------------- Get dilatation coefficient
                if (hasTemp) then
                    call matrth('MASS', &
                                elasID, elasKeyword, zi(jvMater), &
                                hasTemp, tempMoy, alpha)
                else
                    alpha = r8nnem()
                end if

! ------------- Compute stress
                call vdesga(kwgt, nb1, nb2, &
                            vectt, disp, btild, &
                            hasTemp, alpha, tempKpg, &
                            siefKpg)
!
                kpgs = kpgs+1
                do i = 1, 6
                    siefElga(6*((kpgsn-1)*npgLayer*nbLayer+npgLayer*(iLayer-1)+kpgLayer-1)+i) = &
                        siefKpg(i+6*(kpgs-1))
                end do
            end do
        end do
    end do

! - DETERMINATION DES REPERES  LOCAUX DE L'ELEMENT AUX POINTS
! - D'INTEGRATION ET STOCKAGE DE CES REPERES DANS LE VECTEUR .DESR
    k = 0
    do kpgsr = 1, npgsr
        call vectgt(0, nb1, nodeCoor, zero, kpgsr, &
                    zr(lzr), epais, vectn, vectg, vectt)
        do j = 1, 3
            do i = 1, 3
                k = k+1
                zr(lzr+2000+k-1) = vectt(i, j)
            end do
        end do
    end do
!
end subroutine
