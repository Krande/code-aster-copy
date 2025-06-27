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
subroutine cachei(load, model, mesh, valeType, param, keywordFactZ)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: load, mesh, model
    character(len=4), intent(in) :: valeType
    character(len=5), intent(in) :: param
    character(len=*), intent(in) :: keywordFactZ
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads 'PRE_EPSI'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, nchei, ncmp, jvale, jvalv, iocc, nxx, nyy, nzz
    integer(kind=8) :: nxy, nxz, nyz, nex, nky, nkz, nexx, neyy, nexy, nkxx, nkyy, nkxy
    integer(kind=8) :: nbtou, nbma, jma, nepsi, nvectn, nkn1, nkn2
    real(kind=8) :: epxx, epyy, epzz, epxy, epxz, epyz, epx, xky, xkz, xexx
    real(kind=8) :: xeyy, xexy, xkxx, xkyy, xkxy, vect_n(3), xkn1, xkn2
    character(len=8) :: k8b, kepxx, kepyy, kepzz, kepxy, kepxz, kepyz
    character(len=8) :: kepx, kxky, kxkz, kxexx, kxeyy, kxexy, kxkxx, kxkyy, kxkxy
    character(len=8) :: typmcl(2), kxkn1, kxkn2, kxvn1, kxvn2, kxvn3
    character(len=16) :: keywordFact, motcle(2)
    character(len=19) :: carte
    character(len=24) :: mesmai, chepsi
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    keywordFact = keywordFactZ
    call getfac(keywordFact, nchei)
!
    carte = load//'.CHME.'//param
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'EPSI_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'EPSI_F')
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
    call jeveuo(carte//'.VALE', 'E', jvale)
!
    ncmp = 18
!
    vncmp(1) = 'EPXX'
    vncmp(2) = 'EPYY'
    vncmp(3) = 'EPZZ'
    vncmp(4) = 'EPXY'
    vncmp(5) = 'EPXZ'
    vncmp(6) = 'EPYZ'
    vncmp(7) = 'EPX'
    vncmp(8) = 'KY'
    vncmp(9) = 'KZ'
    vncmp(10) = 'EXX'
    vncmp(11) = 'EYY'
    vncmp(12) = 'EXY'
    vncmp(13) = 'KXX'
    vncmp(14) = 'KYY'
    vncmp(15) = 'KXY'
!
    vncmp(16) = 'VECT_2_X'
    vncmp(17) = 'VECT_2_Y'
    vncmp(18) = 'VECT_2_Z'
    if (valeType .eq. 'REEL') then
        do i = 1, ncmp
            zr(jvalv-1+i) = 0.d0
        end do
    else
        do i = 1, ncmp
            zk8(jvalv-1+i) = '&FOZERO'
        end do
    end if
    call nocart(carte, 1, ncmp)
!
    mesmai = '&&CACHEI.MES_MAILLES'
    motcle(1) = 'GROUP_MA'
    motcle(2) = 'MAILLE'
    typmcl(1) = 'GROUP_MA'
    typmcl(2) = 'MAILLE'
!
    do iocc = 1, nchei
        if (valeType .eq. 'REEL') then

            call getvid(keywordFact, 'EPSI', iocc=iocc, scal=chepsi, nbret=nepsi)
            if (nepsi .ne. 0) call utmess('F', 'CHARGES_5')

            call getvr8(keywordFact, 'EPXX', iocc=iocc, scal=epxx, nbret=nxx)
            call getvr8(keywordFact, 'EPYY', iocc=iocc, scal=epyy, nbret=nyy)
            call getvr8(keywordFact, 'EPZZ', iocc=iocc, scal=epzz, nbret=nzz)
            call getvr8(keywordFact, 'EPXY', iocc=iocc, scal=epxy, nbret=nxy)
            call getvr8(keywordFact, 'EPXZ', iocc=iocc, scal=epxz, nbret=nxz)
            call getvr8(keywordFact, 'EPYZ', iocc=iocc, scal=epyz, nbret=nyz)
            call getvr8(keywordFact, 'EPX', iocc=iocc, scal=epx, nbret=nex)
            call getvr8(keywordFact, 'KY', iocc=iocc, scal=xky, nbret=nky)
            call getvr8(keywordFact, 'KZ', iocc=iocc, scal=xkz, nbret=nkz)
            call getvr8(keywordFact, 'EXX', iocc=iocc, scal=xexx, nbret=nexx)
            call getvr8(keywordFact, 'EYY', iocc=iocc, scal=xeyy, nbret=neyy)
            call getvr8(keywordFact, 'EXY', iocc=iocc, scal=xexy, nbret=nexy)
            call getvr8(keywordFact, 'KXX', iocc=iocc, scal=xkxx, nbret=nkxx)
            call getvr8(keywordFact, 'KYY', iocc=iocc, scal=xkyy, nbret=nkyy)
            call getvr8(keywordFact, 'KXY', iocc=iocc, scal=xkxy, nbret=nkxy)
            call getvr8(keywordFact, 'KN1', iocc=iocc, scal=xkn1, nbret=nkn1)
            call getvr8(keywordFact, 'KN2', iocc=iocc, scal=xkn2, nbret=nkn2)
            call getvr8(keywordFact, 'VECT_N', iocc=iocc, nbval=3, vect=vect_n, nbret=nvectn)
!
            do i = 1, ncmp
                zr(jvalv-1+i) = 0.d0
            end do
!
            if (nxx .ne. 0) zr(jvalv-1+1) = epxx
            if (nyy .ne. 0) zr(jvalv-1+2) = epyy
            if (nzz .ne. 0) zr(jvalv-1+3) = epzz
            if (nxy .ne. 0) zr(jvalv-1+4) = epxy
            if (nxz .ne. 0) zr(jvalv-1+5) = epxz
            if (nyz .ne. 0) zr(jvalv-1+6) = epyz
!
            if (nex .ne. 0) zr(jvalv-1+7) = epx
            if (nky .ne. 0) zr(jvalv-1+8) = xky
            if (nkz .ne. 0) zr(jvalv-1+9) = xkz
            if (nkn1 .ne. 0) zr(jvalv-1+8) = xkn1
            if (nkn2 .ne. 0) zr(jvalv-1+9) = xkn2
            if (nexx .ne. 0) zr(jvalv-1+10) = xexx
            if (neyy .ne. 0) zr(jvalv-1+11) = xeyy
            if (nexy .ne. 0) zr(jvalv-1+12) = xexy
            if (nkxx .ne. 0) zr(jvalv-1+13) = xkxx
            if (nkyy .ne. 0) zr(jvalv-1+14) = xkyy
            if (nkxy .ne. 0) zr(jvalv-1+15) = xkxy
            if (nvectn .ne. 0) then
                zr(jvalv-1+16) = vect_n(1)
                zr(jvalv-1+17) = vect_n(2)
                zr(jvalv-1+18) = vect_n(3)
            end if
        else
            call getvid(keywordFact, 'EPXX', iocc=iocc, scal=kepxx, nbret=nxx)
            call getvid(keywordFact, 'EPYY', iocc=iocc, scal=kepyy, nbret=nyy)
            call getvid(keywordFact, 'EPZZ', iocc=iocc, scal=kepzz, nbret=nzz)
            call getvid(keywordFact, 'EPXY', iocc=iocc, scal=kepxy, nbret=nxy)
            call getvid(keywordFact, 'EPXZ', iocc=iocc, scal=kepxz, nbret=nxz)
            call getvid(keywordFact, 'EPYZ', iocc=iocc, scal=kepyz, nbret=nyz)
            call getvid(keywordFact, 'EPX', iocc=iocc, scal=kepx, nbret=nex)
            call getvid(keywordFact, 'KY', iocc=iocc, scal=kxky, nbret=nky)
            call getvid(keywordFact, 'KZ', iocc=iocc, scal=kxkz, nbret=nkz)
            call getvid(keywordFact, 'EXX', iocc=iocc, scal=kxexx, nbret=nexx)
            call getvid(keywordFact, 'EYY', iocc=iocc, scal=kxeyy, nbret=neyy)
            call getvid(keywordFact, 'EXY', iocc=iocc, scal=kxexy, nbret=nexy)
            call getvid(keywordFact, 'KXX', iocc=iocc, scal=kxkxx, nbret=nkxx)
            call getvid(keywordFact, 'KYY', iocc=iocc, scal=kxkyy, nbret=nkyy)
            call getvid(keywordFact, 'KXY', iocc=iocc, scal=kxkxy, nbret=nkxy)
            call getvid(keywordFact, 'KN1', iocc=iocc, scal=kxkn1, nbret=nkn1)
            call getvid(keywordFact, 'KN2', iocc=iocc, scal=kxkn2, nbret=nkn2)
            call getvid(keywordFact, 'VECT_N1', iocc=iocc, scal=kxvn1, nbret=nvectn)
            call getvid(keywordFact, 'VECT_N2', iocc=iocc, scal=kxvn2, nbret=nvectn)
            call getvid(keywordFact, 'VECT_N3', iocc=iocc, scal=kxvn3, nbret=nvectn)
            do i = 1, ncmp
                zk8(jvalv-1+i) = '&FOZERO'
            end do
            if (nxx .ne. 0) zk8(jvalv-1+1) = kepxx
            if (nyy .ne. 0) zk8(jvalv-1+2) = kepyy
            if (nzz .ne. 0) zk8(jvalv-1+3) = kepzz
            if (nxy .ne. 0) zk8(jvalv-1+4) = kepxy
            if (nxz .ne. 0) zk8(jvalv-1+5) = kepxz
            if (nyz .ne. 0) zk8(jvalv-1+6) = kepyz

            if (nex .ne. 0) zk8(jvalv-1+7) = kepx
            if (nky .ne. 0) zk8(jvalv-1+8) = kxky
            if (nkz .ne. 0) zk8(jvalv-1+9) = kxkz
            if (nkn1 .ne. 0) zk8(jvalv-1+8) = kxkn1
            if (nkn2 .ne. 0) zk8(jvalv-1+9) = kxkn2

            if (nexx .ne. 0) zk8(jvalv-1+10) = kxexx
            if (neyy .ne. 0) zk8(jvalv-1+11) = kxeyy
            if (nexy .ne. 0) zk8(jvalv-1+12) = kxexy
            if (nkxx .ne. 0) zk8(jvalv-1+13) = kxkxx
            if (nkyy .ne. 0) zk8(jvalv-1+14) = kxkyy
            if (nkxy .ne. 0) zk8(jvalv-1+15) = kxkxy

            if (nvectn .ne. 0) then
                zk8(jvalv-1+16) = kxvn1
                zk8(jvalv-1+17) = kxvn2
                zk8(jvalv-1+18) = kxvn3
            end if
        end if
!
        call getvtx(keywordFact, 'TOUT', iocc=iocc, scal=k8b, nbret=nbtou)
        if (nbtou .ne. 0) then
            call nocart(carte, 1, ncmp)
        else
            call reliem(model, mesh, 'NU_MAILLE', keywordFact, iocc, 2, motcle, typmcl, &
                        mesmai, nbma)
            if (nbma .eq. 0) goto 20
            call jeveuo(mesmai, 'L', jma)
            call nocart(carte, 3, ncmp, mode='NUM', nma=nbma, limanu=zi(jma))
            call jedetr(mesmai)
        end if
20      continue
    end do
!
    call jedema()
end subroutine
