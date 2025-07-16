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
subroutine cachre(load, model, mesh, geomDime, valeType, &
                  param, keywordFactZ)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterfort/alcart.h"
#include "asterfort/assert.h"
#include "asterfort/char_affe_neum.h"
#include "asterfort/getvc8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
#include "asterfort/tecart.h"
!

    character(len=8), intent(in) :: load, mesh, model
    integer(kind=8), intent(in) :: geomDime
    character(len=4), intent(in) :: valeType
    character(len=5), intent(in) :: param
    character(len=*), intent(in) :: keywordFactZ
!
! --------------------------------------------------------------------------------------------------
!
! Loads affectation
!
! Treatment of loads 'FORCE_CONTOUR' , 'FORCE_INTERNE' , 'FORCE_ARETE'
!                    'FORCE_FACE'    , 'FORCE_POUTRE'  , 'FORCE_COQUE'
!
! --------------------------------------------------------------------------------------------------
!
! In  load             : load
! In  model            : model
! In  mesh             : mesh
! In  geomDime         : dimension of space
! In  valeType         : affected value type (real, complex or function)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, n, nchre, nrep, ncmp, jvalv, iocc, nfx, nfy, nfz
    integer(kind=8) :: nmx, nmy, nmz, nplan, nmgx, nmgy, nmgz
    real(kind=8) :: fx, fy, fz, mx, my, mz, vpre, mgx, mgy, mgz
    complex(kind=8) :: cfx, cfy, cfz, cmx, cmy, cmz, cvpre
    character(len=8) :: kfx, kfy, kfz, kmx, kmy, kmz, typch, plan
    character(len=8) :: kmgx, kmgy, kmgz
    character(len=16) :: keywordFact
    character(len=19) :: carte
    character(len=19) :: cartes(1)
    integer(kind=8) :: ncmps(1)
    character(len=8), pointer :: vncmp(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    keywordFact = keywordFactZ
    call getfac(keywordFact, nchre)
!
    carte = load(1:8)//'.CHME.'//param(1:5)
!
    if (valeType .eq. 'REEL') then
        call alcart('G', carte, mesh, 'FORC_R')
    else if (valeType .eq. 'FONC') then
        call alcart('G', carte, mesh, 'FORC_F')
    else if (valeType .eq. 'COMP') then
        call alcart('G', carte, mesh, 'FORC_C')
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(carte//'.NCMP', 'E', vk8=vncmp)
    call jeveuo(carte//'.VALV', 'E', jvalv)
!
! --- STOCKAGE DE FORCES NULLES SUR TOUT LE MAILLAGE
!     ET REPERE = 0.(SI 'REEL'),REPERE = 'GLOBAL' (SI FONC) ---
!
    vncmp(1) = 'FX'
    vncmp(2) = 'FY'
    vncmp(3) = 'FZ'
    vncmp(4) = 'MX'
    vncmp(5) = 'MY'
    vncmp(6) = 'MZ'
    vncmp(7) = 'REP'
    vncmp(8) = 'PLAN'
    vncmp(9) = 'MGX'
    vncmp(10) = 'MGY'
    vncmp(11) = 'MGZ'
!
    if (valeType(1:4) .eq. 'REEL') then
        do i = 1, 11
            zr(jvalv-1+i) = 0.d0
        end do
    else if (valeType(1:4) .eq. 'COMP') then
        do i = 1, 11
            zc(jvalv-1+i) = dcmplx(0.d0, 0.d0)
        end do
    else if (valeType .eq. 'FONC') then
        do i = 1, 6
            zk8(jvalv-1+i) = '&FOZERO'
        end do
        zk8(jvalv-1+7) = 'GLOBAL'
        zk8(jvalv-1+8) = '&FOZERO'
        zk8(jvalv-1+9) = '&FOZERO'
        zk8(jvalv-1+10) = '&FOZERO'
        zk8(jvalv-1+11) = '&FOZERO'
    else
        ASSERT(.false.)
    end if
    call nocart(carte, 1, 11)
!
! --- STOCKAGE DANS LA CARTE ---
!
    do iocc = 1, nchre
        nrep = 0
        ncmp = 0
        if (keywordFact .eq. 'FORCE_POUTRE') then
            call getvtx(keywordFact, 'TYPE_CHARGE', iocc=iocc, scal=typch, nbret=n)
            if (typch .eq. 'VENT') nrep = 2
        end if
        if (valeType .eq. 'COMP') then
            call getvc8(keywordFact, 'FX', iocc=iocc, scal=cfx, nbret=nfx)
            call getvc8(keywordFact, 'FY', iocc=iocc, scal=cfy, nbret=nfy)
            call getvc8(keywordFact, 'FZ', iocc=iocc, scal=cfz, nbret=nfz)
            if (keywordFact .ne. 'FORCE_INTERNE' .and. &
                keywordFact .ne. 'FORCE_POUTRE' .and. &
                keywordFact .ne. 'FORCE_FACE') then
                call getvc8(keywordFact, 'MX', iocc=iocc, scal=cmx, nbret=nmx)
                call getvc8(keywordFact, 'MY', iocc=iocc, scal=cmy, nbret=nmy)
                call getvc8(keywordFact, 'MZ', iocc=iocc, scal=cmz, nbret=nmz)
            else
                nmx = 0
                nmy = 0
                nmz = 0
            end if
            if (nfx+nfy+nfz+nmx+nmy+nmz .eq. 0) then
                if (keywordFact .eq. 'FORCE_POUTRE') then
                    nrep = 1
                    call getvc8(keywordFact, 'N', iocc=iocc, scal=cfx, nbret=nfx)
                    call getvc8(keywordFact, 'VY', iocc=iocc, scal=cfy, nbret=nfy)
                    call getvc8(keywordFact, 'VZ', iocc=iocc, scal=cfz, nbret=nfz)
                else if (keywordFact .eq. 'FORCE_COQUE') then
                    call getvc8(keywordFact, 'PRES', iocc=iocc, scal=cvpre, nbret=nfz)
                    if (nfz .eq. 0) then
                        nrep = 1
                        call getvc8(keywordFact, 'F1', iocc=iocc, scal=cfx, nbret=nfx)
                        call getvc8(keywordFact, 'F2', iocc=iocc, scal=cfy, nbret=nfy)
                        call getvc8(keywordFact, 'F3', iocc=iocc, scal=cfz, nbret=nfz)
                        call getvc8(keywordFact, 'MF1', iocc=iocc, scal=cmx, nbret=nmx)
                        call getvc8(keywordFact, 'MF2', iocc=iocc, scal=cmy, nbret=nmy)
                        nmz = 0
                    else
                        nrep = 3
                        cfz = cvpre
                        nfx = 0
                        nfy = 0
                        nmx = 0
                        nmy = 0
                        nmz = 0
                    end if
                end if
            end if
            if (nfx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FX'
                zc(jvalv-1+ncmp) = cfx
            end if
            if (nfy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FY'
                zc(jvalv-1+ncmp) = cfy
            end if
            if (nfz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FZ'
                zc(jvalv-1+ncmp) = cfz
            end if
            if (nmx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MX'
                zc(jvalv-1+ncmp) = cmx
            end if
            if (nmy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MY'
                zc(jvalv-1+ncmp) = cmy
            end if
            if (nmz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MZ'
                zc(jvalv-1+ncmp) = cmz
            end if
            if (nrep .ge. 1) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'REP'
                if (nrep .eq. 1) zc(jvalv-1+ncmp) = dcmplx(1.d0, 1.d0)
                if (nrep .eq. 2) zc(jvalv-1+ncmp) = dcmplx(2.d0, 2.d0)
                if (nrep .eq. 1) zc(jvalv-1+ncmp) = 1.d0
                if (nrep .eq. 2) zc(jvalv-1+ncmp) = 2.d0
            end if
        else if (valeType .eq. 'REEL') then
            call getvr8(keywordFact, 'FX', iocc=iocc, scal=fx, nbret=nfx)
            call getvr8(keywordFact, 'FY', iocc=iocc, scal=fy, nbret=nfy)
            call getvr8(keywordFact, 'FZ', iocc=iocc, scal=fz, nbret=nfz)
            if (keywordFact .ne. 'FORCE_INTERNE' .and. keywordFact .ne. 'FORCE_FACE' &
                .and. keywordFact .ne. 'FORCE_CONTOUR') then
                call getvr8(keywordFact, 'MX', iocc=iocc, scal=mx, nbret=nmx)
                call getvr8(keywordFact, 'MY', iocc=iocc, scal=my, nbret=nmy)
                call getvr8(keywordFact, 'MZ', iocc=iocc, scal=mz, nbret=nmz)
            else
                nmx = 0
                nmy = 0
                nmz = 0
            end if
            if (keywordFact .eq. 'FORCE_POUTRE') then
                call getvr8(keywordFact, 'MGX', iocc=iocc, scal=mgx, nbret=nmgx)
                call getvr8(keywordFact, 'MGY', iocc=iocc, scal=mgy, nbret=nmgy)
                call getvr8(keywordFact, 'MGZ', iocc=iocc, scal=mgz, nbret=nmgz)
            else
                nmgx = 0
                nmgy = 0
                nmgz = 0
            end if
            if (nfx+nfy+nfz+nmx+nmy+nmz .eq. 0) then
                if (keywordFact .eq. 'FORCE_POUTRE') then
                    nrep = 1
                    call getvr8(keywordFact, 'N', iocc=iocc, scal=fx, nbret=nfx)
                    call getvr8(keywordFact, 'VY', iocc=iocc, scal=fy, nbret=nfy)
                    call getvr8(keywordFact, 'VZ', iocc=iocc, scal=fz, nbret=nfz)
                    call getvr8(keywordFact, 'MT', iocc=iocc, scal=mx, nbret=nmx)
                    call getvr8(keywordFact, 'MFY', iocc=iocc, scal=my, nbret=nmy)
                    call getvr8(keywordFact, 'MFZ', iocc=iocc, scal=mz, nbret=nmz)
                else if (keywordFact .eq. 'FORCE_COQUE') then
                    nrep = 1
                    call getvr8(keywordFact, 'PRES', iocc=iocc, scal=vpre, nbret=nfz)
                    if (nfz .eq. 0) then
                        call getvr8(keywordFact, 'F1', iocc=iocc, scal=fx, nbret=nfx)
                        call getvr8(keywordFact, 'F2', iocc=iocc, scal=fy, nbret=nfy)
                        call getvr8(keywordFact, 'F3', iocc=iocc, scal=fz, nbret=nfz)
                        call getvr8(keywordFact, 'MF1', iocc=iocc, scal=mx, nbret=nmx)
                        call getvr8(keywordFact, 'MF2', iocc=iocc, scal=my, nbret=nmy)
                        nmz = 0
                    else
                        fz = vpre
                        nfx = 0
                        nfy = 0
                        nmx = 0
                        nmy = 0
                        nmz = 0
                        nrep = 3
                    end if
                end if
            end if
            if (nfx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FX'
                zr(jvalv-1+ncmp) = fx
            end if
            if (nfy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FY'
                zr(jvalv-1+ncmp) = fy
            end if
            if (nfz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FZ'
                zr(jvalv-1+ncmp) = fz
            end if
            if (nmx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MX'
                zr(jvalv-1+ncmp) = mx
            end if
            if (nmy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MY'
                zr(jvalv-1+ncmp) = my
            end if
            if (nmz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MZ'
                zr(jvalv-1+ncmp) = mz
            end if
            if (nmgx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MGX'
                zr(jvalv-1+ncmp) = mgx
            end if
            if (nmgy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MGY'
                zr(jvalv-1+ncmp) = mgy
            end if
            if (nmgz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MGZ'
                zr(jvalv-1+ncmp) = mgz
            end if

            if (nrep .ge. 1) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'REP'
                if (nrep .eq. 1) zr(jvalv-1+ncmp) = 1.d0
                if (nrep .eq. 2) zr(jvalv-1+ncmp) = 2.d0
                if (nrep .eq. 3) zr(jvalv-1+ncmp) = 3.d0
! --           (NREP=3) CAS D UNE PRESSION --> ON PREND L OPPOSE DE
! --           LA VALEUR LUE DANS LE TE
            end if
        else
            call getvid(keywordFact, 'FX', iocc=iocc, scal=kfx, nbret=nfx)
            call getvid(keywordFact, 'FY', iocc=iocc, scal=kfy, nbret=nfy)
            call getvid(keywordFact, 'FZ', iocc=iocc, scal=kfz, nbret=nfz)
            if (keywordFact .ne. 'FORCE_INTERNE' .and. keywordFact .ne. 'FORCE_FACE' &
                .and. keywordFact .ne. 'FORCE_CONTOUR') then
                call getvid(keywordFact, 'MX', iocc=iocc, scal=kmx, nbret=nmx)
                call getvid(keywordFact, 'MY', iocc=iocc, scal=kmy, nbret=nmy)
                call getvid(keywordFact, 'MZ', iocc=iocc, scal=kmz, nbret=nmz)
            else
                nmx = 0
                nmy = 0
                nmz = 0
            end if
            if (keywordFact .eq. 'FORCE_POUTRE') then
                call getvid(keywordFact, 'MGX', iocc=iocc, scal=kmgx, nbret=nmgx)
                call getvid(keywordFact, 'MGY', iocc=iocc, scal=kmgy, nbret=nmgy)
                call getvid(keywordFact, 'MGZ', iocc=iocc, scal=kmgz, nbret=nmgz)
            else
                nmgx = 0
                nmgy = 0
                nmgz = 0
            end if
            if (nfx+nfy+nfz+nmx+nmy+nmz .eq. 0) then
                if (keywordFact .eq. 'FORCE_POUTRE') then
                    nrep = 1
                    call getvid(keywordFact, 'N', iocc=iocc, scal=kfx, nbret=nfx)
                    call getvid(keywordFact, 'VY', iocc=iocc, scal=kfy, nbret=nfy)
                    call getvid(keywordFact, 'VZ', iocc=iocc, scal=kfz, nbret=nfz)
                    call getvid(keywordFact, 'MT', iocc=iocc, scal=kmx, nbret=nmx)
                    call getvid(keywordFact, 'MFY', iocc=iocc, scal=kmy, nbret=nmy)
                    call getvid(keywordFact, 'MFZ', iocc=iocc, scal=kmz, nbret=nmz)
                else if (keywordFact(1:11) .eq. 'FORCE_COQUE') then
                    nrep = 1
                    call getvid(keywordFact, 'PRES', iocc=iocc, scal=kfz, nbret=nfz)
                    if (nfz .eq. 0) then
                        call getvid(keywordFact, 'F1', iocc=iocc, scal=kfx, nbret=nfx)
                        call getvid(keywordFact, 'F2', iocc=iocc, scal=kfy, nbret=nfy)
                        call getvid(keywordFact, 'F3', iocc=iocc, scal=kfz, nbret=nfz)
                        call getvid(keywordFact, 'MF1', iocc=iocc, scal=kmx, nbret=nmx)
                        call getvid(keywordFact, 'MF2', iocc=iocc, scal=kmy, nbret=nmy)
                        nmz = 0
                    else
                        nfx = 0
                        nfy = 0
                        nmx = 0
                        nmy = 0
                        nmz = 0
                        nrep = 3
                    end if
                end if
            end if
            if (nfx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FX'
                zk8(jvalv-1+ncmp) = kfx
            end if
            if (nfy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FY'
                zk8(jvalv-1+ncmp) = kfy
            end if
            if (nfz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'FZ'
                zk8(jvalv-1+ncmp) = kfz
            end if
            if (nmx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MX'
                zk8(jvalv-1+ncmp) = kmx
            end if
            if (nmy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MY'
                zk8(jvalv-1+ncmp) = kmy
            end if
            if (nmz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MZ'
                zk8(jvalv-1+ncmp) = kmz
            end if
            if (nmgx .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MGX'
                zk8(jvalv-1+ncmp) = kmgx
            end if
            if (nmgy .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MGY'
                zk8(jvalv-1+ncmp) = kmgy
            end if
            if (nmgz .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'MGZ'
                zk8(jvalv-1+ncmp) = kmgz
            end if

            if (nrep .ge. 1) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'REP'
                if (nrep .eq. 1) zk8(jvalv-1+ncmp) = 'LOCAL'
                if (nrep .eq. 2) zk8(jvalv-1+ncmp) = 'VENT'
                if (nrep .eq. 3) zk8(jvalv-1+ncmp) = 'LOCAL_PR'
! --           (NREP=3) CAS D UNE PRESSION --> ON PREND L OPPOSE DE
! --           LA VALEUR LUE DANS LE TE
            end if
        end if
        if (ncmp .eq. 0) goto 20
!
        if (keywordFact .eq. 'FORCE_COQUE') then
            call getvtx(keywordFact, 'PLAN', iocc=iocc, scal=plan, nbret=nplan)
            if (nplan .ne. 0) then
                ncmp = ncmp+1
                vncmp(ncmp) = 'PLAN'
                if (valeType .eq. 'REEL') then
                    if (plan .eq. 'MAIL') then
                        zr(jvalv-1+ncmp) = dble(0)
                    else if (plan .eq. 'INF') then
                        zr(jvalv-1+ncmp) = dble(-1)
                    else if (plan .eq. 'SUP') then
                        zr(jvalv-1+ncmp) = dble(1)
                    else if (plan .eq. 'MOY') then
                        zr(jvalv-1+ncmp) = dble(2)
                    end if
                else if (valeType .eq. 'FONC') then
                    zk8(jvalv-1+ncmp) = plan
                end if
            end if
        end if
!
        cartes(1) = carte
        ncmps(1) = ncmp
        call char_affe_neum(model, mesh, geomDime, keywordFact, iocc, 1, &
                            cartes, ncmps)
!
20      continue
    end do
!
    call tecart(carte)
    call jedema()
end subroutine
