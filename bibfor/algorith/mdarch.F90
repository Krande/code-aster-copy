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

subroutine mdarch(typcal, isto1, ipas, disc, nbmode, &
                  iorsto, discst, dt, depger, vitger, &
                  accger, depstr, vitstr, accstr, passto, &
                  nbsym, nomsym, depgec, vitgec, accgec, &
                  depstc, vitstc, accstc)

! aslint: disable=W1504
    implicit none
#include "asterfort/assert.h"
!
!   Obligatory arguments
    character(len=4), intent(in)  :: typcal
    integer(kind=8), intent(in)  :: isto1
    integer(kind=8), intent(in)  :: ipas
    real(kind=8), intent(in)  :: disc
    integer(kind=8), intent(in)  :: nbmode
    integer(kind=8)                        :: iorsto(*)
    real(kind=8)                   :: discst(*)
!
!   Optional arguments
!     typcal = 'TRAN' case
    real(kind=8), optional, intent(in)  :: dt
    real(kind=8), optional, intent(in)  :: depger(nbmode)
    real(kind=8), optional, intent(in)  :: vitger(nbmode)
    real(kind=8), optional, intent(in)  :: accger(nbmode)
    real(kind=8), optional              :: depstr(*)
    real(kind=8), optional              :: vitstr(*)
    real(kind=8), optional              :: accstr(*)
    real(kind=8), optional              :: passto(*)
!     typcal = 'HARM' case
    integer(kind=8), optional, intent(in)  :: nbsym
    character(len=4), optional, intent(in)  :: nomsym(*)
    complex(kind=8), optional, intent(in)  :: depgec(nbmode)
    complex(kind=8), optional, intent(in)  :: vitgec(nbmode)
    complex(kind=8), optional, intent(in)  :: accgec(nbmode)
    complex(kind=8), optional              :: depstc(*)
    complex(kind=8), optional              :: vitstc(*)
    complex(kind=8), optional              :: accstc(*)
!
!-----------------------------------------------------------------------
!
!   ARCHIVAGE DES CHAMPS GENERALISES OBLIGATOIRES POUR LA SD_DYNA_GENE
!
!-----------------------------------------------------------------------
    integer(kind=8) :: im, ind, ich, inom, nbsym2
    character(len=4) :: nomsym2(3)
!-----------------------------------------------------------------------
!
    nbsym2 = 3
    nomsym2 = ['DEPL', 'VITE', 'ACCE']
    ASSERT(ENSEMBLE2(nomsym, nbsym))
    if (present(nbsym)) then
        nbsym2 = nbsym
        do inom = 1, nbsym2
            nomsym2(inom) = nomsym(inom)
        end do
    end if

    ASSERT((typcal(1:4) .eq. 'TRAN') .or. (typcal(1:4) .eq. 'HARM'))
    if (typcal .eq. 'TRAN') then
        ASSERT(present(dt) .and. present(depger) .and. present(vitger))
        ASSERT(present(accger) .and. present(depstr) .and. present(vitstr))
        ASSERT(present(accstr) .and. present(passto))
        ASSERT(absent(nbsym) .and. absent(nomsym))
    else
        ASSERT(present(depgec) .and. present(vitgec) .and. present(accgec))
        ASSERT(present(depstc) .and. present(vitstc) .and. present(accstc))
        ASSERT(absent(dt) .and. absent(passto))
    end if
!
    iorsto(isto1+1) = ipas
    discst(isto1+1) = disc
    ind = nbmode*isto1
!
    if (typcal(1:4) .eq. 'TRAN') then
!
        passto(isto1+1) = dt
        do im = 1, nbmode
            depstr(ind+im) = depger(im)
            vitstr(ind+im) = vitger(im)
            accstr(ind+im) = accger(im)
        end do
    else
        do ich = 1, nbsym
            if (nomsym(ich) (1:4) .eq. 'DEPL') then
                do im = 1, nbmode
                    depstc(ind+im) = depgec(im)
                end do
            else if (nomsym(ich) (1:4) .eq. 'VITE') then
                do im = 1, nbmode
                    vitstc(ind+im) = vitgec(im)
                end do
            else
                do im = 1, nbmode
                    accstc(ind+im) = accgec(im)
                end do
            end if
        end do
!
    end if
!
end subroutine
