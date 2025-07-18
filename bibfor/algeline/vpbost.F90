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
subroutine vpbost(typres, nbmode, nbvect, omeshi, valpro, &
                  nvpro, vpinf, vpmax, precdc, method, &
                  omecor)
!     RECTIFIE LES VALEURS PROPRES
!-----------------------------------------------------------------------
!     IN  : TYPRES  : TYPE DE RESULTAT (DYNAMIQUE OU FLAMBEMENT)
!     IN  : NBMODE  : NOMBRE DE MODE DEMANDES
!     IN  : NBVECT  : NOMBRE DE VECTEURS UTILISES AU COURS DU CALCUL
!     IN  : OMESHI  : DECALAGE UTILISE POUR LE CALCUL
!     IN  : OMECOR  : OMEGA2 DE CORPS RIGIDE
!     IN  : VALPRO  : VALEURS PROPRES
!     IN  : NVPRO   : DIMENSION DU VECTEUR VALPRO
!     IN  : METHOD  : TYPE DE METHODE
!     IN  : PRECDC  : POURCENTAGE D'AUGMENTATION DES BORNES
!     OUT : VPINF : PLUS PETIT OMEGA2 CALCULE ET RETENU
!     OUT : VPMAX : PLUS GRAND OMEGA2 CALCULE ET RETENU
!----------------------------------------------------------------------
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/freqom.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nbmode, nbvect, nvpro
    real(kind=8) :: valpro(nvpro), precdc, omeshi, omecor, vpinf, vpmax
    character(len=8) :: method
    character(len=16) :: typres
!     ------------------------------------------------------------------
    real(kind=8) :: vpinf2, vpmax2, tole
    real(kind=8) :: valr(2)
    aster_logical :: loginf, logmax
    integer(kind=8) :: niv, ifm, i
!     ------------------------------------------------------------------
!
!     ------------------------------------------------------------------
!     --------  RECTIFICATION DES FREQUENCES DUE AU SHIFT  -------------
!     --------     DETERMINATION DE LA POSITION MODALE     -------------
!     ------------------------------------------------------------------
!
!     ---RECUPERATION DU NIVEAU D'IMPRESSION----
    call infniv(ifm, niv)
!     -------------------------------------------
!
    do i = 1, nbvect
        valpro(i) = valpro(i)+omeshi
    end do
!
    vpinf = valpro(1)
    vpmax = valpro(1)
    do i = 2, nbmode
        if (valpro(i) .lt. vpinf) then
            vpinf = valpro(i)
        end if
        if (valpro(i) .gt. vpmax) then
            vpmax = valpro(i)
        end if
    end do
    if (niv .ge. 1) then
        write (ifm, 1600)
        if (typres .eq. 'DYNAMIQUE') then
            valr(1) = freqom(vpinf)
            valr(2) = freqom(vpmax)
            call utmess('I', 'ALGELINE6_16', nr=2, valr=valr)
        else
            valr(1) = -vpmax
            valr(2) = -vpinf
            call utmess('I', 'ALGELINE6_17', nr=2, valr=valr)
        end if
    end if
!
    if (method .eq. 'SORENSEN') then
        if (abs(vpmax) .le. omecor) then
            vpmax = omecor
        end if
        if (abs(vpinf) .le. omecor) then
            vpinf = -omecor
        end if
        vpinf = vpinf*(1.d0-sign(precdc, vpinf))
        vpmax = vpmax*(1.d0+sign(precdc, vpmax))
    end if
!
!     -----POUR LES OPTIONS JACOBI ET LANCZOS---
!
    loginf = .false.
    logmax = .false.
    if (method .ne. 'SORENSEN') then
        do i = nbmode+1, nbvect
            if (valpro(i) .le. vpinf) then
                if (.not. loginf) then
                    loginf = .true.
                    vpinf2 = valpro(i)
                end if
            end if
            if (valpro(i) .ge. vpmax) then
                if (.not. logmax) then
                    logmax = .true.
                    vpmax2 = valpro(i)
                end if
            end if
        end do
!
!     ----ON REGARDE L'ECART QU'IL Y A ENTRE FREQMIN ET LA
!         FREQUENCE PRECEDENTE, PUIS ON RECALCULE FREQMIN-----
!
        if (loginf) then
            if (vpinf2 .lt. vpinf) then
                if (vpmax .gt. r8prem()) then
                    tole = (abs(vpinf2-vpinf)/vpinf)
                    if (tole .lt. precdc) then
                        call utmess('A', 'ALGELINE3_58')
                        valr(1) = freqom(vpinf2)
                        call utmess('A', 'ALGELINE4_66', sr=valr(1))
                        vpinf = vpinf*(1.d0-sign(precdc, vpinf))
                    end if
                else
                    tole = abs(vpinf2-vpinf)
                    if (tole .lt. precdc) then
                        call utmess('A', 'ALGELINE3_58')
                        valr(1) = freqom(vpinf2)
                        call utmess('A', 'ALGELINE4_66', sr=valr(1))
                        vpinf = vpinf*(1.d0-sign(precdc, vpinf))
                    end if
                end if
                vpinf = 0.5d0*(vpinf+vpinf2)
            else
                vpinf = vpinf*(1.d0-sign(precdc, vpinf))
            end if
        else
            vpinf = vpinf*(1.d0-sign(precdc, vpinf))
        end if
!
!     -----ON FAIT LES MEMES CALCULS AVEC FREQMAX------
!
        if (logmax) then
            if (vpmax2 .gt. vpmax) then
                if (vpinf .gt. r8prem()) then
                    tole = (abs(vpmax2-vpmax)/vpmax)
                    if (tole .lt. precdc) then
                        call utmess('A', 'ALGELINE3_58')
                        valr(1) = freqom(vpmax2)
                        call utmess('A', 'ALGELINE4_66', sr=valr(1))
                        vpmax = vpmax*(1.d0+sign(precdc, vpmax))
                    end if
                else
                    tole = abs(vpmax2-vpmax)
                    if (tole .lt. precdc) then
                        call utmess('A', 'ALGELINE3_58')
                        valr(1) = freqom(vpmax2)
                        call utmess('A', 'ALGELINE4_66', sr=valr(1))
                        vpmax = vpmax*(1.d0+sign(precdc, vpmax))
                    end if
                end if
                vpmax = 0.5d0*(vpmax+vpmax2)
            else
                vpmax = vpmax*(1.d0+sign(precdc, vpmax))
            end if
        else
            vpmax = vpmax*(1.d0+sign(precdc, vpmax))
        end if
    end if
!
!     -----DETERMINATION DE FREQMIN ET FREQMAX-----
!
    if (abs(vpmax) .le. omecor) then
        vpmax = omecor
    end if
    if (abs(vpinf) .le. omecor) then
        vpinf = -omecor
    end if
!
!      -----IMPRESSIONS-----
!
    if (loginf) then
        if (niv .ge. 1) then
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_18', sr=freqom(vpinf2))
            else
                call utmess('I', 'ALGELINE6_20', sr=vpinf2)
            end if
        end if
!
    end if
    if (logmax) then
        if (niv .ge. 1) then
            if (typres .eq. 'DYNAMIQUE') then
                call utmess('I', 'ALGELINE6_19', sr=freqom(vpmax2))
            else
                call utmess('I', 'ALGELINE6_21', sr=vpmax2)
            end if
        end if
    end if
!
    if (niv .ge. 1) then
        write (ifm, 1600)
    end if
!
1600 format(72('-'))
!
end subroutine
