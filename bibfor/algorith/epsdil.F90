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
subroutine epsdil(npi, ipoids, ipoid2, ivf, ivf2,&
                  idfde, idfde2, geom, dimdef, dimuel,&
                  ndim, nddls, nddlm, nno, nnos,&
                  nnom, interp, axi, regula, deplp,&
                  defgep)
! ======================================================================
! person_in_charge: romeo.fernandes at edf.fr
! aslint: disable=W1306,W1504
    implicit none
#include "asterf_types.h"
#include "asterfort/cabrp0.h"
#include "asterfort/cabrp1.h"
#include "asterfort/cabrsl.h"
    aster_logical :: axi
    integer :: npi, ipoids, ipoid2, ivf, ivf2, idfde, idfde2, dimdef, dimuel
    integer :: ndim, nddls, nddlm, nno, nnos, nnom, regula(6)
    real(kind=8) :: geom(ndim, *), deplp(dimuel), defgep(npi*dimdef)
    character(len=2) :: interp
! ======================================================================
! --- BUT : CALCUL DE EPSI_ELGA -----------------------------------
! ======================================================================
    integer :: kpi, i, n
    real(kind=8) :: poids, poids2, b(dimdef, dimuel)
! ======================================================================
! --- BOUCLE SUR LES POINTS D'INTEGRATION ------------------------------
! ======================================================================
    do kpi = 1, npi
! ======================================================================
! --- DEFINITION DE L'OPERATEUR B (DEFINI PAR E=B.U) -------------------
! ======================================================================
        if (interp .eq. 'P0') then
            call cabrp0(kpi, ipoids, ipoid2, ivf, ivf2,&
                        idfde, idfde2, geom, dimdef, dimuel,&
                        ndim, nddls, nddlm, nno, nnos,&
                        nnom, axi, regula, b, poids,&
                        poids2)
        else if (interp.eq.'SL') then
            call cabrsl(kpi, ipoids, ipoid2, ivf, ivf2,&
                        idfde, idfde2, geom, dimdef, dimuel,&
                        ndim, nddls, nddlm, nno, nnos,&
                        nnom, axi, regula, b, poids,&
                        poids2)
        else if (interp.eq.'P1') then
            call cabrp1(kpi, ipoids, ipoid2, ivf, ivf2,&
                        idfde, idfde2, geom, dimdef, dimuel,&
                        ndim, nddls, nddlm, nno, nnos,&
                        nnom, axi, regula, b, poids,&
                        poids2)
        endif
! ======================================================================
! --- CALCUL DES DEFORMATIONS GENERALISEES E=B.U -----------------------
! ======================================================================
        do i = 1, dimdef
            defgep((kpi-1)*dimdef+i)=0.0d0
            do n = 1, dimuel
                defgep((kpi-1)*dimdef+i) = defgep( (kpi-1)*dimdef+i)+ b(i,n)*deplp(n)
            end do
        end do
    end do
! ======================================================================
end subroutine
