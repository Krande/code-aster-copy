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

subroutine mdfedy(nbpal, nbmode, numpas, dt, dtsto,&
                  tcf, vrotat, dplmod, depgen, vitgen,&
                  fexgen, typal, finpal, cnpal, prdeff,&
                  conv, fsauv)
! aslint: disable=W1306
    implicit none
#include "asterf_types.h"
#include "asterfort/envdep.h"
#include "asterfort/recfor.h"
#include "asterfort/nlget.h"
    integer :: numpas, nbmode, nbpal
    real(kind=8) :: depgen(*), vitgen(*), fexgen(*)
    real(kind=8) :: dplmod(nbpal, nbmode, *), dt, dtsto, tcf
    real(kind=8) :: vrotat, conv
    aster_logical :: prdeff
!
! person_in_charge: nicolas.greffet at edf.fr
! ======================================================================
!
!              RECUPERATION DES FORCES VENANT D'EDYOS
!                 ET ENVOI DES CHAMPS CINEMATIQUES
!     ------------------------------------------------------------------
! IN  : NBMODE : NOMBRE DE MODES
! IN  : DEPGEN : DEPLACEMENTS GENERALISES
! IN  : VITGEN : VITESSES GENERALISEES
! VAR : FEXGEN : FORCES GENERALISEES
! IN  : NBPAL  : NOMBRE DE PALIERS
! IN  : DPLMOD : TABLEAU DES DEPL MODAUX AUX NOEUDS DE CHOC
!
! IN  : TCF    : INSTANT DE CALCUL
! IN  : DT     : PAS DE TEMPS
! IN  : DTSTO  : INSTANT DE STOCKAGE
! IN  : VROTAT : VITESSE DE ROTATION
! IN  : NUMDDL : NUMEROTATION DDL
!
! OUT  : CONV   : INDICATEUR DE CONVERGENCE EDYOS
! ----------------------------------------------------------------------
    integer :: i, j, k, l, palmax
    real(kind=8) :: dep(nbpal, 6), vit(nbpal, 6), force(nbpal, 3)
!-----------------------------------------------------------------------
    parameter (palmax=20)
    character(len=3) :: finpal(palmax)
    character(len=6) :: typal(palmax)
    character(len=8) :: cnpal(palmax), sd_nl
    real(kind=8) :: fsauv(palmax, 3)

    real(kind=8)          , pointer :: modal_depl_no1(:) => null()
!
    sd_nl = '&&OP29NL'
    do j = 1, nbpal
        do l = 1, 6
            dep(j,l)= 0.d0
            vit(j,l)= 0.d0
        end do
        do i = 1, nbmode
            do k = 1, 6
                dep(j,k)= dep(j,k)+ dplmod(j,i,k)*depgen(i)
                vit(j,k)= vit(j,k)+ dplmod(j,i,k)*vitgen(i)
            end do
        end do
    end do
    prdeff = .true.
!   ENVOI DES CHAMPS CINEMTATIQUES A EDYOS
    call envdep(numpas, nbpal, dt, dtsto, tcf,&
                dep, vit, vrotat, finpal, prdeff)
!   RECEPTION DES EFFORTS VENANT D'EDYOS
    call recfor(numpas, nbpal, force, typal, finpal,&
                cnpal, prdeff, conv)
!   COMBINAISON DES EFFORTS GENERALISES
    if (conv .gt. 0.0d0) then
        do j = 1, nbpal
            fsauv(j,1)=force(j,1)
            fsauv(j,2)=force(j,2)
            fsauv(j,2)=force(j,3)
            do i = 1, nbmode
                call nlget(sd_nl, _MODAL_DEPL_NO1, iocc=j, vr=modal_depl_no1)
                fexgen(i) = fexgen(i)+dplmod(j, i, 1)*force(j, 1)+&
                                      dplmod(j, i, 2)*force(j, 2)+&
                                      dplmod(j, i, 3)*force(j, 3)
            end do
        end do
    else
!      NON CONVERGENCE EDYOS : ON UTILISE FSAUV
        do j = 1, nbpal
            fsauv(j,1)=force(j,1)
            fsauv(j,2)=force(j,2)
            fsauv(j,2)=force(j,3)
            do i = 1, nbmode
                call nlget(sd_nl, _MODAL_DEPL_NO1, iocc=j, vr=modal_depl_no1)
                fexgen(i) = fexgen(i)+dplmod(j, i, 1)*fsauv(j, 1)+&
                                      dplmod(j, i, 2)*fsauv(j, 2)+&
                                      dplmod(j, i, 3)*fsauv(j, 3)
            end do
        end do
    endif
!
end subroutine
