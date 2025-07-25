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

subroutine ef0154(nomte)
!
! --------------------------------------------------------------------------------------------------
!
!     Calcul  EFGE_ELNO pour MECA_BARRE MECA_2D_BARRE
!
! --------------------------------------------------------------------------------------------------
!
    implicit none
    character(len=16) :: nomte
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/get_value_mode_local.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/verift.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codres(1)
    integer(kind=8) :: i, jdepl, jeffo
    integer(kind=8) :: lmater, lorien, iret, nc, nno

    real(kind=8) :: pgl(3, 3), klc(6, 6)
    real(kind=8) :: ugr(6), ulr(6), flr(6)
    real(kind=8) :: a, epsth, e, r8bid, xfl1, xfl4, xl, xrig, val(1)

    character(len=4)    :: fami
    character(len=16)   :: ch16
!
    aster_logical   :: lteimp
!
    real(kind=8)        :: valr(2)
    character(len=8)    :: valp(2)
!
! --------------------------------------------------------------------------------------------------
!
    lteimp = ASTER_FALSE
    nno = 2
    nc = 3
    fami = 'RIGI'
!
    if ((nomte .ne. 'MECA_BARRE') .and. &
        (nomte .ne. 'MECA_2D_BARRE') .and. &
        (nomte .ne. 'MECABL2')) then
        ch16 = nomte
        call utmess('F', 'ELEMENTS2_42', sk=ch16)
    end if
!
!   RECUPERATION DES CARACTERISTIQUES MATERIAUX ---
    call jevech('PMATERC', 'L', lmater)
!
    call verift(fami, 1, 1, '+', zi(lmater), epsth_=epsth)
!
    r8bid = 0.0d0
    call rcvalb(fami, 1, 1, '+', zi(lmater), ' ', 'ELAS', 0, ' ', [r8bid], &
                1, 'E', val, codres, 1)
    e = val(1)
    if (epsth .ne. 0.d0) lteimp = ASTER_TRUE
!
!   Longueur de l'élément
!   Caracteristiques de la section
    if (nomte .eq. 'MECA_BARRE') then
        xl = lonele()
        valp(1) = 'A1'
        call get_value_mode_local('PCAGNBA', valp, valr, iret, nbpara_=1)
    else if (nomte .eq. 'MECA_2D_BARRE') then
        xl = lonele(dime=2)
        valp(1) = 'A1'
        call get_value_mode_local('PCAGNBA', valp, valr, iret, nbpara_=1)
    else if (nomte .eq. 'MECABL2') then
        xl = lonele()
        valp(1) = 'A1'
        call get_value_mode_local('PCACABL', valp, valr, iret, nbpara_=1)
    else
        xl = 0.0d0
        ASSERT(ASTER_FALSE)
    end if
    a = valr(1)
!
!   RECUPERATION DES ORIENTATIONS ALPHA,BETA,GAMMA ---
    call jevech('PCAORIE', 'L', lorien)
!   MATRICE DE ROTATION PGL
    call matrot(zr(lorien), pgl)
!
!   RECUPERATION DES DEPLACEMENTS OU DES VITESSES
    ugr(:) = 0.d0
!
!   ON RECUPERE DES DEPLACEMENTS
    call jevech('PDEPLAR', 'L', jdepl)
    if ((nomte .eq. 'MECA_BARRE') .or. (nomte .eq. 'MECABL2')) then
        do i = 1, 6
            ugr(i) = zr(jdepl+i-1)
        end do
    else if (nomte .eq. 'MECA_2D_BARRE') then
        ugr(1) = zr(jdepl+1-1)
        ugr(2) = zr(jdepl+2-1)
        ugr(4) = zr(jdepl+3-1)
        ugr(5) = zr(jdepl+4-1)
    end if
!
!   VECTEUR DANS REPERE LOCAL  ULR = PGL * UGR
    call utpvgl(nno, nc, pgl, ugr, ulr)
!
!   RIGIDITE ELEMENTAIRE
    klc(:, :) = 0.0d0
!
    xrig = e*a/xl
    klc(1, 1) = xrig
    klc(1, 4) = -xrig
    klc(4, 1) = -xrig
    klc(4, 4) = xrig
!
!   VECTEUR EFFORT LOCAL  FLR = KLC * ULR
    call pmavec('ZERO', 6, klc, ulr, flr)
!
!   TENIR COMPTE DES EFFORTS DUS A LA DILATATION ---
    if (lteimp) then
!       CALCUL DES FORCES INDUITES ---
        xfl1 = -epsth*e*a
        xfl4 = -xfl1
        flr(1) = flr(1)-xfl1
        flr(4) = flr(4)-xfl4
    end if
!
    call jevech('PEFFORR', 'E', jeffo)
    zr(jeffo) = -flr(1)
    zr(jeffo+1) = flr(4)
!
end subroutine
