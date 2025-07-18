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

subroutine te0154(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
#include "asterfort/ptenci.h"
#include "asterfort/ptenpo.h"
#include "asterfort/ptenth.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "asterfort/verift.h"
    character(len=*) :: option, nomte
! ----------------------------------------------------------------------
!     CALCUL
!       - DU VECTEUR ELEMENTAIRE EFFORT GENERALISE,
!       - DU VECTEUR ELEMENTAIRE CONTRAINTE
!       - DE L'ENERGIE DE DEFORMATION
!       - DE L'ENERGIE CINETIQUE
!     POUR LES ELEMENTS DE BARRE
! ----------------------------------------------------------------------
! IN  OPTION : K16 : NOM DE L'OPTION A CALCULER
!        'SIEF_ELGA'   : CALCUL DU VECTEUR EFFORT GENERALISE
!        'EPSI_ELGA'   : CALCUL DU VECTEUR DEFORMATION
!        'EPOT_ELEM'   : CALCUL DE L'ENERGIE DE DEFORMATION
!        'ECIN_ELEM'   : CALCUL DE L'ENERGIE CINETIQUE
! IN  NOMTE  : K16 : NOM DU TYPE ELEMENT
!        'MECA_BARRE'   : BARRE
!        'MECA_2D_BARRE'   : BARRE
!
!
    real(kind=8) :: pgl(3, 3), klc(6, 6), enerth
    real(kind=8) :: ugr(6), ulr(12), flr(6)
    integer(kind=8) :: codres(1)
    character(len=3) :: stopz
    character(len=4) :: fami
    character(len=16) :: ch16
    aster_logical :: lteimp
    real(kind=8) :: aire, epsth, e(1), r8bid, rho(1), xfl1, xfl4, xl, xmas, xrig
    integer(kind=8) :: ii, iif, itype, j, jdepl, jeffo, jende, jfreq, jdefo, kanl
    integer(kind=8) :: lmater, lorien, lsect, iret, nc, nno
    integer(kind=8) :: jvite
!     ------------------------------------------------------------------
!
    r8bid = 0.d0
    lteimp = .false.
    nno = 2
    nc = 3
    fami = 'RIGI'
!
    if ((nomte .ne. 'MECA_BARRE') .and. (nomte .ne. 'MECA_2D_BARRE')) then
        ch16 = nomte
        call utmess('F', 'ELEMENTS2_42', sk=ch16)
    end if
!
!     --- RECUPERATION DES CARACTERISTIQUES MATERIAUX ---
    call jevech('PMATERC', 'L', lmater)
!
    call verift(fami, 1, 1, '+', zi(lmater), epsth_=epsth)
!
    r8bid = 0.0d0
    call rcvalb(fami, 1, 1, '+', zi(lmater), ' ', 'ELAS', 0, ' ', [r8bid], 1, 'E', e, codres, 1)
    if (epsth .ne. 0.d0) lteimp = .true.
!
!   Longueur de l'élément
    if (nomte .eq. 'MECA_BARRE') then
        xl = lonele()
    else if (nomte .eq. 'MECA_2D_BARRE') then
        xl = lonele(dime=2)
    end if
!
!     --- RECUPERATION DES CARACTERISTIQUES GENERALES DES SECTIONS ---
    aire = 0.0d0
    if (option .ne. 'EPSI_ELGA') then
        call jevech('PCAGNBA', 'L', lsect)
        aire = zr(lsect)
    end if
!
!     --- RECUPERATION DES ORIENTATIONS ALPHA,BETA,GAMMA ---
    call jevech('PCAORIE', 'L', lorien)
!     --- MATRICE DE ROTATION PGL
    call matrot(zr(lorien), pgl)
!
!     --- RECUPERATION DES DEPLACEMENTS OU DES VITESSES ----
    do ii = 1, 6
        ugr(ii) = 0.d0
    end do
!
    if (option .ne. 'ECIN_ELEM') then
!       ON RECUPERE DES DEPLACEMENTS
        call jevech('PDEPLAR', 'L', jdepl)
        if (nomte .eq. 'MECA_BARRE') then
            do ii = 1, 6
                ugr(ii) = zr(jdepl+ii-1)
            end do
        else if (nomte .eq. 'MECA_2D_BARRE') then
            ugr(1) = zr(jdepl+1-1)
            ugr(2) = zr(jdepl+2-1)
            ugr(4) = zr(jdepl+3-1)
            ugr(5) = zr(jdepl+4-1)
        end if
    else
        stopz = 'ONO'
        call tecach(stopz, 'PVITESR', 'L', iret, iad=jvite)
!       IRET NE PEUT VALOIR QUE 0 (TOUT EST OK) OU 2 (CHAMP NON FOURNI)
        if (iret .eq. 0) then
!           ON RECUPERE DES VITESSES
            if (nomte .eq. 'MECA_BARRE') then
                do ii = 1, 6
                    ugr(ii) = zr(jvite+ii-1)
                end do
            else if (nomte .eq. 'MECA_2D_BARRE') then
                ugr(1) = zr(jvite+1-1)
                ugr(2) = zr(jvite+2-1)
                ugr(4) = zr(jvite+3-1)
                ugr(5) = zr(jvite+4-1)
            end if
        else
!           ON RECUPERE DES DEPLACEMENTS
            call tecach(stopz, 'PDEPLAR', 'L', iret, iad=jdepl)
            if (iret .eq. 0) then
                if (nomte .eq. 'MECA_BARRE') then
                    do ii = 1, 6
                        ugr(ii) = zr(jdepl+ii-1)
                    end do
                else if (nomte .eq. 'MECA_2D_BARRE') then
                    ugr(1) = zr(jdepl+1-1)
                    ugr(2) = zr(jdepl+2-1)
                    ugr(4) = zr(jdepl+3-1)
                    ugr(5) = zr(jdepl+4-1)
                end if
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
!
        end if
!
    end if
!
!     --- VECTEUR DANS REPERE LOCAL  ULR = PGL * UGR
    call utpvgl(nno, nc, pgl, ugr, ulr)
!
!     --- RIGIDITE ELEMENTAIRE ---
    do ii = 1, 6
        do j = 1, 6
            klc(ii, j) = 0.d0
        end do
    end do
!
!     --- ENERGIE DE DEFORMATION ----
    if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', jende)
        xrig = e(1)*aire/xl
        klc(1, 1) = xrig
        klc(1, 4) = -xrig
        klc(4, 1) = -xrig
        klc(4, 4) = xrig
        iif = 0
        call ptenpo(6, ulr, klc, zr(jende), iif, iif)
!
        if (lteimp) then
            call ptenth(ulr, xl, epsth, 6, klc, enerth)
            zr(jende) = zr(jende)-enerth
        end if
!
    else if (option .eq. 'ECIN_ELEM') then
    call rcvalb(fami, 1, 1, '+', zi(lmater), ' ', 'ELAS', 0, ' ', [r8bid], 1, 'RHO', rho, codres, 1)
        call jevech('PENERCR', 'E', jende)
        call jevech('POMEGA2', 'L', jfreq)
        xmas = rho(1)*aire*xl/6.d0
        klc(1, 1) = xmas*2.d0
        klc(2, 2) = xmas*2.d0
        klc(3, 3) = xmas*2.d0
        klc(4, 4) = xmas*2.d0
        klc(5, 5) = xmas*2.d0
        klc(6, 6) = xmas*2.d0
        klc(1, 4) = xmas
        klc(4, 1) = xmas
        klc(2, 5) = xmas
        klc(5, 2) = xmas
        klc(3, 6) = xmas
        klc(6, 3) = xmas
        iif = 0
        itype = 50
        kanl = 1
        call ptenci(6, ulr, klc, zr(jfreq), zr(jende), itype, kanl, iif)
!
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jdefo)
        zr(jdefo-1+1) = (ulr(4)-ulr(1))/xl
    else
        xrig = e(1)*aire/xl
        klc(1, 1) = xrig
        klc(1, 4) = -xrig
        klc(4, 1) = -xrig
        klc(4, 4) = xrig
!        --- VECTEUR EFFORT LOCAL  FLR = KLC * ULR
        call pmavec('ZERO', 6, klc, ulr, flr)
!
!        --- TENIR COMPTE DES EFFORTS DUS A LA DILATATION ---
        if (lteimp) then
!              --- CALCUL DES FORCES INDUITES ---
            xfl1 = -epsth*e(1)*aire
            xfl4 = -xfl1
            flr(1) = flr(1)-xfl1
            flr(4) = flr(4)-xfl4
        end if
!
        if (option .eq. 'SIEF_ELGA') then
            call jevech('PCONTRR', 'E', jeffo)
            zr(jeffo) = -flr(1)
        else
!           OPTION NON PROGRAMMEE
            ASSERT(.false.)
        end if
    end if
!
end subroutine
