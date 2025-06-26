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
subroutine mnlcst(parcho, adime, ninc, nd, nchoc, &
                  h, hf, xcst)
    implicit none
!
!
!     MODE_NON_LINE PARTIE CONSTANTE G0
!     -    -               ---
! ----------------------------------------------------------------------
!
! REGROUPE LES TERMES CONSTANT DU PROBLEME A RESOUDRE
! ----------------------------------------------------------------------
! IN  PARCHO : K14  : NOM DE LA SD PARAMETRE DES CONTACTEURS
! IN  ADIME  : K14  : SD PARAMETRE POUR ADIMENSIONNEMENT
! IN  NINC   : I    : NOMBRE D INCONNUES DU SYSTEME
! IN  ND     : I    : NOMBRE DE DEGRES DE LIBERTE
! IN  NCHOC  : I    : NOMBRE DE CONTACTEURS
! IN  H      : I    : NOMBRE D'HARMONIQUES de X
! IN  HF     : I    : NOMBRE D'HARMONIQUES POUR F
! OUT XCST   : K14  : NOM DU VECTEUR DES TERMES CONSTANTS
! ----------------------------------------------------------------------
!
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "blas/dscal.h"
! ----------------------------------------------------------------------
! --- DECLARATION DES ARGUMENTS DE LA ROUTINE
! ----------------------------------------------------------------------
    integer(kind=8) :: ninc, nd, nchoc, h, hf
    character(len=14) :: parcho, adime, xcst
! ----------------------------------------------------------------------
! --- DECLARATION DES VARIABLES LOCALES
! ----------------------------------------------------------------------
    real(kind=8) :: alpha, eta, jeu
    integer(kind=8) :: iadim, icst, neqs, i
    real(kind=8), pointer :: jeumax(:) => null()
    real(kind=8), pointer :: raid(:) => null()
    character(len=8), pointer :: type(:) => null()
    real(kind=8), pointer :: orig(:) => null()
    real(kind=8), pointer :: vjeu(:) => null()
    integer(kind=8), pointer :: vneqs(:) => null()
    real(kind=8), pointer :: reg(:) => null()
    blas_int :: b_incx, b_n
!
    call jemarq()
! ----------------------------------------------------------------------
! --- RECUPERATION POINTEUR DE XVECT, L(XVECT) ET XCDL
! ----------------------------------------------------------------------
    call jeveuo(adime, 'L', iadim)
    call jeveuo(xcst, 'E', icst)
    b_n = to_blas_int(ninc-1)
    b_incx = to_blas_int(1)
    call dscal(b_n, 0.d0, zr(icst), b_incx)
! ----------------------------------------------------------------------
! --- EQUATION DE LA DYNAMIQUE
! ----------------------------------------------------------------------
!
! ----------------------------------------------------------------------
! --- EQUATION SUPPLEMENTAIRE POUR DEFINIR LA FORCE NON-LINEAIRE
! ----------------------------------------------------------------------
    call jeveuo(parcho//'.RAID', 'L', vr=raid)
    call jeveuo(parcho//'.REG', 'L', vr=reg)
    call jeveuo(parcho//'.NEQS', 'L', vi=vneqs)
    call jeveuo(parcho//'.TYPE', 'L', vk8=type)
    call jeveuo(parcho//'.ORIG', 'L', vr=orig)
    call jeveuo(parcho//'.JEU', 'L', vr=vjeu)
    call jeveuo(parcho//'.JEUMAX', 'L', vr=jeumax)
    neqs = 0
    do i = 1, nchoc
        alpha = raid(i)/zr(iadim)
        eta = reg(i)
        jeu = vjeu(i)/jeumax(1)
        if (type(i) (1:7) .eq. 'CERCLE') then
! ---     -ORIG1^2 - ORIG2^2
            zr(icst+nd*(2*h+1)+(neqs+2)*(2*hf+1)) = -( &
                                                    orig(1+3*(i-1))/jeu)**2-(orig(1+3*(i-1)+1 &
                                                                                  )/jeu &
                                                                             )**2
! ---     ETA
            zr(icst+nd*(2*h+1)+(neqs+3)*(2*hf+1)) = -eta
        else if (type(i) (1:6) .eq. 'PLAN') then
! ---     ETA
            zr(icst+nd*(2*h+1)+neqs*(2*hf+1)) = -eta
        end if
        neqs = neqs+vneqs(i)
    end do
! ----------------------------------------------------------------------
! --- AUTRES EQUATIONS
! ----------------------------------------------------------------------
!
    call jedema()
!
end subroutine
