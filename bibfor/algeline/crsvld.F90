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

subroutine crsvld(motfac, solveu, istop, nprec, &
                  epsmat, mixpre, kellag, kxfem)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jevtbl.h"
!
    integer(kind=8) :: istop, nprec
    real(kind=8) :: epsmat
    character(len=3) :: mixpre, kellag
    character(len=8) :: kxfem
    character(len=16) :: motfac
    character(len=19) :: solveu
!  BUT : REMPLISSAGE SD_SOLVEUR LDLT
!
! IN K19 SOLVEU  : NOM DU SOLVEUR DONNE EN ENTREE
! OUT    SOLVEU  : LE SOLVEUR EST CREE ET INSTANCIE
! IN  IN ISTOP   : PARAMETRE LIE AUX MOT-CLE STOP_SINGULIER
! IN  IN NPREC   :                           NPREC
! IN  R8 EPSMAT  :                           FILTRAGE_MATRICE
! IN  K3 MIXPRE  :                           MIXER_PRECISION
! IN  K3 KELLAG  :                           ELIM_LAGR
! IN  K8 KXFEM   :                           PRE_COND_XFEM
! ----------------------------------------------------------
!
!
!
!
    integer(kind=8) :: ibid
    character(len=8) :: renum
    integer(kind=8), pointer :: slvi(:) => null()
    real(kind=8), pointer :: slvr(:) => null()
    character(len=24), pointer :: slvk(:) => null()
!
!------------------------------------------------------------------
    call jemarq()
!
! --- LECTURES PARAMETRES DEDIES AU SOLVEUR
    call getvtx(motfac, 'RENUM', iocc=1, scal=renum, nbret=ibid)
    ASSERT(ibid .eq. 1)
!
! --- ON REMPLIT LA SD_SOLVEUR
    call jeveuo(solveu//'.SLVK', 'E', vk24=slvk)
    call jeveuo(solveu//'.SLVR', 'E', vr=slvr)
    call jeveuo(solveu//'.SLVI', 'E', vi=slvi)
!
    slvk(1) = 'LDLT'
    slvk(2) = 'XXXX'
    slvk(3) = 'XXXX'
    slvk(4) = renum
    slvk(5) = 'XXXX'
    slvk(6) = 'XXXX'
    slvk(7) = 'XXXX'
    slvk(8) = 'XXXX'
    slvk(9) = 'XXXX'
    slvk(10) = 'XXXX'
    slvk(11) = 'XXXX'
    slvk(12) = 'XXXX'
    slvk(13) = kellag
    slvk(14) = kxfem
!
    slvr(1) = 0.d0
    slvr(2) = 0.d0
    slvr(3) = jevtbl('TAILLE_BLOC')
    slvr(4) = 0.d0
    slvr(5) = 0.d0
!
    slvi(1) = nprec
    slvi(2) = -9999
    slvi(3) = istop
    slvi(4) = -9999
    slvi(5) = -9999
    slvi(6) = -9999
    slvi(7) = -9999
    slvi(8) = 0
!
    call jedema()
end subroutine
