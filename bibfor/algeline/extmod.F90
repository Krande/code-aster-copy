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
subroutine extmod(basemo, numddl, nume, nbnumo, dmode, &
                  nbeq, nbnoe, iddl, nbddl)
    implicit none
! EXTRAIRE D'UN CONCEPT MODE_MECA LA DEFORMEE POUR UN OU PLUSIEURS DDL
! LES LAGRANGES SONT SUPPRIMES.
!-----------------------------------------------------------------------
! IN  :BASEMO : CONCEPT DE TYPE MODE_MECA
! IN  :NUMDDL : PERMET D'ACCEDER AU PROFIL DU CHAMP_NO EXTRAIT
! IN  :NUME   : LISTE DES NUMEROS D'ORDRE DES MODES CONSIDERES
! IN  :NBNUMO : NB DE MODES CONSIDERES
! IN  :NBEQ   : NB D'EQUATIONS
! IN  :NBNOE  : NB DE NOEUDS DU MAILLAGE
! IN  :IDDL   : LISTE DES INDICES DES DDL A EXTRAIRE
! IN  :NBDDL  : NB DE DDLS A EXTRAIRE
! OUT :DMODE  : VECTEUR => DEFORMEES MODALES
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
    integer(kind=8) :: nbddl, nbnumo, nbnoe, nume(nbnumo), iddl(nbddl)
    real(kind=8) :: dmode(nbddl*nbnoe*nbnumo)
    character(len=8) :: basemo
    character(len=14) :: numddl
    character(len=24) :: deeq, nomcha
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadmod, icm, ideeq, inumo, ipm, iret
    integer(kind=8) :: j, k, nbeq
!-----------------------------------------------------------------------
    call jemarq()
    deeq = numddl//'.NUME.DEEQ'
    call jeveuo(deeq, 'L', ideeq)
    ipm = 0
    icm = 0
    do i = 1, nbnumo
        inumo = nume(i)
        call rsexch('F', basemo, 'DEPL', inumo, nomcha, &
                    iret)
        nomcha = nomcha(1:19)//'.VALE'
        call jeveuo(nomcha, 'L', iadmod)
        ipm = ipm+icm
        icm = 0
        do j = 1, nbeq
            do k = 1, nbddl
                if (zi(ideeq+(2*j)-1) .eq. iddl(k)) then
                    icm = icm+1
                    dmode(ipm+icm) = zr(iadmod+j-1)
                    goto 22
                end if
            end do
22          continue
        end do
    end do
!
    call jedema()
end subroutine
