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

subroutine majdva(numedd, sdnume, sddyna, valinc, solalg)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndgrot.h"
#include "asterfort/nmchex.h"
    character(len=24) :: numedd
    character(len=19) :: sddyna, sdnume
    character(len=19) :: solalg(*), valinc(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - UTILITAIRE - DYNAMIQUE)
!
! MET A JOUR LES ACCELERATIONS/VITESSES/ROTATIONS DANS
! LE CAS DES POUTRES EN GRANDES ROTATIONS
!
! ----------------------------------------------------------------------
!
!
! IN  NUMEDD : NUME_DDL
! IN  SDNUME : SD NUMEROTATION
! IN  SDDYNA : SD DYNAMIQUE
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
!
!
!
!
    character(len=19) :: vitplu, accplu
    character(len=19) :: ddepla, dvitla, daccla
    character(len=19) :: romk
    integer(kind=8) :: i, icomp, iran(3)
    integer(kind=8) :: neq
    real(kind=8) :: theta1(3), theta2(3), deldet(3)
    character(len=19) :: depplu, depdel
    character(len=19) :: depkm1, vitkm1, acckm1, romkm1
    real(kind=8), pointer :: acckm(:) => null()
    real(kind=8), pointer :: accp(:) => null()
    real(kind=8), pointer :: dacce(:) => null()
    real(kind=8), pointer :: ddepl(:) => null()
    real(kind=8), pointer :: depde(:) => null()
    real(kind=8), pointer :: depkm(:) => null()
    real(kind=8), pointer :: depp(:) => null()
    real(kind=8), pointer :: dvite(:) => null()
    real(kind=8), pointer :: romkm(:) => null()
    real(kind=8), pointer :: vromk(:) => null()
    real(kind=8), pointer :: vitkm(:) => null()
    real(kind=8), pointer :: vitp(:) => null()
    integer(kind=8), pointer :: ndro(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- POUTRES EN GRANDES ROTATIONS
!
    call jeveuo(sdnume//'.NDRO', 'L', vi=ndro)
!
! --- INITIALISATIONS
!
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
!
! --- DECOMPACTION VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPPLU', depplu)
    call nmchex(valinc, 'VALINC', 'VITPLU', vitplu)
    call nmchex(valinc, 'VALINC', 'ACCPLU', accplu)
    call nmchex(valinc, 'VALINC', 'DEPKM1', depkm1)
    call nmchex(valinc, 'VALINC', 'VITKM1', vitkm1)
    call nmchex(valinc, 'VALINC', 'ACCKM1', acckm1)
    call nmchex(valinc, 'VALINC', 'ROMKM1', romkm1)
    call nmchex(valinc, 'VALINC', 'ROMK  ', romk)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
    call nmchex(solalg, 'SOLALG', 'DDEPLA', ddepla)
    call nmchex(solalg, 'SOLALG', 'DVITLA', dvitla)
    call nmchex(solalg, 'SOLALG', 'DACCLA', daccla)
!
! --- RECUPERATION DES ADRESSES
!
    call jeveuo(ddepla(1:19)//'.VALE', 'L', vr=ddepl)
    call jeveuo(dvitla(1:19)//'.VALE', 'L', vr=dvite)
    call jeveuo(daccla(1:19)//'.VALE', 'L', vr=dacce)
    call jeveuo(depdel(1:19)//'.VALE', 'E', vr=depde)
    call jeveuo(depplu(1:19)//'.VALE', 'E', vr=depp)
    call jeveuo(vitplu(1:19)//'.VALE', 'E', vr=vitp)
    call jeveuo(accplu(1:19)//'.VALE', 'E', vr=accp)
    call jeveuo(depkm1(1:19)//'.VALE', 'E', vr=depkm)
    call jeveuo(vitkm1(1:19)//'.VALE', 'E', vr=vitkm)
    call jeveuo(acckm1(1:19)//'.VALE', 'E', vr=acckm)
    call jeveuo(romkm1(1:19)//'.VALE', 'E', vr=romkm)
    call jeveuo(romk(1:19)//'.VALE', 'L', vr=vromk)
!
! --- MISE A JOUR DEPL/VITE/ACCE
!
    icomp = 0
    do i = 1, neq
        if (ndro(i) .eq. 0) then
            depde(i) = depde(i)+ddepl(i)
            depp(i) = depp(i)+ddepl(i)
            vitp(i) = vitp(i)+dvite(i)
            accp(i) = accp(i)+dacce(i)
        else if (ndro(i) .eq. 1) then
            depkm(i) = depp(i)
            vitkm(i) = vitp(i)
            acckm(i) = accp(i)
            romkm(i) = vromk(i)
            icomp = icomp+1
            iran(icomp) = i
            deldet(icomp) = ddepl(i)
            theta1(icomp) = depp(i)
            theta2(icomp) = vromk(i)
            if (icomp .eq. 3) then
                icomp = 0
                call ndgrot(sddyna, valinc, solalg, deldet, theta1, &
                            theta2, iran)
            end if
        else
            ASSERT(.false.)
        end if
    end do
!
    call jedema()
end subroutine
