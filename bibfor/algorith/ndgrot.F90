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
subroutine ndgrot(sddyna, valinc, solalg, deldet, theta1, &
                  theta2, iran)
!
! person_in_charge: mickael.abbas at edf.fr
!
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/marota.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmchex.h"
#include "asterfort/promat.h"
#include "asterfort/proqua.h"
#include "asterfort/quavro.h"
#include "asterfort/transp.h"
#include "asterfort/vroqua.h"
    real(kind=8) :: theta2(3), theta1(3), deldet(3)
    character(len=19) :: sddyna
    character(len=19) :: solalg(*), valinc(*)
    integer(kind=8) :: iran(3)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - UTILITAIRE - DYNAMIQUE)
!
! MISE A JOUR DES VITESSES/ACCELERATIONS EN GRANDES ROTATIONS
! POUR POU_D_GD
!
! ----------------------------------------------------------------------
!
!
! IN  SDDYNA : SD DYNAMIQUE
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  THETA2 : VALEUR DE LA ROTATION PRECEDENTE
! IN  DELDET : INCREMENT DE ROTATION
! IN  IRAN   : NUMEROS ABSOLUS D'EQUATION DES DDL DE ROTATION DANS LES
!                 CHAM_NO
!
!
!
!
    character(len=19) :: depmoi, depplu, vitplu, accplu
    character(len=19) :: depkm1, vitkm1, acckm1, romk, romkm1
    integer(kind=8) :: ic
    real(kind=8) :: quapro(4), quarot(4), delqua(4)
    real(kind=8) :: qim(3), qikm1(3), qik(3), omkm1(3), ompkm1(3), delrot(3)
    real(kind=8) :: vect1(3), vect2(3), vect3(3), vect4(3), rotm(3, 3)
    real(kind=8) :: rotkm(3, 3), rotk(3, 3), rotmt(3, 3), rotkmt(3, 3)
    real(kind=8) :: coevit, coeacc
    character(len=19) :: depdel
    real(kind=8), pointer :: acckm(:) => null()
    real(kind=8), pointer :: accp(:) => null()
    real(kind=8), pointer :: depde(:) => null()
    real(kind=8), pointer :: depkm(:) => null()
    real(kind=8), pointer :: depm(:) => null()
    real(kind=8), pointer :: depp(:) => null()
    real(kind=8), pointer :: romkm(:) => null()
    real(kind=8), pointer :: vromk(:) => null()
    real(kind=8), pointer :: vitkm(:) => null()
    real(kind=8), pointer :: vitp(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- COEFFICIENTS
!
    coevit = ndynre(sddyna, 'COEF_VITE')
    coeacc = ndynre(sddyna, 'COEF_ACCE')
!
! --- DECOMPACTION VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPMOI', depmoi)
    call nmchex(valinc, 'VALINC', 'DEPPLU', depplu)
    call nmchex(valinc, 'VALINC', 'VITPLU', vitplu)
    call nmchex(valinc, 'VALINC', 'ACCPLU', accplu)
    call nmchex(valinc, 'VALINC', 'DEPKM1', depkm1)
    call nmchex(valinc, 'VALINC', 'VITKM1', vitkm1)
    call nmchex(valinc, 'VALINC', 'ACCKM1', acckm1)
    call nmchex(valinc, 'VALINC', 'ROMKM1', romkm1)
    call nmchex(valinc, 'VALINC', 'ROMK  ', romk)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
!
! --- RECUPERATION DES ADRESSES
!
    call jeveuo(depmoi(1:19)//'.VALE', 'L', vr=depm)
    call jeveuo(depplu(1:19)//'.VALE', 'E', vr=depp)
    call jeveuo(depdel(1:19)//'.VALE', 'E', vr=depde)
    call jeveuo(vitplu(1:19)//'.VALE', 'E', vr=vitp)
    call jeveuo(accplu(1:19)//'.VALE', 'E', vr=accp)
    call jeveuo(depkm1(1:19)//'.VALE', 'L', vr=depkm)
    call jeveuo(vitkm1(1:19)//'.VALE', 'L', vr=vitkm)
    call jeveuo(acckm1(1:19)//'.VALE', 'L', vr=acckm)
    call jeveuo(romkm1(1:19)//'.VALE', 'L', vr=romkm)
    call jeveuo(romk(1:19)//'.VALE', 'E', vr=vromk)
!
! --- QUATERNION DE L'INCREMENT DE ROTATION
!
    call vroqua(deldet, delqua)
!
! --- QUATERNION DE LA ROTATION PRECEDENTE
!
    call vroqua(theta1, quarot)
!
! --- CALCUL DE LA NOUVELLE ROTATION
!
    call proqua(delqua, quarot, quapro)
    call quavro(quapro, theta1)
!
! --- MISE A JOUR DES DEPLACEMENTS
!
    do ic = 1, 3
        depp(1+iran(ic)-1) = theta1(ic)
        depde(1+iran(ic)-1) = theta1(ic)
    end do
!
! --- QUATERNION DE LA ROTATION PRECEDENTE
!
    call vroqua(theta2, quarot)
!
! --- CALCUL DE LA NOUVELLE ROTATION
!
    call proqua(delqua, quarot, quapro)
    call quavro(quapro, theta2)
!
! --- MISE A JOUR DE LA ROTATION
!
    do ic = 1, 3
        vromk(1+iran(ic)-1) = theta2(ic)
    end do
!
! --- CALCUL DES INCREMENTS DE ROTATION
!
    do ic = 1, 3
        qim(ic) = depm(1+iran(ic)-1)
        qikm1(ic) = depkm(1+iran(ic)-1)
        qik(ic) = depp(1+iran(ic)-1)
        omkm1(ic) = vitkm(1+iran(ic)-1)
        ompkm1(ic) = acckm(1+iran(ic)-1)
    end do
!
! --- CALCUL DE L'INCREMENT DE ROTATION TOTALE
!
    do ic = 1, 3
        delrot(ic) = vromk(1+iran(ic)-1)-romkm(1+iran(ic)-1)
    end do
!
! --- CALCUL DES MATRICES DE ROTATION
!
    call marota(qim, rotm)
    call marota(qikm1, rotkm)
    call marota(qik, rotk)
    call transp(rotm, 3, 3, 3, rotmt, &
                3)
    call transp(rotkm, 3, 3, 3, rotkmt, &
                3)
!
! --- CALCUL DE LA VITESSE ANGULAIRE
!
    call promat(rotmt, 3, 3, 3, delrot, &
                3, 3, 1, vect3)
    call promat(rotk, 3, 3, 3, vect3, &
                3, 3, 1, vect2)
    call promat(rotkmt, 3, 3, 3, omkm1, &
                3, 3, 1, vect3)
    call promat(rotk, 3, 3, 3, vect3, &
                3, 3, 1, vect1)
    do ic = 1, 3
        vitp(1+iran(ic)-1) = vect1(ic)+coevit*vect2(ic)
    end do
!
! --- CALCUL DE L'ACCELERATION ANGULAIRE
!
    call promat(rotkmt, 3, 3, 3, ompkm1, &
                3, 3, 1, vect4)
    call promat(rotk, 3, 3, 3, vect4, &
                3, 3, 1, vect3)
    do ic = 1, 3
        accp(1+iran(ic)-1) = vect3(ic)+coeacc*vect2(ic)
    end do
!
    call jedema()
end subroutine
