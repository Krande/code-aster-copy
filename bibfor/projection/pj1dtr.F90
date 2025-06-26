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
subroutine pj1dtr(corrMeshTemp, corrMesh, cellListType, cellListCode)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/indiis.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=16), intent(in) :: corrMesh, corrMeshTemp
    character(len=8), intent(in) :: cellListCode(MT_NTYMAX)
    integer(kind=8), intent(in) :: cellListType(MT_NTYMAX)
!
! --------------------------------------------------------------------------------------------------
!
!  but :
!    transformer corrMeshTemp en corrMesh en utilisant les fonc. de forme
!    des mailles du maillage1 (en 1d isoparametrique)
!
! --------------------------------------------------------------------------------------------------
!
!  in/jxin   corrMeshTemp   k16 : nom du corresp_2_mailla fait avec les tria3
!  in/jxout  corrMesh   k16 : nom du corresp_2_mailla final
!  in        cellListType i  : numeros des types de mailles 1d
!  in        cellListCode k8 : noms des types de mailles 1d
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: mesh1, mesh2, elrefa
    real(kind=8) :: crrefe(3, MT_NNOMAX), ff(MT_NNOMAX)
    real(kind=8) :: ksi
    real(kind=8) :: x1
    real(kind=8) :: x(1)
    integer(kind=8) :: i2cocf, i2coco
    integer(kind=8) :: nbCell1, nbCell2, nbNode1, nbNode2
    integer(kind=8) :: i2com1, i2conb, j2xxk1, i2conu
    integer(kind=8) :: ideca1, ideca2, ilcnx1
    integer(kind=8) :: ima1, ino, iNode2, kk, itypm, nbno, nno, itr
    integer(kind=8) :: nuno, nutm
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: pjef_tr(:) => null()
    real(kind=8), pointer :: pjef_cf(:) => null()
    character(len=24), pointer :: pjxx_k1(:) => null()
    integer(kind=8), pointer :: seg2(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - General parameters
    call jeveuo(corrMeshTemp//'.PJXX_K1', 'L', vk24=pjxx_k1)
    call jeveuo(corrMeshTemp//'.PJEF_CF', 'L', vr=pjef_cf)
    call jeveuo(corrMeshTemp//'.PJEF_TR', 'L', vi=pjef_tr)
    mesh1 = pjxx_k1(1) (1:8)
    mesh2 = pjxx_k1(2) (1:8)
    call dismoi('NB_NO_MAILLA', mesh1, 'MAILLAGE', repi=nbNode1)
    call dismoi('NB_NO_MAILLA', mesh2, 'MAILLAGE', repi=nbNode2)
    call dismoi('NB_MA_MAILLA', mesh1, 'MAILLAGE', repi=nbCell1)
    call dismoi('NB_MA_MAILLA', mesh2, 'MAILLAGE', repi=nbCell2)
    call jeveuo('&&PJXXCO.SEG2', 'L', vi=seg2)
    call jeveuo(mesh1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    call jeveuo(mesh1//'.TYPMAIL', 'L', vi=typmail)

! - Allocate corrMesh
    call wkvect(corrMesh//'.PJXX_K1', 'V V K24', 5, j2xxk1)
    zk24(j2xxk1-1+1) = mesh1
    zk24(j2xxk1-1+2) = mesh2
    zk24(j2xxk1-1+3) = 'COLLOCATION'

! - Create .pjef_nb and .pjef_m1
    call wkvect(corrMesh//'.PJEF_NB', 'V V I', nbNode2, i2conb)
    call wkvect(corrMesh//'.PJEF_M1', 'V V I', nbNode2, i2com1)
    ideca2 = 0
    do iNode2 = 1, nbNode2
!       ITR : SEG2 ASSOCIE A INO2
        itr = pjef_tr(iNode2)
        if (itr .eq. 0) cycle
!       IMA1 : MAILLE DE M1 ASSOCIEE AU SEG2 ITR
        ima1 = seg2(1+3*(itr-1)+3)
        nbno = zi(ilcnx1+ima1)-zi(ilcnx1-1+ima1)
        zi(i2conb-1+iNode2) = nbno
        zi(i2com1-1+iNode2) = ima1
        ideca2 = ideca2+nbno
    end do
    if (ideca2 .eq. 0) then
        call utmess('F', 'CALCULEL3_97')
    end if

! - Create .pjef_nu .pjef_cf .pjef_co
    call wkvect(corrMesh//'.PJEF_NU', 'V V I', ideca2, i2conu)
    call wkvect(corrMesh//'.PJEF_CF', 'V V R', ideca2, i2cocf)
    call wkvect(corrMesh//'.PJEF_CO', 'V V R', 3*nbNode2, i2coco)
    ideca1 = 0
    ideca2 = 0
    do iNode2 = 1, nbNode2
!       ITR : SEG2 ASSOCIE A INO2
        itr = pjef_tr(iNode2)
        if (itr .eq. 0) cycle
!       IMA1 : MAILLE DE M1 ASSOCIE AU SEG2 ITR
        ima1 = seg2(1+3*(itr-1)+3)
!       ITYPM : TYPE DE LA MAILLE IMA1
        itypm = typmail(ima1)
        nutm = indiis(cellListType, itypm, 1, MT_NTYMAX)
        ASSERT(nutm .ne. 0)
        elrefa = cellListCode(nutm)
        nbno = zi(ilcnx1+ima1)-zi(ilcnx1-1+ima1)
        call elrfno(elrefa, nno=nno, nodeCoor=crrefe)
        ASSERT(nbno .eq. nno)

!       determination des coordonnees de iNode2 dans l'element de reference
        ksi = 0.d0
        do kk = 1, 2
            x1 = crrefe(1, kk)
            ksi = ksi+pjef_cf(ideca1+kk)*x1
        end do
        x(1) = ksi
        zr(i2coco-1+3*(iNode2-1)+1) = x(1)

!       2.2.2 :
!       CALCUL DES F. DE FORME AUX NOEUDS POUR LE POINT KSI
!       -------------------------------------------------------
        call elrfvf(elrefa, x, ff)
        do ino = 1, nbno
            nuno = connex(1+zi(ilcnx1-1+ima1)-2+ino)
            zi(i2conu-1+ideca2+ino) = nuno
            zr(i2cocf-1+ideca2+ino) = ff(ino)
        end do

        ideca1 = ideca1+2
        ideca2 = ideca2+nbno

    end do

    call jedema()
end subroutine
