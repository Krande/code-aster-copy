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
subroutine pj3dtr(corrMeshTemp, corrMesh, &
                  cellListType, cellListCode, &
                  geom1, geom2, &
                  dala, &
                  listInterc_, nbInterc_)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/indiis.h"
#include "asterfort/inslri.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/pjeflo.h"
#include "asterfort/pjefmi.h"
#include "asterfort/pjloin.h"
#include "asterfort/reereg.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/int_to_char8.h"
!
    character(len=16), intent(in) :: corrMesh, corrMeshTemp
    character(len=8), intent(in) :: cellListCode(MT_NTYMAX)
    integer(kind=8), intent(in) :: cellListType(MT_NTYMAX)
    real(kind=8), intent(in) :: geom1(*), geom2(*)
    real(kind=8), intent(in) :: dala
    character(len=16), optional, intent(in)  :: listInterc_
    integer(kind=8), optional, intent(in)  :: nbInterc_
!
! --------------------------------------------------------------------------------------------------
!
!  but :
!    transformer corrMeshTemp en corrMesh en utilisant les fonc. de forme
!    des mailles du maillage1 (en 3d isoparametrique)
!
! --------------------------------------------------------------------------------------------------
!
!  in/jxin   corrMeshTemp    k16 : nom du corresp_2_mailla fait avec les tetr4
!  in/jxout  corrMesh    k16 : nom du corresp_2_mailla final
!  in        cellListType i  : numeros des types de mailles
!  in        cellListCode k8 : noms des types de mailles
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lext
    character(len=8) :: mesh1, mesh2, elrefa, nodeName2
    integer(kind=8) :: cntetr(4, 1), cnpent(4, 3), cnhexa(4, 6), cnpyra(4, 2)
    real(kind=8) :: crrefe(3, MT_NNOMAX), ff(MT_NNOMAX), cooele(3*MT_NNOMAX)
    real(kind=8) :: ksi, eta, dzeta
    real(kind=8) :: x1, x2, x3
    real(kind=8) :: xr1(3), xr2(3), xr3(3)
    integer(kind=8) :: nbCell1, nbCell2, nbNode1, nbNode2
    integer(kind=8) :: i2cocf, i2coco
    integer(kind=8) :: i2com1, i2conb, j2xxk1, i2conu
    integer(kind=8) :: ideca1, ideca2, ilcnx1
    integer(kind=8) :: ima1, ino, iNode2, kk, ityp, itypm, nbno, nno, itr
    integer(kind=8) :: iret, kdim, ndim
    integer(kind=8) :: nuno, nutm
    integer(kind=8), parameter :: nbmax = 5
    integer(kind=8) :: tino2m(nbmax), nbnod, nbnodm
    real(kind=8) :: tdmin2(nbmax), disprj, distv
    aster_logical :: loin2
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: pjef_tr(:) => null()
    real(kind=8), pointer :: pjef_cf(:) => null()
    character(len=24), pointer :: pjxx_k1(:) => null()
    integer(kind=8), pointer :: tetr4(:) => null()
    integer(kind=8), pointer :: vinterc(:) => null()
    integer(kind=8), pointer :: lino_loin(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

!   0. DECOUPAGE DES ELEMENTS 3D EN TETRA (VOIR PJ3DC0) :
!   ----------------------------------------------------

!     0.1 : TETRAEDRE :
!     -----------------
    cntetr(1, 1) = 1
    cntetr(2, 1) = 2
    cntetr(3, 1) = 3
    cntetr(4, 1) = 4

!     0.2 : PENTAEDRE :
!     -----------------
    cnpent(1, 1) = 1
    cnpent(2, 1) = 3
    cnpent(3, 1) = 6
    cnpent(4, 1) = 5

    cnpent(1, 2) = 1
    cnpent(2, 2) = 6
    cnpent(3, 2) = 4
    cnpent(4, 2) = 5

    cnpent(1, 3) = 1
    cnpent(2, 3) = 3
    cnpent(3, 3) = 5
    cnpent(4, 3) = 2

!     0.3 : HEXAEDRE :
!     -----------------
    cnhexa(1, 1) = 1
    cnhexa(2, 1) = 4
    cnhexa(3, 1) = 8
    cnhexa(4, 1) = 6

    cnhexa(1, 2) = 1
    cnhexa(2, 2) = 8
    cnhexa(3, 2) = 6
    cnhexa(4, 2) = 5

    cnhexa(1, 3) = 1
    cnhexa(2, 3) = 4
    cnhexa(3, 3) = 6
    cnhexa(4, 3) = 2

    cnhexa(1, 4) = 2
    cnhexa(2, 4) = 4
    cnhexa(3, 4) = 8
    cnhexa(4, 4) = 7

    cnhexa(1, 5) = 2
    cnhexa(2, 5) = 8
    cnhexa(3, 5) = 6
    cnhexa(4, 5) = 7

    cnhexa(1, 6) = 2
    cnhexa(2, 6) = 4
    cnhexa(3, 6) = 7
    cnhexa(4, 6) = 3

!     0.4 : PYRAMIDE :
!     -----------------
    cnpyra(1, 1) = 1
    cnpyra(2, 1) = 2
    cnpyra(3, 1) = 3
    cnpyra(4, 1) = 5

    cnpyra(1, 2) = 1
    cnpyra(2, 2) = 3
    cnpyra(3, 2) = 4
    cnpyra(4, 2) = 5

!   General parameters
    call jeveuo(corrMeshTemp//'.PJXX_K1', 'L', vk24=pjxx_k1)
    call jeveuo(corrMeshTemp//'.PJEF_CF', 'L', vr=pjef_cf)
    call jeveuo(corrMeshTemp//'.PJEF_TR', 'L', vi=pjef_tr)
    mesh1 = pjxx_k1(1) (1:8)
    mesh2 = pjxx_k1(2) (1:8)
    call dismoi('NB_NO_MAILLA', mesh1, 'MAILLAGE', repi=nbNode1)
    call dismoi('NB_NO_MAILLA', mesh2, 'MAILLAGE', repi=nbNode2)
    call dismoi('NB_MA_MAILLA', mesh1, 'MAILLAGE', repi=nbCell1)
    call dismoi('NB_MA_MAILLA', mesh2, 'MAILLAGE', repi=nbCell2)
    call jeveuo('&&PJXXCO.TETR4', 'L', vi=tetr4)
    call jeveuo(mesh1//'.CONNEX', 'L', vi=connex)
    call jeveuo(jexatr(mesh1//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
    call jeveuo(mesh1//'.TYPMAIL', 'L', vi=typmail)

!   -- l'objet lino_loin contiendra la liste des noeuds projetes un peu loin
    AS_ALLOCATE(vi=lino_loin, size=nbNode2)

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
!       ITR : TETR4 ASSOCIE A INO2
        itr = pjef_tr(iNode2)
        if (itr .eq. 0) cycle
!       IMA1 : MAILLE DE M1 ASSOCIE AU TETR4 ITR
        ima1 = tetr4(1+6*(itr-1)+5)
        nbno = zi(ilcnx1+ima1)-zi(ilcnx1-1+ima1)
        zi(i2conb-1+iNode2) = nbno
        zi(i2com1-1+iNode2) = ima1
        ideca2 = ideca2+nbno
    end do
    if (ideca2 .eq. 0) then
        call utmess('F', 'CALCULEL3_97')
    end if

    loin2 = .false.
    nbnod = 0
    tdmin2 = 0.d0
    tino2m = 0
    nbnodm = 0

! - Create .pjef_nu .pjef_cf .pjef_co
    call wkvect(corrMesh//'.PJEF_NU', 'V V I', ideca2, i2conu)
    call wkvect(corrMesh//'.PJEF_CF', 'V V R', ideca2, i2cocf)
    call wkvect(corrMesh//'.PJEF_CO', 'V V R', 3*nbNode2, i2coco)
    ideca1 = 0
    ideca2 = 0
    do iNode2 = 1, nbNode2
!       ITR : TETR4 ASSOCIE A INO2
        itr = pjef_tr(iNode2)
        if (itr .eq. 0) cycle
!       IMA1 : MAILLE DE M1 ASSOCIE AU TETR4 ITR
        ima1 = tetr4(1+6*(itr-1)+5)
!       ITYPM : TYPE DE LA MAILLE IMA1
        itypm = typmail(ima1)
        nutm = indiis(cellListType, itypm, 1, MT_NTYMAX)
        ASSERT(nutm .ne. 0)
        elrefa = cellListCode(nutm)
        ityp = tetr4(1+6*(itr-1)+6)
        nbno = zi(ilcnx1+ima1)-zi(ilcnx1-1+ima1)

        call elrfno(elrefa, ndim=ndim, nno=nno, nodeCoor=crrefe)
        ASSERT(nbno .eq. nno)

!       determination des coordonnees de iNode2 dans l'element de reference
        ksi = 0.d0
        eta = 0.d0
        dzeta = 0.d0

        if (elrefa .eq. 'TE4' .or. elrefa .eq. 'T10') then
            do kk = 1, 4
                x1 = crrefe(1, cntetr(kk, ityp))
                x2 = crrefe(2, cntetr(kk, ityp))
                x3 = crrefe(3, cntetr(kk, ityp))
                ksi = ksi+pjef_cf(ideca1+kk)*x1
                eta = eta+pjef_cf(ideca1+kk)*x2
                dzeta = dzeta+pjef_cf(ideca1+kk)*x3
            end do

        else if (elrefa .eq. 'PE6' .or. elrefa .eq. 'P15' .or. elrefa .eq. 'P18') then
            do kk = 1, 4
                x1 = crrefe(1, cnpent(kk, ityp))
                x2 = crrefe(2, cnpent(kk, ityp))
                x3 = crrefe(3, cnpent(kk, ityp))
                ksi = ksi+pjef_cf(ideca1+kk)*x1
                eta = eta+pjef_cf(ideca1+kk)*x2
                dzeta = dzeta+pjef_cf(ideca1+kk)*x3
            end do

        else if (elrefa .eq. 'HE8' .or. elrefa .eq. 'H20' .or. &
                 elrefa .eq. 'H27' .or. elrefa .eq. 'HE9') then
            do kk = 1, 4
                x1 = crrefe(1, cnhexa(kk, ityp))
                x2 = crrefe(2, cnhexa(kk, ityp))
                x3 = crrefe(3, cnhexa(kk, ityp))
                ksi = ksi+pjef_cf(ideca1+kk)*x1
                eta = eta+pjef_cf(ideca1+kk)*x2
                dzeta = dzeta+pjef_cf(ideca1+kk)*x3
            end do

        else if (elrefa .eq. 'PY5' .or. elrefa .eq. 'P13') then
            do kk = 1, 4
                x1 = crrefe(1, cnpyra(kk, ityp))
                x2 = crrefe(2, cnpyra(kk, ityp))
                x3 = crrefe(3, cnpyra(kk, ityp))
                ksi = ksi+pjef_cf(ideca1+kk)*x1
                eta = eta+pjef_cf(ideca1+kk)*x2
                dzeta = dzeta+pjef_cf(ideca1+kk)*x3
            end do

        else
            ASSERT(ASTER_FALSE)
        end if
        xr1(1) = ksi
        xr1(2) = eta
        xr1(3) = dzeta

!       -- on essaye d'ameliorer la precision de xr1(*) en utilisant reereg :
        do ino = 1, nbno
            nuno = connex(1+zi(ilcnx1-1+ima1)-2+ino)
            do kdim = 1, ndim
                cooele(ndim*(ino-1)+kdim) = geom1(3*(nuno-1)+kdim)
            end do
        end do
        call reereg('C', elrefa, nno, cooele, geom2(3*(iNode2-1)+1), &
                    ndim, xr2, iret)

!       -- on regarde si iNode2 est exterieur a ima1 :
        call pjeflo(elrefa, ndim, iret, xr2, disprj)
        lext = (disprj .gt. 1.0d-02)

!       -- on choisit la meilleure approximation entre xr1 et xr2:
        call pjefmi(elrefa, nno, cooele, geom2(3*(iNode2-1)+1), ndim, &
                    xr1, xr2, lext, xr3, distv)

        if (distv .lt. dala) then
            lext = .false.
        else
            if (nint(disprj) .eq. 999) then
                loin2 = .true.
            else if (disprj .gt. 1.0d-01) then
!               on regarde si le noeud est deja projete par une autre
!               occurrence de VIS_A_VIS
                if (present(nbInterc_)) then
                    if (nbInterc_ .ne. 0) then
                        call jeveuo(listInterc_, 'L', vi=vinterc)
                    end if
                    do ino = 1, nbInterc_
                        if (iNode2 .eq. vinterc(ino)) then
                            zi(i2conb-1+iNode2) = 0
                            nodeName2 = int_to_char8(iNode2)
                            call utmess('A', 'CALCULEL5_47', si=vinterc(nbInterc_+1), sk=nodeName2)
                            exit
                        end if
                    end do
                end if
                if (zi(i2conb-1+iNode2) .ne. 0) then
                    loin2 = .true.
                    nbnodm = nbnodm+1
                    lino_loin(nbnodm) = iNode2
                    call inslri(nbmax, nbnod, tdmin2, tino2m, distv, iNode2)
                else
                    ideca1 = ideca1+4
                    cycle
                end if
            end if
        end if

        zr(i2coco-1+3*(iNode2-1)+1) = xr3(1)
        zr(i2coco-1+3*(iNode2-1)+2) = xr3(2)
        zr(i2coco-1+3*(iNode2-1)+3) = xr3(3)

        call elrfvf(elrefa, xr3, ff)
        do ino = 1, nbno
            nuno = connex(1+zi(ilcnx1-1+ima1)-2+ino)
            zi(i2conu-1+ideca2+ino) = nuno
            zr(i2cocf-1+ideca2+ino) = ff(ino)
        end do

        ideca1 = ideca1+4
        ideca2 = ideca2+nbno
    end do

! - Alarm message
    if (loin2) then
        call pjloin(nbnod, nbnodm, mesh2, geom2, nbmax, tino2m, tdmin2, lino_loin)
    end if

    AS_DEALLOCATE(vi=lino_loin)

    call jedema()
end subroutine
