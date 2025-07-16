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
subroutine pj2dtr(corrMeshTemp, corrMesh, &
                  cellListType, cellListCode, &
                  geom1, geom2, &
                  spacedim, dala, &
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
#include "asterfort/pj4d2d.h"
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
    integer(kind=8), intent(in) :: cellListType(MT_NTYMAX), spacedim
    real(kind=8), intent(in) :: geom1(*), geom2(*)
    real(kind=8), intent(in) :: dala
    character(len=16), optional, intent(in)  :: listInterc_
    integer(kind=8), optional, intent(in)  :: nbInterc_
!
! --------------------------------------------------------------------------------------------------
!
!  but :
!    transformer corrMeshTemp en corrMesh en utilisant les fonc. de forme
!    des mailles du maillage1 (en 2d isoparametrique)
!
! --------------------------------------------------------------------------------------------------
!
!  in/jxin   corrMeshTemp   k16 : nom du corresp_2_mailla fait avec les tria3
!  in/jxout  corrMesh   k16 : nom du corresp_2_mailla final
!  in        cellListType i  : numeros des types de mailles 2d
!  in        cellListCode k8 : noms des types de mailles 2d
!  in        geom1        : geometrie des noeuds du maillage 1
!  in        geom2        : geometrie des noeuds du maillage 2
!  in        spacedim     : dimension de l'espace
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lext
    character(len=8) :: mesh1, mesh2, elrefa, nodeName2
    integer(kind=8) :: cnquad(3, 2)
    real(kind=8) :: crrefe(3, MT_NNOMAX), ff(MT_NNOMAX), cooele(3*MT_NNOMAX)
    integer(kind=8) :: nunos(MT_NNOMAX)
    real(kind=8) :: ksi, eta
    real(kind=8) :: x1, x2
    real(kind=8) :: xr1(2), xr2(2), xr3(2)
    real(kind=8) :: xg(2)
    integer(kind=8) :: nbCell1, nbCell2, nbNode1, nbNode2
    integer(kind=8) :: i2cocf, i2coco
    integer(kind=8) :: i2com1, i2conb, j2xxk1, i2conu
    integer(kind=8) :: ideca1, ideca2, ilcnx1
    integer(kind=8) :: ima1, ino, iNode2, kk, itypm, nbno, nno, itr
    integer(kind=8) :: iret, kdim, ndim
    integer(kind=8) :: nuno, nutm, nuno2
    integer(kind=8), parameter :: nbmax = 5
    integer(kind=8) :: tino2m(nbmax), nbnod, nbnodm
    real(kind=8) :: tdmin2(nbmax), disprj, distv
    aster_logical :: loin2
    integer(kind=8), pointer :: typmail(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: pjef_tr(:) => null()
    real(kind=8), pointer :: pjef_cf(:) => null()
    character(len=24), pointer :: pjxx_k1(:) => null()
    integer(kind=8), pointer :: tria3(:) => null()
    integer(kind=8), pointer :: vinterc(:) => null()
    integer(kind=8), pointer :: lino_loin(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

!   0. decoupage des quadrangles en 2 triangles (voir pj2dco)
!   ----------------------------------------------------------

    cnquad(1, 1) = 1
    cnquad(2, 1) = 2
    cnquad(3, 1) = 3

    cnquad(1, 2) = 1
    cnquad(2, 2) = 3
    cnquad(3, 2) = 4

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
    call jeveuo('&&PJXXCO.TRIA3', 'L', vi=tria3)
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
!       ITR : TRIA3 ASSOCIE A INO2
        itr = pjef_tr(iNode2)
        if (itr .eq. 0) cycle
!       IMA1 : MAILLE DE M1 ASSOCIE AU TRIA3 ITR
        ima1 = tria3(1+4*(itr-1)+4)
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
    nbnodm = 0

! - Create .pjef_nu .pjef_cf .pjef_co
    call wkvect(corrMesh//'.PJEF_NU', 'V V I', ideca2, i2conu)
    call wkvect(corrMesh//'.PJEF_CF', 'V V R', ideca2, i2cocf)
    call wkvect(corrMesh//'.PJEF_CO', 'V V R', 3*nbNode2, i2coco)
    ideca1 = 0
    ideca2 = 0
    do iNode2 = 1, nbNode2
!       ITR : TRIA3 ASSOCIE A INO2
        itr = pjef_tr(iNode2)
        if (itr .eq. 0) cycle
!       IMA1 : MAILLE DE M1 ASSOCIE AU TRIA3 ITR
        ima1 = tria3(1+4*(itr-1)+4)
!       ITYPM : TYPE DE LA MAILLE IMA1
        itypm = typmail(ima1)
        nutm = indiis(cellListType, itypm, 1, MT_NTYMAX)
        ASSERT(nutm .ne. 0)
        elrefa = cellListCode(nutm)
        nbno = zi(ilcnx1+ima1)-zi(ilcnx1-1+ima1)

        call elrfno(elrefa, ndim=ndim, nno=nno, nodeCoor=crrefe)
        ASSERT(nbno .eq. nno)

!       determination des coordonnees de iNode2 dans l'element de reference
        ksi = 0.d0
        eta = 0.d0
!       -- numero du 2eme noeud de ima1 : nuno2
        nuno2 = connex(1+zi(ilcnx1-1+ima1)-2+2)
!       si nuno2 est identique au 2eme noeud du tria3
!       c'est que le tria3 est en "dessous" :

        if (elrefa .eq. 'TR3' .or. elrefa .eq. 'TR6' .or. elrefa .eq. 'TR7') then
            do kk = 1, 3
                x1 = crrefe(1, kk)
                x2 = crrefe(2, kk)
                ksi = ksi+pjef_cf(ideca1+kk)*x1
                eta = eta+pjef_cf(ideca1+kk)*x2
            end do

        else if (elrefa .eq. 'QU4' .or. elrefa .eq. 'QU8' .or. elrefa .eq. 'QU9') then
            if (nuno2 .eq. tria3(1+4*(itr-1)+2)) then
!               SI 1ER TRIANGLE :
                do kk = 1, 3
                    x1 = crrefe(1, cnquad(kk, 1))
                    x2 = crrefe(2, cnquad(kk, 1))
                    ksi = ksi+pjef_cf(ideca1+kk)*x1
                    eta = eta+pjef_cf(ideca1+kk)*x2
                end do
            else
!               SI 2EME TRIANGLE :
                do kk = 1, 3
                    x1 = crrefe(1, cnquad(kk, 2))
                    x2 = crrefe(2, cnquad(kk, 2))
                    ksi = ksi+pjef_cf(ideca1+kk)*x1
                    eta = eta+pjef_cf(ideca1+kk)*x2
                end do
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
        xr1(1) = ksi
        xr1(2) = eta

!       -- on essaye d'ameliorer la precision de xr1(*) en utilisant reereg :
        if (spacedim == 3) then
!           -- en 3d on se place dans le plan de la maille
            itr = pjef_tr(iNode2)
            do ino = 1, nbno
                nunos(ino) = connex(1+zi(ilcnx1-1+ima1)-2+ino)
            end do
            call pj4d2d(tria3, itr, geom1, geom2(3*(iNode2-1)+1), nbno, &
                        nunos, cooele, xg)
        else
            do ino = 1, nbno
                nuno = connex(1+zi(ilcnx1-1+ima1)-2+ino)
                do kdim = 1, ndim
                    cooele(ndim*(ino-1)+kdim) = geom1(3*(nuno-1)+kdim)
                end do
            end do
            xg(1:2) = geom2(3*(iNode2-1)+1:3*(iNode2-1)+2)
        end if

        call reereg('C', elrefa, nno, cooele, xg, &
                    ndim, xr2, iret)

!       -- on regarde si iNode2 est exterieur a ima1 :
        call pjeflo(elrefa, ndim, iret, xr2, disprj)
        lext = (disprj .gt. 1.0d-02)

!       -- on choisit la meilleure approximation entre xr1 et xr2:
        call pjefmi(elrefa, nno, cooele, xg, ndim, &
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
                    ideca1 = ideca1+3
                    cycle
                end if
            end if
        end if

        zr(i2coco-1+3*(iNode2-1)+1) = xr3(1)
        zr(i2coco-1+3*(iNode2-1)+2) = xr3(2)

        call elrfvf(elrefa, xr3, ff)
        do ino = 1, nbno
            nuno = connex(1+zi(ilcnx1-1+ima1)-2+ino)
            zi(i2conu-1+ideca2+ino) = nuno
            zr(i2cocf-1+ideca2+ino) = ff(ino)
        end do

        ideca1 = ideca1+3
        ideca2 = ideca2+nbno
    end do

! - Alarm message
    if (loin2) then
        call pjloin(nbnod, nbnodm, mesh2, geom2, nbmax, tino2m, tdmin2, lino_loin)
    end if

    AS_DEALLOCATE(vi=lino_loin)

    call jedema()
end subroutine
