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

subroutine vtcreb(field_nodez, base, type_scalz, &
                  nume_ddlz, &
                  meshz, nume_equaz, idx_gdz, nb_equa_inz, &
                  nb_equa_outz, nbz, vchamz)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/wkvect.h"
#include "asterfort/dismoi.h"
#include "asterfort/gnomsd.h"
#include "asterfort/copisd.h"
#include "asterfort/exisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/codent.h"
#include "asterfort/idensd.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/asmpi_comm_vect.h"
!
!
    character(len=*), intent(in) :: field_nodez
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: type_scalz
    character(len=*), optional, intent(in) :: nume_ddlz
    character(len=*), optional, intent(in) :: meshz
    character(len=*), optional, intent(in) :: nume_equaz
    integer(kind=8), optional, intent(in) :: nb_equa_inz
    integer(kind=8), optional, intent(in) :: idx_gdz
    integer(kind=8), optional, intent(out) :: nb_equa_outz
    integer(kind=8), optional, intent(in) :: nbz
    character(len=24), optional, intent(in) :: vchamz
!
! --------------------------------------------------------------------------------------------------
!
! Field utility
!
! Create NODE field
!
! --------------------------------------------------------------------------------------------------
!
! In  field_node    : name of field
! In  base          : JEVEUX base to create field
! In  type_scal     : type of GRANDEUR (real or complex)
! With numbering:
!   In  nume_ddl    : name of numbering
! With complete informations:
!   In  mesh        : name of mesh
!   In  nume_equa   : name of NUME_EQUA
!   In  idx_gd      : index of GRANDEUR
!   In  nb_equa_in  : number of equations
! Create simultaneously nbz NODE fields of name vchamz(1)...vchamz(nbz)
!   In  nbz         : number of fields
!   In  vchamz      : vector of names of fields
!
!   Out nb_equa_out : number of equations
!
! --------------------------------------------------------------------------------------------------
!
    character(len=3) :: type_scal, type_scal2
    character(len=8) :: mesh, nomgd
    character(len=19) :: nume_equa, field_node, chamno, nume_equa_tmp
    character(len=24) :: obj_refe, obj_vale, noojb
    character(len=24), pointer :: p_refe(:) => null()
    integer(kind=8) :: idx_gd, nb_equa, j_vale, ideb, ifin, i, jvcham, nb_equa_gl, jrefn, prev, iexi
    aster_logical :: l_pmesh
!
! --------------------------------------------------------------------------------------------------
!
    field_node = field_nodez
    type_scal = type_scalz
    if (present(nbz) .and. present(vchamz)) then
        ideb = 1
        ifin = nbz
        ASSERT(nbz .gt. 1)
    else if (present(nbz) .and. .not. present(vchamz)) then
        ASSERT(.False.)
    else if (.not. present(nbz) .and. present(vchamz)) then
        ASSERT(.False.)
    else
        ideb = 1
        ifin = 1
    end if
!
! - Get parameters from NUME_DDL
!
    if (present(nume_ddlz)) then
        call dismoi('NUM_GD_SI', nume_ddlz, 'NUME_DDL', repi=idx_gd)
        call dismoi('NB_EQUA', nume_ddlz, 'NUME_DDL', repi=nb_equa)
        call dismoi('NOM_MAILLA', nume_ddlz, 'NUME_DDL', repk=mesh)
        call dismoi('NUME_EQUA', nume_ddlz, 'NUME_DDL', repk=nume_equa)
    else
        idx_gd = idx_gdz
        nb_equa = nb_equa_inz
        nume_equa = nume_equaz
        mesh = meshz
    end if
!
    call dismoi('TYPE_SCA', nume_equa, 'NUME_EQUA', repk=type_scal2)
    if (type_scal .ne. type_scal2) then
        call dismoi('NOM_GD', nume_equa, 'NUME_EQUA', repk=nomgd)
        nomgd(6:6) = type_scal(1:1)
        noojb = '12345678.NUME000000.PRNO'
        call gnomsd(field_node, noojb, 14, 19)
        noojb(1:8) = field_node(1:8)
        nume_equa_tmp = noojb(1:19)
        call copisd("NUME_EQUA", base, nume_equa, nume_equa_tmp)
        nume_equa = nume_equa_tmp
        call jeveuo(nume_equa//".REFN", 'E', jrefn)
        zk24(jrefn+1) = nomgd
        read (nume_equa(14:19), '(i6)') prev
        if (prev > 0) then
            prev = max(0, prev-1)
            nume_equa_tmp = nume_equa
            call codent(prev, "D0", nume_equa_tmp(14:19))
            call exisd('NUME_EQUA', nume_equa_tmp, iexi)
            if (iexi .gt. 0) then
                if (idensd('NUME_EQUA', nume_equa, nume_equa_tmp)) then
                    call detrsd("NUME_EQUA", nume_equa)
                    nume_equa = nume_equa_tmp
                end if
            end if
        end if
    end if
!
    l_pmesh = isParallelMesh(mesh)
    nb_equa_gl = nb_equa
    ! J'enlève la vérif car il y a des deadlock sinon
    ! if(l_pmesh) then
    !   call asmpi_comm_vect("MPI_SUM", "I", sci=nb_equa_gl)
    ! end if

    if (.not. l_pmesh .and. nb_equa_gl == 0) then
        ASSERT(ASTER_FALSE)
    end if

    if (ideb .eq. ifin) then
        obj_refe = field_node(1:19)//'.REFE'
        obj_vale = field_node(1:19)//'.VALE'
!
! Create only one node FIELD
! - Object .REFE
        call wkvect(obj_refe, base//' V K24', 4, vk24=p_refe)
        call jeecra(obj_refe, 'DOCU', cval='CHNO')
        p_refe(2) = nume_equa
! - Object .VALE
        call wkvect(obj_vale, base//' V '//type_scal, max(1, nb_equa), j_vale)
        call jeecra(obj_vale, "LONUTI", nb_equa)
        if (present(nb_equa_outz)) then
            nb_equa_outz = nb_equa
        end if
!
    else
!
! ! Create at the same time the (ifin-ideb+1) node FIELDs of name chamno(i)
        call jeveuo(vchamz, 'L', jvcham)
        do i = ideb, ifin
            chamno = zk24(jvcham+i-1) (1:19)
            obj_refe = chamno(1:19)//'.REFE'
            obj_vale = chamno(1:19)//'.VALE'
! - Object .REFE
            call wkvect(obj_refe, base//' V K24', 4, vk24=p_refe)
            call jeecra(obj_refe, 'DOCU', cval='CHNO')
            p_refe(2) = nume_equa
! - Object .VALE
            call wkvect(obj_vale, base//' V '//type_scal, max(1, nb_equa), j_vale)
            call jeecra(obj_vale, "LONUTI", nb_equa)
        end do
! DIVERS MUTUALISE
        if (present(nb_equa_outz)) then
            nb_equa_outz = nb_equa
        end if
    end if
!
end subroutine
