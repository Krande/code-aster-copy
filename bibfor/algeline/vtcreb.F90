! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine vtcreb(field_nodez , base      , type_scalz,&
                  nume_ddlz   ,&
                  meshz       , prof_chnoz, idx_gdz, nb_equa_inz,&
                  nb_equa_outz, nbz, vchamz)
!
implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/wkvect.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeveuo.h"
#include "asterfort/sdchgd.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/asmpi_comm_vect.h"
!
!
    character(len=*), intent(in) :: field_nodez
    character(len=1), intent(in) :: base
    character(len=*), intent(in) :: type_scalz
    character(len=*), optional, intent(in) :: nume_ddlz
    character(len=*), optional, intent(in) :: meshz
    character(len=*), optional, intent(in) :: prof_chnoz
    integer, optional, intent(in) :: nb_equa_inz
    integer, optional, intent(in) :: idx_gdz
    integer, optional, intent(out) :: nb_equa_outz
    integer, optional, intent(in) :: nbz
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
!   In  prof_chno   : name of PROF_CHNO
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
    character(len=3) :: type_scal
    character(len=8) :: mesh
    character(len=19) :: prof_chno, field_node, chamno
    character(len=24) :: obj_refe, obj_vale, obj_desc
    character(len=24), pointer :: p_refe(:) => null()
    integer :: idx_gd, nb_equa, j_vale, ideb, ifin, i, pdesc_save, jvcham, nb_equa_gl
    integer, pointer :: p_desc(:) => null()
    aster_logical :: lchange, l_pmesh
!
! --------------------------------------------------------------------------------------------------
!
    field_node = field_nodez
    type_scal  = type_scalz
    if (present(nbz).and.present(vchamz)) then
      ideb=1
      ifin=nbz
      ASSERT(nbz.gt.1)
    else if (present(nbz).and..not.present(vchamz)) then
      ASSERT(.False.)
    else if (.not.present(nbz).and.present(vchamz)) then
      ASSERT(.False.)
    else
      ideb=1
      ifin=1
    endif
!
! - Get parameters from NUME_DDL
!
    if (present(nume_ddlz)) then
        call dismoi('NUM_GD_SI' , nume_ddlz, 'NUME_DDL', repi=idx_gd)
        call dismoi('NB_EQUA'   , nume_ddlz, 'NUME_DDL', repi=nb_equa)
        call dismoi('NOM_MAILLA', nume_ddlz, 'NUME_DDL', repk=mesh)
        call dismoi('PROF_CHNO' , nume_ddlz, 'NUME_DDL', repk=prof_chno)
    else
        idx_gd    = idx_gdz
        nb_equa   = nb_equa_inz
        prof_chno = prof_chnoz
        mesh      = meshz
    endif
!
    l_pmesh = isParallelMesh(mesh)
    nb_equa_gl = nb_equa
    ! J'enl??ve la v??rif car il y a des deadlock sinon
    ! if(l_pmesh) then
    !   call asmpi_comm_vect("MPI_SUM", "I", sci=nb_equa_gl)
    ! end if

    if(.not.l_pmesh .and. nb_equa_gl == 0) then
      ASSERT(ASTER_FALSE)
    end if

    if (ideb.eq.ifin) then
      obj_refe = field_node(1:19)//'.REFE'
      obj_vale = field_node(1:19)//'.VALE'
      obj_desc = field_node(1:19)//'.DESC'
!
! Create only one node FIELD
! - Object .REFE
      call wkvect(obj_refe, base//' V K24', 4, vk24 = p_refe)
      p_refe(1) = mesh
      p_refe(2) = prof_chno
! - Object .DESC
      call wkvect(obj_desc, base//' V I', 2, vi = p_desc)
      call jeecra(obj_desc, 'DOCU', cval='CHNO')
      p_desc(1) = idx_gd
      p_desc(2) = 1
! - Object .VALE
      call wkvect(obj_vale, base//' V '//type_scal, max(1, nb_equa), j_vale)
      call jeecra(obj_vale, "LONUTI", nb_equa)
      if (present(nb_equa_outz)) then
        nb_equa_outz = nb_equa
      endif
! - Change GRANDEUR
      if (type_scal.eq.'R'.or.type_scal.eq.'C'.or.type_scal.eq.'F') then
        call sdchgd(field_node, type_scal)
      endif
!
    else
!
! ! Create at the same time the (ifin-ideb+1) node FIELDs of name chamno(i)
      lchange=(type_scal.eq.'R'.or.type_scal.eq.'C'.or.type_scal.eq.'F')
      call jeveuo(vchamz,'L',jvcham)
      do i=ideb,ifin
        chamno=zk24(jvcham+i-1)(1:19)
        obj_refe = chamno(1:19)//'.REFE'
        obj_vale = chamno(1:19)//'.VALE'
        obj_desc = chamno(1:19)//'.DESC'
! - Object .REFE
        call wkvect(obj_refe, base//' V K24', 4, vk24 = p_refe)
        p_refe(1) = mesh
        p_refe(2) = prof_chno
! - Object .DESC
        call wkvect(obj_desc, base//' V I', 2, vi = p_desc)
        call jeecra(obj_desc, 'DOCU', cval='CHNO')
        p_desc(1) = idx_gd
        p_desc(2) = 1
! - Object .VALE
        call wkvect(obj_vale, base//' V '//type_scal, max(1, nb_equa), j_vale)
        call jeecra(obj_vale, "LONUTI", nb_equa)
! - Change GRANDEUR
        if (lchange) then
          if (i.eq.ideb) then
! CETTE ROUTINE CHANGE EVENTUELLEMENT LA VALEUR DE PDESC(1). ON LA SAUVEGARDE POUR LA TRANSFERER
! AUX AUTRES CHAMNO SANS SE POSER LES MEMES QUESTIONS QUI CONDUIRONT A LA MEME REPONSE: PDESC_SAVE
            call sdchgd(chamno, type_scal)
            pdesc_save=p_desc(1)
          else
            p_desc(1)=pdesc_save
          endif
        endif
      enddo
! DIVERS MUTUALISE
      if (present(nb_equa_outz)) then
        nb_equa_outz = nb_equa
      endif
    endif
!
end subroutine
