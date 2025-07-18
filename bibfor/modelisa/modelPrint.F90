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
subroutine modelPrint(model)
!
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/lxl_find.h"
#include "asterfort/int_to_char8.h"
!
    character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_MODELE
!
! Print model informations
!
! --------------------------------------------------------------------------------------------------
!
! In  model        : name of the model
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ifm, niv
    integer(kind=8) :: numvec, nbgrel, igrel, long_grel, nume_elem
    integer(kind=8) :: nume_type_poi1, jc, j, k, ibegin, isharp
    integer(kind=8) :: nume_type_elem, nume_type_geom
    integer(kind=8) :: iexi, nb_type_elem
    integer(kind=8) :: nb_elem_grel
    character(len=8) :: type_geom, name_entity, tabmai(8)
    character(len=19) :: ligrel_model
    character(len=8) :: mesh
    character(len=16) :: type_elem, modelisa, valk(8), formul
    character(len=32) :: phemod
    integer(kind=8), pointer :: p_mesh_typmai(:) => null()
    integer(kind=8), pointer :: p_nb_elem(:) => null()
    character(len=24) :: model_liel
    integer(kind=8), pointer :: p_model_liel(:) => null()
    character(len=8), pointer :: p_type_geom(:) => null()
    character(len=16), pointer :: p_modeli(:) => null()
    character(len=16), pointer :: p_formul(:) => null()
    character(len=16), pointer :: p_type_elem(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    numvec = 1
    modelisa = ' '
!
! - Access to mesh
!
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call jeveuo(mesh//'.TYPMAIL', 'L', vi=p_mesh_typmai)
!
! - Access to catalogs
!
    call jelira('&CATA.TE.NOMTE', 'NOMMAX', nb_type_elem)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), nume_type_poi1)
!
! - Access to model
!
    ligrel_model = model//'.MODELE'
    model_liel = model//'.MODELE    .LIEL'
!
! - Allocate
!
    AS_ALLOCATE(vi=p_nb_elem, size=nb_type_elem)
    AS_ALLOCATE(vk16=p_modeli, size=nb_type_elem)
    AS_ALLOCATE(vk16=p_formul, size=nb_type_elem)
    AS_ALLOCATE(vk8=p_type_geom, size=nb_type_elem)
    AS_ALLOCATE(vk16=p_type_elem, size=nb_type_elem)
!
! - Head
!
    call utmess('I', 'MODELE1_20')
!
! - Access to LIGREL
!
    call jeexin(ligrel_model//'.LIEL', iexi)
    if (iexi .ne. 0) then
        call jeveuo(model_liel, 'L', vi=p_model_liel)
        call jelira(model_liel, 'NMAXOC', nbgrel)
    end if
!
! - Access to "late" elements
!
!
! - Each type: counting
!
    do igrel = 1, nbgrel
        call jelira(jexnum(ligrel_model//'.LIEL', igrel), 'LONMAX', long_grel)
        nb_elem_grel = long_grel-1
        if (nb_elem_grel .ne. 0) then
            nume_type_elem = p_model_liel(numvec+nb_elem_grel)
            if (nume_type_elem .eq. 0) then
                goto 10
            end if
            ASSERT(nume_type_elem .gt. 0 .and. nume_type_elem .le. nb_type_elem)
            p_nb_elem(nume_type_elem) = p_nb_elem(nume_type_elem)+nb_elem_grel
            if (p_modeli(nume_type_elem) .eq. ' ') then
                call jenuno(jexnum('&CATA.TE.NOMTE', nume_type_elem), type_elem)
                if (type_elem .eq. 'MECA_HEXS8') then
                    call utmess('A', 'MODELE1_7')
                end if
                if (type_elem .eq. 'MECA_PYRAM13') then
                    call utmess('I', 'MODELE1_15')
                end if
                nume_elem = p_model_liel(numvec)
                if (nume_elem .lt. 0) then
                    nume_type_geom = nume_type_poi1
                else
                    nume_type_geom = p_mesh_typmai(nume_elem)
                end if
                call jenuno(jexnum('&CATA.TM.NOMTM', nume_type_geom), type_geom)
                call dismoi('PHEN_MODE', type_elem, 'TYPE_ELEM', repk=phemod)
                modelisa = phemod(17:32)
                if (phemod(1:10) .eq. '#PLUSIEURS') then
                    modelisa = ' '
                end if
                call dismoi('FORMULATION', type_elem, 'TYPE_ELEM', repk=formul)
                isharp = lxl_find(modelisa, '#')
                if (modelisa(1:10) .eq. '#PLUSIEURS') then
                    modelisa = ' '
                else
                    if (isharp .eq. 0) then
                        formul = ' '
                    else
                        modelisa = modelisa(1:isharp-1)
                    end if
                end if
! ------------- To have "right" format in logger
                if (type_elem .eq. ' ') type_elem = '_'
                if (modelisa .eq. ' ') modelisa = '_'
                if (formul .eq. ' ') formul = '_'
                if (type_geom .eq. ' ') type_geom = '_'
                p_type_elem(nume_type_elem) = type_elem
                p_modeli(nume_type_elem) = modelisa
                p_formul(nume_type_elem) = formul
                p_type_geom(nume_type_elem) = type_geom
            end if
        end if
        numvec = numvec+long_grel
10      continue
    end do
!
! - Each type: printing
!
    do igrel = 1, nb_type_elem
        nb_elem_grel = p_nb_elem(igrel)
        if (nb_elem_grel .ne. 0) then
            valk(1) = p_modeli(igrel)
            valk(2) = p_formul(igrel)
            valk(3) = p_type_geom(igrel)
            valk(4) = p_type_elem(igrel)
            call utmess('I', 'MODELE1_21', nk=4, valk=valk, si=nb_elem_grel)
        end if
    end do
!
! - Level 2 printing
!
    numvec = 1
    modelisa = ' '
    if (niv .eq. 2) then
        do igrel = 1, nbgrel
            call jelira(jexnum(ligrel_model//'.LIEL', igrel), 'LONMAX', long_grel)
            nb_elem_grel = long_grel-1
            if (nb_elem_grel .ne. 0) then
                nume_elem = p_model_liel(numvec)
                if (nume_elem .lt. 0) then
                    nume_type_geom = nume_type_poi1
                else
                    nume_type_geom = p_mesh_typmai(nume_elem)
                end if
                nume_type_elem = p_model_liel(numvec+nb_elem_grel)
                call jenuno(jexnum('&CATA.TM.NOMTM', nume_type_geom), type_geom)
                call jenuno(jexnum('&CATA.TE.NOMTE', nume_type_elem), type_elem)
                call dismoi('MODELISATION', type_elem, 'TYPE_ELEM', repk=modelisa)
                isharp = lxl_find(modelisa, '#')
                if (modelisa(1:10) .eq. '#PLUSIEURS') then
                    modelisa = ' '
                else
                    if (isharp .eq. 0) then
                        formul = ' '
                    else
                        modelisa = modelisa(1:isharp-1)
                        call dismoi('FORMULATION', type_elem, 'TYPE_ELEM', repk=formul)
                    end if
                end if
                valk(1) = modelisa
                valk(2) = formul
                valk(3) = type_geom
                valk(4) = type_elem
                if (nume_elem .lt. 0) then
                    call utmess('I', 'MODELE1_8')
                else
                    call utmess('I', 'MODELE1_9')
                end if
                call utmess('I', 'MODELE1_21', nk=4, valk=valk, si=long_grel-1)
                jc = 0
                ibegin = numvec
                do k = 1, nb_elem_grel
                    j = ibegin-1+k
                    jc = jc+1
                    nume_elem = p_model_liel(j)
                    ASSERT(nume_elem .gt. 0)
                    name_entity = int_to_char8(nume_elem)
                    if (jc .le. 8) then
                        tabmai(jc) = name_entity
                    end if
                    if (jc .eq. 8) then
                        call utmess('I', 'MODELE1_38', nk=8, valk=tabmai)
                        jc = 0
                    end if
                    if (k .eq. nb_elem_grel .and. jc .gt. 0) then
                        valk(1:8) = ' '
                        if (jc .le. 7) then
                            valk(1:jc) = tabmai(1:jc)
                        else
                            ASSERT(.false.)
                        end if
                        call utmess('I', 'MODELE1_38', nk=8, valk=valk)
                    end if
                end do
            end if
            numvec = numvec+long_grel
        end do
    end if
!
! - Deallocate
!
    AS_DEALLOCATE(vi=p_nb_elem)
    AS_DEALLOCATE(vk16=p_modeli)
    AS_DEALLOCATE(vk16=p_formul)
    AS_DEALLOCATE(vk8=p_type_geom)
    AS_DEALLOCATE(vk16=p_type_elem)
!
end subroutine
