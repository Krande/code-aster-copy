! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
subroutine cmhho(mesh_in, mesh_out, nb_list_elem, list_elem, prefix, ndinit)
!
use crea_maillage_module
!
implicit none
!
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
integer, intent(in) :: ndinit, nb_list_elem, list_elem(nb_list_elem)
character(len=8), intent(in) :: mesh_in, mesh_out
character(len=8), intent(in) :: prefix
!
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES POUR HHO
! ----------------------------------------------------------------------
! IN        mesh_in   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  mesh_out  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO ET TYPE DES MAILLES A TRAITER
! IN        PREFIX K8  PREFIXE DU NOM DES NOEUDS CREES (EX: N, NO, ...)
! IN        NDINIT  I  NUMERO INITIAL DES NOEUDS CREES
! ----------------------------------------------------------------------
!
type(Mmesh) :: mesh_hho
character(len=8) :: conv_type(2)

!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Create new mesh
!
    call mesh_hho%init(mesh_in)
!
! - Add conversion
!
    conv_type = ["SEG2", "SEG3"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA3", "TRIA7"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA6", "TRIA7"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD4", "QUAD9"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD8", "QUAD9"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA4 ", "TETRA15"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA10", "TETRA15"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA8 ", "HEXA27"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA20", "HEXA27"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PYRAM5 ", "PYRAM19"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PYRAM13", "PYRAM19"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA6 ", "PENTA21"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA15", "PENTA21"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA18", "PENTA21"]
    call mesh_hho%converter%add_conversion(conv_type(1), conv_type(2))
!
! - Convert cells
!
    call mesh_hho%convert_cells(nb_list_elem, list_elem, prefix, ndinit)
!
! - Copy mesh
!
    call mesh_hho%copy_mesh(mesh_out)
!
! - Create joints for parallel mesh
!
    call mesh_hho%create_joints(mesh_out)
!
! - Cleaning
!
    call mesh_hho%clean()
!
    call jedema()
end subroutine
