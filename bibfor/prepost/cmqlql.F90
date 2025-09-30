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

subroutine cmqlql(main, maout, nbma, lima)
    use crea_maillage_module
!
    implicit none
#include "asterf_types.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
    integer(kind=8) :: nbma, lima(nbma)
    character(len=8) :: main, maout
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES QUADRATIQUES -> LINEAIRE
!-----------------------------------------------------------------------
!               SEG3            --> SEG2
!               TRIA6           --> TRIA3,
!               QUAD8,QUAD9     --> QUAD4,
!               TETRA10         --> TETRA4
!               PYRAM13         --> PYRMA5
!               PENTA15,PENTA18 --> PENTA6
!               HEXA20,HEXA27   --> HEXA8
! ----------------------------------------------------------------------
! IN        MAIN   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  MAOUT  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO DES MAILLES A TRAITER
! ----------------------------------------------------------------------
!
!
    type(Mmesh) :: mesh_conv
    character(len=8) :: conv_type(2)
!
    call jemarq()
!
! - Create new mesh
!
    call mesh_conv%init(main, convert_max=ASTER_FALSE)
!
! - Add conversion
!
    conv_type = ["SEG3", "SEG2"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA6", "TRIA3"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA7", "TRIA3"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD8", "QUAD4"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD9", "QUAD4"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA10", "TETRA4 "]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA15", "PENTA6 "]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA18", "PENTA6 "]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PYRAM13", "PYRAM5 "]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA20", "HEXA8 "]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA27", "HEXA8 "]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
!
! - Convert cells
!
    call mesh_conv%convert_cells(nbma, lima)
!
! - Copy mesh
!
    call mesh_conv%copy_mesh(maout)
!
! - Cleaning
!
    call mesh_conv%clean()
!
    call jedema()
!
end subroutine
