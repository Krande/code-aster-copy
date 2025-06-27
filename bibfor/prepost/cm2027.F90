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

subroutine cm2027(main, maout, nbma, lima)
!
    use crea_maillage_module
!
    implicit none
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
!
    integer(kind=8) :: nbma, lima(nbma)
    character(len=8) :: main, maout
!
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES HEXA20 HEXA27
! ----------------------------------------------------------------------
! IN        MAIN   K8  NOM DU MAILLAGE INITIAL
! IN/JXOUT  MAOUT  K8  NOM DU MAILLAGE TRANSFORME
! IN        NBMA    I  NOMBRE DE MAILLES A TRAITER
! IN        LIMA    I  NUMERO ET TYPE DES MAILLES A TRAITER
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
    call mesh_conv%init(main)
!
! - Add conversion
!
    conv_type = ["QUAD8", "QUAD9"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA20", "HEXA27"]
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
end subroutine
