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

subroutine cmlqlq(main, maout, nbma, lima)
!
    use crea_maillage_module
!
    implicit none
#include "asterfort/jemarq.h"
#include "asterfort/jedema.h"
#include "jeveux.h"
!
    integer(kind=8) :: nbma, lima(nbma)
    character(len=8) :: main, maout
! ----------------------------------------------------------------------
!         TRANSFORMATION DES MAILLES LINEAIRES -> QUADRATIQUES
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
!     A PARTIR DU CATALOGUE TYPE_MAILLE__  :
!     REFERENCE     -->  NOUVEAU TYPE              NB DE NOEUDS
!
!                      REFTYP                         NBREF
!
!     1   POI1      -->  1                               1
!     2   SEG2      -->  4    ( SEG2 EN SEG3 )           3
!     3   SEG22     -->  3                               4
!     4   SEG3      -->  4                               3
!     5   SEG33     -->  5                               6
!     6   SEG4      -->  6                               4
!     7   TRIA3     -->  9    ( TRIA3 EN TRIA6 )         6
!     8   TRIA33    -->  8                               6
!     9   TRIA6     -->  9                               6
!     10  TRIA66    -->  10                             12
!     11  TRIA7     -->  11                              7
!     12  QUAD4     -->  14   ( QUAD4 EN QUAD8 )         8
!     13  QUAD44    -->  13                              8
!     14  QUAD8     -->  14                              8
!     15  QUAD88    -->  15                             16
!     16  QUAD9     -->  16                              9
!     17  QUAD99    -->  17                             18
!     18  TETRA4    -->  19   ( TETRA4 EN TETRA10 )     10
!     19  TETRA10   -->  19                             10
!     20  PENTA6    -->  21   ( PENTA6 EN PENTA15 )     15
!     21  PENTA15   -->  21                             15
!     22  PENTA18   -->  22                             18
!     23  PYRAM5    -->  24   ( PYRAM5 EN PYRAM13 )     13
!     24  PYRAM13   -->  24                             13
!     25  HEXA8     -->  26   ( HEXA8 EN HEXA20 )       20
!     26  HEXA20    -->  26                             20
!     27  HEXA27    -->  27                             27
!
! ----------------------------------------------------------------------
    call jemarq()
!
! - Create new mesh
!
    call mesh_conv%init(main)
!
! - Add conversion
!
    conv_type = ["SEG2", "SEG3"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TRIA3", "TRIA6"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["QUAD4", "QUAD8"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["TETRA4 ", "TETRA10"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PENTA6 ", "PENTA15"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["PYRAM5 ", "PYRAM13"]
    call mesh_conv%converter%add_conversion(conv_type(1), conv_type(2))
    conv_type = ["HEXA8 ", "HEXA20"]
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
