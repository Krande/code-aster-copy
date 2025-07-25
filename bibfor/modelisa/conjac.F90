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

subroutine conjac(i0, i1, i2, i3, macoc, &
                  nbcoc, mailla)
    implicit none
!
! DESCRIPTION : VERIFICATION D'UNE MAILLE FISSURE
!
!               HYPOTHESE (CF MODI_MAILLAGE) :
!               LA MAILLE DE FISSURE PEUT ETRE DU TYPE
!               - QUAD4 OU QUAD8 EN 2D
!               - PENTA6 OU PENTA15 OU HEXA8 OU HEXA20 EN 3D
!
! IN     : MAILLA : CHARACTER*8 , SCALAIRE
!                   NOM DU CONCEPT MAILLAGE
! IN     : NBCOC  : INTEGER , SCALAIRE
!                   NOMBRE DE NOEUDS DE LA MAILLE DE FISSURE
! IN/OUT : MACOC  : CHARACTER*8 , VECTEUR DE DIMENSION NBCOC+2
!                   MACOC(1) : NOM DE LA MAILLE DE FISSURE
!                   MACOC(2) : NOM DU TYPE DE LA MAILLE DE FISSURE
!                   MACOC(2+1) A MACOC(2+NBCOC) : NOMS DES NOEUDS DE
!                   LA MAILLE DE FISSURE, DANS L'ORDRE DEFINI PAR LA
!                   CONNECTIVITE. EN SORTIE L'ORDRE PEUT ETRE CHANGE
!                   SI LA MAILLE A ETE REORIENTEE.
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "jeveux.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: i0, i1, i2, i3
    character(len=8) :: mailla
    integer(kind=8) :: nbcoc
    character(len=8) :: macoc(2+nbcoc)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: jcoor, no0, no1, no2, no3, niv, ifm
    character(len=24) :: coorno
    real(kind=8) :: x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, scal
    real(kind=8) :: vnx, vny, vnz
    real(kind=8) :: vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3
!
! FONCTIONS EXTERNES
! ------------------
!     EXTERNAL      JEXNOM, JEXNUM, R8DOT, R8NRM2
!
! ROUTINES EXTERNES
! -----------------
!     EXTERNAL      JEDEMA, JEMARQ, JENONU, JEVEUO,
!                   DISMOI.
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
!-----------------------------------------------------------------------
!     INITIALISATIONS - ACCES AUX OBJETS DU CONCEPT MAILLAGE
!-----------------------------------------------------------------------
!
    coorno = mailla//'.COORDO    .VALE'
    call jeveuo(coorno, 'L', jcoor)
!
!     DETERMINATION D'UN VECTEUR NORMAL A LA MAILLE DE FISSURE
!
    no0 = char8_to_int(macoc(2+i0))
    x0 = zr(jcoor+3*(no0-1)+0)
    y0 = zr(jcoor+3*(no0-1)+1)
    z0 = zr(jcoor+3*(no0-1)+2)
!
    no1 = char8_to_int(macoc(2+i1))
    x1 = zr(jcoor+3*(no1-1)+0)
    y1 = zr(jcoor+3*(no1-1)+1)
    z1 = zr(jcoor+3*(no1-1)+2)
    vx1 = x1-x0
    vy1 = y1-y0
    vz1 = z1-z0
!
    no2 = char8_to_int(macoc(2+i2))
    x2 = zr(jcoor+3*(no2-1)+0)
    y2 = zr(jcoor+3*(no2-1)+1)
    z2 = zr(jcoor+3*(no2-1)+2)
    vx2 = x2-x0
    vy2 = y2-y0
    vz2 = z2-z0
!
    vnx = vy1*vz2-vy2*vz1
    vny = vz1*vx2-vz2*vx1
    vnz = vx1*vy2-vx2*vy1
!     SI I3 EST NUL NOUS SOMMES EN PRESENCE D'UN QUADRILATERE
    if (i3 .eq. 0) then
        vx3 = 0
        vy3 = 0
        vz3 = -1
    else
        no3 = char8_to_int(macoc(2+i3))
        x3 = zr(jcoor+3*(no3-1)+0)
        y3 = zr(jcoor+3*(no3-1)+1)
        z3 = zr(jcoor+3*(no3-1)+2)
        vx3 = x3-x0
        vy3 = y3-y0
        vz3 = z3-z0
    end if
!
!     VERIFICATION DU SIGNE DU JACOBIEN
!
    scal = vnx*vx3+vny*vy3+vnz*vz3
    if (niv .eq. 2) then
        write (ifm, *) 'SOMMET ', macoc(2+i0), ' DE LA MAILLE ', macoc(1), &
            ' VERS LES SOMMETS :'
        if (i3 .eq. 0) then
            write (ifm, *) macoc(2+i1), macoc(2+i2)
        else
            write (ifm, *) macoc(2+i1), macoc(2+i2), macoc(2+i3)
        end if
        write (ifm, *) 'JACOBIEN ', scal
    end if
    if (scal .lt. 0) then
        call utmess('E', 'MODELISA4_35')
    end if
!
    call jedema()
!
end subroutine
