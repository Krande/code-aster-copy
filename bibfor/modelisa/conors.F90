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
subroutine conors(i1, i2, i3, macoc, nbcoc, &
                  macor, nbcor, loreor, mailla)
    implicit none
!
! DESCRIPTION : REORIENTATION D'UNE MAILLE DE FISSURE
! -----------   LA NORMALE A CHAQUE BORD EST SORTANTE
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
! OUT    : LOREOR : LOGICAL , SCALAIRE
!                   INDICATEUR DE REORIENTATION DE LA MAILLE DE FISSURE
!                   SI LOREOR =.TRUE. LA CONNECTIVITE DOIT ETRE MODIFIEE
!                   AFIN D'OBTENIR UNE NORMALE SORTANTE
!
!-------------------   DECLARATION DES VARIABLES   ---------------------
!
!
! ARGUMENTS
! ---------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
!
    integer(kind=8) :: i1, i2, i3
    character(len=8) :: mailla
    integer(kind=8) :: nbcoc
    character(len=8) :: macoc(2+nbcoc)
    aster_logical :: loreor
    integer(kind=8) :: nbcor
    character(len=8) :: macor(2+nbcor)
!
! VARIABLES LOCALES
! -----------------
    integer(kind=8) :: jcoor, no, no1, no2, no3, inor
    character(len=24) :: coorno
    real(kind=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, scal
    real(kind=8) :: xg, yg, zg, vnx, vny, vnz, vrx, vry, vrz
    real(kind=8) :: vx1, vy1, vz1, vx2, vy2, vz2
!
! FONCTIONS EXTERNES
! ------------------
!     EXTERNAL      JEXNOM, JEXNUM, R8DOT, R8NRM2
!
! ROUTINES EXTERNES
! -----------------
!     EXTERNAL      JEDEMA, JEMARQ, JENONU, JEVEUO,
!                   DISMOI
!
!-------------------   DEBUT DU CODE EXECUTABLE    ---------------------
!
    call jemarq()
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
    no1 = char8_to_int(macoc(2+i1))
    x1 = zr(jcoor+3*(no1-1)+0)
    y1 = zr(jcoor+3*(no1-1)+1)
    z1 = zr(jcoor+3*(no1-1)+2)
!
    no2 = char8_to_int(macoc(2+i2))
    x2 = zr(jcoor+3*(no2-1)+0)
    y2 = zr(jcoor+3*(no2-1)+1)
    z2 = zr(jcoor+3*(no2-1)+2)
    vx1 = x2-x1
    vy1 = y2-y1
    vz1 = z2-z1
!     SI I3 EST NUL NOUS SOMMES EN PRESENCE D'UN QUADRILATERE
    if (i3 .eq. 0) then
        vnx = -vy1
        vny = vx1
        vnz = 0
    else
        no3 = char8_to_int(macoc(2+i3))
        x3 = zr(jcoor+3*(no3-1)+0)
        y3 = zr(jcoor+3*(no3-1)+1)
        z3 = zr(jcoor+3*(no3-1)+2)
        vx2 = x3-x1
        vy2 = y3-y1
        vz2 = z3-z1
        vnx = vy1*vz2-vy2*vz1
        vny = vz1*vx2-vz2*vx1
        vnz = vx1*vy2-vx2*vy1
    end if
!
!     CENTRE DE GRAVITE DE LA MAILLE DE REFERENCE
!
    xg = 0.d0
    yg = 0.d0
    zg = 0.d0
    do inor = 1, nbcor
        no = char8_to_int(macor(2+inor))
        xg = xg+zr(jcoor+3*(no-1)+0)
        yg = yg+zr(jcoor+3*(no-1)+1)
        zg = zg+zr(jcoor+3*(no-1)+2)
    end do
    xg = xg/nbcor
    yg = yg/nbcor
    zg = zg/nbcor
!
!     VECTEUR NOEUD 1 - CENTRE DE GRAVITE
!
    vrx = xg-x1
    vry = yg-y1
    vrz = zg-z1
!
!     VERIFICATION QUE LA NORMALE EST SORTANTE
!
    scal = vnx*vrx+vny*vry+vnz*vrz
!
    if (scal .eq. 0) then
        call utmess('E', 'MODELISA4_36', sk=macoc(1))
    end if
!
    loreor = scal .gt. 0
!
    call jedema()
!
end subroutine
