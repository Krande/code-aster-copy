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

subroutine irgmpv(ifi, lresu, nomcon, chamsy, nbordr, &
                  para, nocmp, nbel, scal, vect, &
                  tens, versio)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterfort/jenonu.h"
#include "asterfort/jexnom.h"
#include "asterfort/lxlgut.h"
    integer(kind=8) :: ifi, nbordr, lch, ich, versio
    real(kind=8) :: para(*)
    aster_logical :: lresu, scal, vect, tens
    character(len=8) :: nocmp
    character(len=*) :: nomcon, chamsy
!     NBRE POUR CHAQUE TYPE D'ELEMENT
    integer(kind=8) :: nbel(*)
!
!     BUT :   ECRITURE D'UN RESULTAT AU FORMAT GMSH
!
!     ------------------------------------------------------------------
    integer(kind=8) :: ior
    character(len=8) :: nomsd
    character(len=50) :: k50b
    integer(kind=8) :: nbpoi, nbseg, nbtri, nbtet, nbqua, nbpyr, nbpri, nbhex
    integer(kind=8) :: typpoi, typseg, typtri, typtet, typqua
    integer(kind=8) :: typpyr, typpri, typhex
!     ------------------------------------------------------------------
!
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), typpoi)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), typseg)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA3'), typtri)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD4'), typqua)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA4'), typtet)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM5'), typpyr)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA6'), typpri)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA8'), typhex)
    nbpoi = nbel(typpoi)
    nbseg = nbel(typseg)
    nbtri = nbel(typtri)
    nbqua = nbel(typqua)
    nbtet = nbel(typtet)
    nbpyr = nbel(typpyr)
    nbpri = nbel(typpri)
    nbhex = nbel(typhex)
!
    write (ifi, 1000) '$View'
!
!     ECRITURE DE LIGNE 1 (VIEW_NAME NB_TIME_STEPS)
!
    nomsd = nomcon(1:8)
!
    if (lresu) then
!
        lch = lxlgut(nomsd)
        k50b(1:lch) = nomsd(1:lch)
        ich = lch+1
!
        k50b(ich:ich) = '_'
        lch = lxlgut(chamsy)
        k50b(ich+1:ich+lch) = chamsy(1:lch)
        ich = ich+lch+1
!
        k50b(ich:ich) = '_'
        lch = lxlgut(nocmp)
        k50b(ich+1:ich+lch) = nocmp(1:lch)
!
        k50b(ich+lch+1:ich+lch+1) = ' '
        write (ifi, 1020) k50b(1:ich+lch+1), nbordr
!
    else
!
        lch = lxlgut(nomsd)
        k50b(1:lch) = nomsd(1:lch)
!
        ich = lch+1
        k50b(ich:ich) = '_'
        lch = lxlgut(nocmp)
        k50b(ich+1:ich+lch) = nocmp(1:lch)
!
        k50b(ich+lch+1:ich+lch+1) = ' '
        write (ifi, 1022) k50b(1:ich+lch+1), nbordr
!
    end if
!
!     ECRITURE DE LA LIGNE 2 A 4 (nb elt par de type de maille)
!
    if (scal) then
        write (ifi, 1030) nbpoi, 0, 0
        write (ifi, 1030) nbseg, 0, 0
        write (ifi, 1030) nbtri, 0, 0
        if (versio .eq. 2) then
            write (ifi, 1030) nbqua, 0, 0
        end if
        write (ifi, 1030) nbtet, 0, 0
        if (versio .eq. 2) then
            write (ifi, 1030) nbhex, 0, 0
            write (ifi, 1030) nbpri, 0, 0
            write (ifi, 1030) nbpyr, 0, 0
        end if
    else if (vect) then
        write (ifi, 1030) 0, nbpoi, 0
        write (ifi, 1030) 0, nbseg, 0
        write (ifi, 1030) 0, nbtri, 0
        if (versio .eq. 2) then
            write (ifi, 1030) 0, nbqua, 0
        end if
        write (ifi, 1030) 0, nbtet, 0
        if (versio .eq. 2) then
            write (ifi, 1030) 0, nbhex, 0
            write (ifi, 1030) 0, nbpri, 0
            write (ifi, 1030) 0, nbpyr, 0
        end if
    else if (tens) then
        write (ifi, 1030) 0, 0, nbpoi
        write (ifi, 1030) 0, 0, nbseg
        write (ifi, 1030) 0, 0, nbtri
        if (versio .eq. 2) then
            write (ifi, 1030) 0, 0, nbqua
        end if
        write (ifi, 1030) 0, 0, nbtet
        if (versio .eq. 2) then
            write (ifi, 1030) 0, 0, nbhex
            write (ifi, 1030) 0, 0, nbpri
            write (ifi, 1030) 0, 0, nbpyr
        end if
    else
    end if
!
    if (versio .eq. 2) then
        write (ifi, 1050) 0, 0, 0, 0
    end if
!
!     ECRITURE DE LA LIGNE 5 (time_step_values)
!
    write (ifi, 1040) (para(ior), ior=1, nbordr)
!
1000 format(a5)
1020 format(a, 1x, i4)
1022 format(a, 1x, i4)
1030 format(3(i8, 1x))
1040 format(1p, 10(e15.8, 1x))
1050 format(4(i6, 1x))
!
end subroutine
