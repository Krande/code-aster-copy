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
subroutine cargeo(meshZ)
!
    use mesh_module
    implicit none
!
#include "jeveux.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/detrsd.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ltcrsd.h"
#include "asterfort/ltnotb.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
!
    character(len=*), intent(in) :: meshZ
!
! --------------------------------------------------------------------------------------------------
!
! CALCULER DES CARACTERISTIQUES DU MAILLAGE
!
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbNode, jvale, jdime, iNode, nbParaTable, ibid
    integer(kind=8) :: iret
    real(kind=8) :: xmax, ymax, zmax, xmin, ymin, zmin
    real(kind=8) :: armin, armax
    complex(kind=8) :: c16b
    character(len=8) :: k8b, mesh
    character(len=1) :: bas1
    character(len=19) :: nomt19
    character(len=24) :: nodime, connex, coordo
    aster_logical :: l_pmesh

    integer(kind=8), parameter :: nbPara = 8
    character(len=8), parameter :: paraName(nbPara) = (/'X_MIN ', 'X_MAX ', &
                                                        'Y_MIN ', 'Y_MAX ', &
                                                        'Z_MIN ', 'Z_MAX ', &
                                                        'AR_MIN', 'AR_MAX'/)
    character(len=8), parameter :: paraType(nbPara) = (/'R', 'R', &
                                                        'R', 'R', &
                                                        'R', 'R', &
                                                        'R', 'R'/)
    real(kind=8) :: paraVale(nbPara)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    c16b = (0.d0, 0.d0)
    ibid = 0
!
    mesh = meshZ
    nodime = mesh//'.DIME           '
    connex = mesh//'.CONNEX         '
    coordo = mesh//'.COORDO    .VALE'
    l_pmesh = isParallelMesh(mesh)
!
    call jelira(nodime, 'CLAS', cval=bas1)
    ASSERT(bas1 .eq. 'G' .or. bas1 .eq. 'V')
!
    call jeveuo(nodime, 'L', jdime)
    nbNode = zi(jdime)

! - Global coordinates of mesh
    call jeveuo(coordo, 'L', jvale)
    xmax = zr(jvale)
    ymax = zr(jvale+1)
    zmax = zr(jvale+2)
    xmin = xmax
    ymin = ymax
    zmin = zmax
    do iNode = 2, nbNode
        xmax = max(xmax, zr(jvale+3*(iNode-1)))
        xmin = min(xmin, zr(jvale+3*(iNode-1)))
        ymax = max(ymax, zr(jvale+3*(iNode-1)+1))
        ymin = min(ymin, zr(jvale+3*(iNode-1)+1))
        zmax = max(zmax, zr(jvale+3*(iNode-1)+2))
        zmin = min(zmin, zr(jvale+3*(iNode-1)+2))
    end do
    if (l_pmesh) then
        call asmpi_comm_vect("MPI_MIN", 'R', scr=xmin)
        call asmpi_comm_vect("MPI_MIN", 'R', scr=ymin)
        call asmpi_comm_vect("MPI_MIN", 'R', scr=zmin)
        call asmpi_comm_vect("MPI_MAX", 'R', scr=xmax)
        call asmpi_comm_vect("MPI_MAX", 'R', scr=ymax)
        call asmpi_comm_vect("MPI_MAX", 'R', scr=zmax)
    end if
!
    paraVale(1) = xmin
    paraVale(2) = xmax
    paraVale(3) = ymin
    paraVale(4) = ymax
    paraVale(5) = zmin
    paraVale(6) = zmax
    paraVale(7) = 0
    paraVale(8) = 0
!
    call jeexin(connex, iret)
    if (iret .eq. 0) then
        nbParaTable = 6
        goto 100
    else
        nbParaTable = nbPara
    end if

    call compMinMaxEdges(mesh, armin, armax)
    paraVale(7) = armin
    paraVale(8) = armax
!
100 continue
! - Write in table
    nomt19 = ' '
    call jeexin(mesh//'           .LTNT', iret)
    if (iret .ne. 0) then
        call ltnotb(meshZ, 'CARA_GEOM', nomt19)
        call detrsd('TABLE', nomt19)
    else
        call ltcrsd(meshZ, bas1)
    end if
    call ltnotb(meshZ, 'CARA_GEOM', nomt19)
!
    call jeexin(nomt19//'.TBBA', iret)
    if (iret .ne. 0) call detrsd('TABLE', nomt19)
!
    call tbcrsd(nomt19, bas1)
    call tbajpa(nomt19, nbParaTable, paraName, paraType)
    call tbajli(nomt19, nbParaTable, paraName, [ibid], paraVale, &
                [c16b], k8b, 0)
!
    call jedema()
end subroutine
