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

function typmat(nbmat, tlimat)
    implicit none
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeexin.h"
#include "asterfort/redetr.h"
#include "asterfort/isParallelMesh.h"
    character(len=*) :: tlimat(*)
    integer(kind=8) :: nbmat
    character(len=1) :: typmat
!-----------------------------------------------------------------------
!
!  BUT : CETTE FONCTION RETOURNE LE TYPE SYMETRIQUE : 'S'
!                                     OU PLEINE     : 'N'
!        DE LA MATRICE GLOBALE RESULTANTE DE L'ASSEMBLAGE
!        DES MATR_ELEM TLIMAT
!  ATTENTION :
!    SI ON EST EN PARALLELE MPI, LA REPONSE EST GLOBALE.
!    CETTE ROUTINE FAIT DU MPI_ALLREDUCE. IL FAUT DONC L'APPELER
!    DE LA MEME FACON SUR TOUS LES PROCS.
!
!-----------------------------------------------------------------------
! --- DESCRIPTION DES PARAMETRES
! IN  I  NBMAT  : NOMBRE DE MATR_ELEM DE LA LISTE TLIMAT
! IN  K* TLIMAT : LISTE DES MATR_ELEM
! ----------------------------------------------------------------------
!----------------------------------------------------------------------
    character(len=8) :: sym, zero, mesh
    character(len=19) :: matel
    integer(kind=8) :: i, itymat
    integer(kind=8) :: iexi
    aster_logical :: l_pmesh
!----------------------------------------------------------------------
!     ITYMAT =  0 -> SYMETRIQUE
!            =  1 -> NON-SYMETRIQUE
!
!     --  PAR DEFAUT LE TYPE DE MATRICE EST SYMETRIQUE
    itymat = 0
    l_pmesh = ASTER_FALSE
!
    do i = 1, nbmat
        matel = tlimat(i)
        call jeexin(matel//'.RELR', iexi)
        iexi = min(1, abs(iexi))
        if (iexi .ne. 0) then
!
!       -- LA LOGIQUE CI-DESSOUS N'EST VALABLE QUE SI LE MATR_ELEM
!          A ETE EXPURGE DE SES RESUELEM NULS => CALL REDETR()
            call redetr(matel)
!
            call dismoi('TYPE_MATRICE', matel, 'MATR_ELEM', repk=sym)
            if (sym .eq. 'NON_SYM') then
                call dismoi('ZERO', matel, 'MATR_ELEM', repk=zero)
                if (zero .eq. 'NON') then
                    itymat = 1
                end if
            end if

            call dismoi('NOM_MAILLA', matel, 'MATR_ELEM', repk=mesh)
            l_pmesh = isParallelMesh(mesh)
!
            if (.not. l_pmesh) then
!
! --- Il faut communiquer entre proc pour sortir tous en même temps
                call asmpi_comm_vect('MPI_MAX', 'I', sci=itymat)
            end if
!
            if (itymat .eq. 1) then
                exit
            end if
        end if
!
    end do
!
    if (l_pmesh) then
        call asmpi_comm_vect('MPI_MAX', 'I', sci=itymat)
    end if
!
    if (itymat .eq. 0) then
        typmat = 'S'
    else
        typmat = 'N'
    end if
end function
