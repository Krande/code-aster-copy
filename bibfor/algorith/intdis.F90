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
subroutine intdis(coint, nnoint, noddli, ddlsst, nbsst)
    implicit none
!    M. CORUS     DATE 05/02/10
!-----------------------------------------------------------------------
!
!  BUT:      < DETERMINATION DES PARTIES D'INTERFACES DISJOINTES >
!
!-----------------------------------------------------------------------
!  IN  : COINT  : DEFINITION DE LA CONNECTIVITE DE L'INTERFACE
!  IN  : NNOINT  : NOMBRE DE NOEUD A L'INTERFACE
!  IN  : NODDLI : DEFINITION DES DDL PORTES PAR LES NOEUDS D'INTERFACE
!  OUT : DDLSST   : DEFINITION DES DDL POUR CHAQUE PARTIE D'INTERFACE
!  OUT : NBSST    : NOMBRE DE PARTIE D'INTERFACE DISJOINTES
!-----------------------------------------------------------------------
!
!
!
!
!
!     ------------------------------------------------------------------
!
!-- VARIABLES EN ENTREES / SORTIE
#include "jeveux.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nnoint, nbsst
    character(len=24) :: coint, noddli, ddlsst
!
!-- VARIABLES DE LA ROUTINE
    integer(kind=8) :: i1, j1, k1, n1, l1, lconnc
    integer(kind=8) :: nz0, nz1, lindin, decal, nbno, nbvois, no, lnddli
    real(kind=8), pointer :: defi_ss_lib(:) => null()
    integer(kind=8), pointer :: numero_noeuds(:) => null()
    integer(kind=8), pointer :: vect_ind_mat(:) => null()
    integer(kind=8), pointer :: vect_indsst(:) => null()
    real(kind=8), pointer :: vect_temp(:) => null()
    integer(kind=8), pointer :: ind_noeud(:) => null()
!
!-----------C
!--       --C
!-- DEBUT --C
!--       --C
!-----------C
!
    call jemarq()
!
!-- CONSTRUCTION DE LA CONNECTIVITE REDUITE
!
    call jeveuo('&&MOIN93.IND_NOEUD', 'L', vi=ind_noeud)
    AS_ALLOCATE(vr=defi_ss_lib, size=nnoint**2)
    call jeveuo(coint, 'L', lconnc)
!
    do i1 = 1, nnoint
        nbvois = zi(lconnc+i1-1)
        do k1 = 1, nbvois
            no = zi(lconnc+nnoint*k1+i1-1)
            j1 = ind_noeud(no)
            defi_ss_lib(1+(j1-1)*nnoint+i1-1) = 1.d0
            defi_ss_lib(1+(j1-1)*nnoint+j1-1) = 1.d0
            defi_ss_lib(1+(i1-1)*nnoint+i1-1) = 1.d0
            defi_ss_lib(1+(i1-1)*nnoint+j1-1) = 1.d0
        end do
    end do
!
    AS_ALLOCATE(vi=numero_noeuds, size=nnoint)
!
    do i1 = 1, nnoint
        numero_noeuds(i1) = i1
    end do
!
    AS_ALLOCATE(vr=vect_temp, size=nnoint)
    AS_ALLOCATE(vi=vect_ind_mat, size=nnoint)
    AS_ALLOCATE(vi=vect_indsst, size=nnoint)
!
!-- INITIALISATION
!
    decal = 0
    nbsst = 0
    nbno = 0
    vect_indsst(1) = 1
!
!-- RECHERCHE DES PARTIES DISJOINTES
!
!      DO WHILE (NBNO .LT. NNOINT)
666 continue
    nz0 = 0
    k1 = 1
!        DO WHILE (NZ0 .EQ. 0)
667 continue
    if (numero_noeuds(k1) .gt. 0) then
        nz0 = 1
        vect_ind_mat(decal+1) = k1
    end if
    k1 = k1+1
    if (nz0 .eq. 0) then
        goto 667
    end if
!        END DO
!
    nz1 = 1
!        DO WHILE (NZ1 .GT. NZ0)
668 continue
    nz0 = nz1
    do j1 = 1, nz1
        do i1 = 1, nnoint
            vect_temp(i1) = vect_temp(i1)+defi_ss_lib(1+(vect_ind_mat(1+decal+j1- &
                                                                      1)-1)*nnoint+i1-1)
        end do
    end do
!
    nz1 = 0
    do i1 = 1, nnoint
        if (vect_temp(i1) .gt. 0.d0) then
            nz1 = nz1+1
            vect_ind_mat(1+decal+nz1-1) = i1
            numero_noeuds(i1) = 0
            vect_temp(i1) = 0.d0
        end if
    end do
!
    if (nz1 .gt. nz0) then
        goto 668
    end if
!        END DO
!
    nbsst = nbsst+1
    decal = decal+nz1
    nbno = nbno+nz1
    vect_indsst(1+2*nbsst-1) = decal
    vect_indsst(1+2*nbsst) = nbno+1
!
    if (nbno .lt. nnoint) then
        goto 666
    end if
!      END DO
!
    call jeveuo(noddli, 'L', lnddli)
    call wkvect(ddlsst, 'V V I', nbsst*6*nnoint, lindin)
    do i1 = 1, nbsst
        k1 = vect_indsst(1+2*(i1-1))
        l1 = vect_indsst(1+2*(i1-1)+1)
!
        do j1 = k1, l1
            do n1 = 1, 6
                zi(lindin+6*nnoint*(i1-1)+6*(vect_ind_mat(j1)-1)+n1-1) &
                    = 1
            end do
        end do
    end do
!
!----------------------------------------C
!--                                    --C
!-- DESTRUCTION DES OBJETS TEMPORAIRES --C
!--                                    --C
!----------------------------------------C
!
    AS_DEALLOCATE(vr=defi_ss_lib)
    AS_DEALLOCATE(vi=numero_noeuds)
    AS_DEALLOCATE(vr=vect_temp)
    AS_DEALLOCATE(vi=vect_ind_mat)
    AS_DEALLOCATE(vi=vect_indsst)
!
!---------C
!--     --C
!-- FIN --C
!--     --C
!---------C
!
    call jedema()
end subroutine
