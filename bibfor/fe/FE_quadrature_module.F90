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
! person_in_charge: mickael.abbas at edf.fr
!
module FE_quadrature_module
!
    use FE_topo_module
!
    implicit none
!
    private
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! FE - generic
!
! Module to generate quadratures used for FE
!
! --------------------------------------------------------------------------------------------------
!
    type FE_Quadrature
        character(len=8)                    :: fami = " "
        integer(kind=8)                             :: nbQuadPoints = 0
        real(kind=8), dimension(3, MAX_QP)  :: points_param = 0.d0
        real(kind=8), dimension(MAX_QP)     :: weights_param = 0.d0
        real(kind=8), dimension(3, MAX_QP)  :: points = 0.d0
        real(kind=8), dimension(MAX_QP)     :: weights = 0.d0
        real(kind=8), dimension(3, 3, MAX_QP) :: jacob = 0.d0
! ----- member functions
    contains
        procedure, private, pass :: FE_rules
        procedure, public, pass :: print => FEQuadPrint
        procedure, public, pass :: initCell => FEinitCellFamiQ
        procedure, public, pass :: initFace => FEinitFaceFamiQ
    end type
!
!===================================================================================================
!
!===================================================================================================
!
    public   :: FE_Quadrature
    private  :: FE_transfo, FE_rules, FEQuadPrint, &
                FEinitCellFamiQ, FEinitFaceFamiQ
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FE_transfo(coorno, nbnodes, ndim, l_skin, &
                          basis, dbasis, coorac, jaco, jacob)
!
        implicit none
!
        integer(kind=8), intent(in)                             :: nbnodes, ndim
        real(kind=8), dimension(3, nbnodes), intent(in) :: coorno
        aster_logical, intent(in)                       :: l_skin
        real(kind=8), intent(in)                        :: basis(*), dbasis(*)
        real(kind=8), dimension(3), intent(out)         :: coorac
        real(kind=8), intent(out)                       :: jacob
        real(kind=8), dimension(3, 3), intent(out)      :: jaco

!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   From reference element to current element
!   In coorno       : coordinates of the nodes
!   In coorref      : coordinates in the reference conf
!   Out coorac      : coordinates in the current conf
!   Out jacob       : determiant of the jacobienne of the transformation
!
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: normal(3)
        integer(kind=8) :: i, iadzi, iazk24, ind
!
        coorac = 0.d0
!
        do i = 1, nbnodes
            coorac(1:3) = coorac(1:3)+coorno(1:3, i)*basis(i)
        end do
!
! ---  Compute the jacobienne
        jaco = 0.d0
        if (ndim == 3) then
            do i = 1, nbnodes
                ind = 3*(i-1)
                jaco(1, 1:3) = jaco(1, 1:3)+coorno(1:3, i)*dbasis(ind+1)
                jaco(2, 1:3) = jaco(2, 1:3)+coorno(1:3, i)*dbasis(ind+2)
                jaco(3, 1:3) = jaco(3, 1:3)+coorno(1:3, i)*dbasis(ind+3)
            end do
!
            ASSERT(.not. l_skin)
            jacob = jaco(1, 1)*jaco(2, 2)*jaco(3, 3)+jaco(1, 3)*jaco(2, 1)*jaco(3, 2) &
                    +jaco(3, 1)*jaco(1, 2)*jaco(2, 3)-jaco(3, 1)*jaco(2, 2)*jaco(1, 3) &
                    -jaco(3, 3)*jaco(2, 1)*jaco(1, 2)-jaco(1, 1)*jaco(2, 3)*jaco(3, 2)
        elseif (ndim == 2) then
            do i = 1, nbnodes
                ind = 2*(i-1)
                jaco(1, 1:3) = jaco(1, 1:3)+coorno(1:3, i)*dbasis(ind+1)
                jaco(2, 1:3) = jaco(2, 1:3)+coorno(1:3, i)*dbasis(ind+2)
            end do
            if (l_skin) then
                normal(1) = jaco(1, 2)*jaco(2, 3)-jaco(1, 3)*jaco(2, 2)
                normal(2) = jaco(1, 3)*jaco(2, 1)-jaco(1, 1)*jaco(2, 3)
                normal(3) = jaco(1, 1)*jaco(2, 2)-jaco(1, 2)*jaco(2, 1)
                jacob = sqrt(normal(1)*normal(1)+normal(2)*normal(2)+normal(3)*normal(3))
            else
                jaco(3, 3) = 1.d0
                jacob = jaco(1, 1)*jaco(2, 2)-jaco(2, 1)*jaco(1, 2)
            end if
        elseif (ndim == 1) then
            do i = 1, nbnodes
                jaco(1, 1:2) = jaco(1, 1:2)+coorno(1:2, i)*dbasis(i)
            end do
            if (l_skin) then
                jacob = sqrt(jaco(1, 1)*jaco(1, 1)+jaco(1, 2)*jaco(1, 2))
            else
                jacob = jaco(1, 1)
                jaco(2, 2) = 1.d0
                jaco(3, 3) = 1.d0
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (abs(jacob) .le. 1.d0/r8gaem()) then
            call tecael(iadzi, iazk24)
            call utmess('F', 'ALGORITH2_59', sk=zk24(iazk24-1+3) (1:8))
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FE_rules(this, coorno, nbnodes, l_skin)
!
        implicit none
!
        real(kind=8), dimension(3, *), intent(in)  :: coorno
        integer(kind=8), intent(in)                        :: nbnodes
        aster_logical, intent(in)                  :: l_skin
        class(FE_quadrature), intent(inout)        :: this
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Get the quadrature rules for an edge
!   In coorno       : coordinates of the nodes
!   In measure      : length of the edge
!   In barycenter   : barycenter
!   Out this        : FE quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nbpg, ipg, jpoids, jcoopg, idim, dimp, jdfde, jvf, ind
        real(kind=8) :: coorac(3), jaco
!
!------ get quadrature points
        call elrefe_info(fami=this%fami, ndim=dimp, npg=nbpg, jpoids=jpoids, jcoopg=jcoopg, &
                         jvf=jvf, jdfde=jdfde)
!
! ----- fill FEQuad
        ASSERT(nbpg <= MAX_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            ind = jcoopg-1+dimp*(ipg-1)
            do idim = 1, dimp
                this%points_param(idim, ipg) = zr(ind+idim)
            end do
            this%weights_param(ipg) = zr(jpoids-1+ipg)
!
            call FE_transfo(coorno, nbnodes, dimp, l_skin, &
                            zr(jvf+nbnodes*(ipg-1)), zr(jdfde+dimp*nbnodes*(ipg-1)), &
                            coorac, this%jacob(1:3, 1:3, ipg), jaco)
!
            this%points(1:3, ipg) = coorac
            this%weights(ipg) = abs(jaco)*zr(jpoids-1+ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEinitCellFamiQ(this, FECell, fami)
!
        implicit none
!
        type(FE_cell), intent(in)          :: FECell
        character(len=*), intent(in)       :: fami
        class(FE_quadrature), intent(out)  :: this
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Get the quadrature rules from a familly definied in the catalogue of code_aster
!   In FECell  : FECell
!   In npg      : number of quadrature points
!   In axis     : axisymetric ? multpiply by r the weith if True
!   Out this    : FE quadrature
!
! --------------------------------------------------------------------------------------------------
        integer(kind=8) :: ipg
!
        this%fami = fami
!
        call this%FE_rules(FECell%coorno(1:3, 1:FECell%nbnodes), FECell%nbnodes, ASTER_FALSE)
!
        if (lteatt('AXIS', 'OUI')) then
            do ipg = 1, this%nbQuadPoints
                this%weights(ipg) = this%weights(ipg)*this%points(1, ipg)
            end do
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEinitFaceFamiQ(this, FESkin, fami)
!
        implicit none
!
        type(FE_Skin), intent(in)          :: FESkin
        character(len=*), intent(in)       :: fami
        class(FE_quadrature), intent(out)  :: this
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Get the quadrature rules from a familly definied in the catalogue of code_aster
!   In FESkin  : FESkin
!   In npg      : number of quadrature points
!   In axis     : axisymetric ? multpiply by r the weith if True
!   Out this    : FE quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ipg
!
        this%fami = fami
!
        call this%FE_rules(FESkin%coorno(1:3, 1:FESkin%nbnodes), FESkin%nbnodes, ASTER_TRUE)
!
        if (lteatt('AXIS', 'OUI')) then
            do ipg = 1, this%nbQuadPoints
                this%weights(ipg) = this%weights(ipg)*this%points(1, ipg)
            end do
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEQuadPrint(this)
!
        implicit none
!
        class(FE_quadrature), intent(in)  :: this
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Print quadrature informations
!   In this    : FE quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ipg
!
        write (6, *) "QUADRATURE INFORMATIONS"
        write (6, *) "familly: ", this%fami
        write (6, *) "number of qp: ", this%nbQuadPoints
!
        do ipg = 1, this%nbQuadPoints
            write (6, *) "coordo qp ", ipg, ":", this%points(1:3, ipg)
            write (6, *) "weight qp ", ipg, ":", this%weights(ipg)
        end do
        write (6, *) "END QUADRATURE INFORMATIONS"
!
    end subroutine
!
end module
