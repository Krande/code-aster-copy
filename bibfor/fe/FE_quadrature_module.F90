! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elraga.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfno.h"
#include "asterfort/elrfvf.h"
#include "asterfort/fe_module.h"
#include "asterfort/lteatt.h"
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
        integer                             :: order = 0
        integer                             :: nbQuadPoints = 0
        real(kind=8), dimension(3, MAX_QP)  :: points_param = 0.d0
        real(kind=8), dimension(MAX_QP)     :: weights_param = 0.d0
        real(kind=8), dimension(3, MAX_QP)  :: points = 0.d0
        real(kind=8), dimension(MAX_QP)     :: weights = 0.d0
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
    private  :: FE_transfo, FE_rules, FE_order_rule, &
                FEQuadPrint, &
                FEinitCellFamiQ, FEinitFaceFamiQ, check_order, FESelectOrder

!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine check_order(order, maxAutorized)
!
        implicit none
!
        integer, intent(in) :: order
        integer, intent(in) :: maxAutorized
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   control the order of integration of order has to be included in [0, maxAutorized]
!   In order  : order of integration
!   In maxAutorized : maximum autorized order
!
! --------------------------------------------------------------------------------------------------
!
        if (order > maxAutorized) then
            ASSERT(ASTER_FALSE)
        end if
!
        if (order < 0) then
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!===================================================================================================
!
!===================================================================================================
!
    function FE_order_rule(typema, order) result(rule)
!
        implicit none
!
        integer, intent(in) :: order
        character(len=8), intent(in) :: typema
        character(len=8) :: rule
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   control the order of integration of order has to be included in [0, maxAutorized]
!   In order  : order of integration
!   In maxAutorized : maximum autorized order
!
! --------------------------------------------------------------------------------------------------
!
        integer :: max_order
        character(len=8), dimension(0:7) ::rules

        if (typema(1:2) == "SE") then
            rules = (/'FPG1', 'FPG1', 'FPG2', 'FPG2', 'FPG3', 'FPG3', 'FPG4', 'FPG4'/)
            max_order = 7
        elseif (typema(1:2) == "QU") then
            rules = (/'FPG1 ', 'FPG1 ', 'FPG4 ', 'FPG4 ', 'FPG9 ', 'FPG9 ', 'FPG16', 'FPG16'/)
            max_order = 7
        elseif (typema(1:2) == "TR") then
            rules = (/'FPG1 ', 'FPG1 ', 'FPG3 ', 'FPG4 ', 'FPG6 ', 'FPG7 ', 'FPG12', 'FPG13'/)
            max_order = 7
        elseif (typema(1:3) == "HE8" .or. typema(1:3) == "H20" .or. typema(1:3) == "H27") then
            rules = (/'FPG1 ', 'FPG1 ', 'FPG8 ', 'FPG8 ', 'FPG27', 'FPG27', 'FPG64', 'FPG64'/)
            max_order = 7
        elseif (typema(1:3) == "TE4" .or. typema(1:3) == "T10" .or. typema(1:3) == "T15") then
            rules = (/'FPG1 ', 'FPG1 ', 'FPG4 ', 'FPG5 ', 'FPG11', 'FPG15', 'FPG23', 'XXXXX'/)
            max_order = 6
        elseif (typema(1:3) == "PY5" .or. typema(1:3) == "P13" .or. typema(1:3) == "T19") then
            rules = (/'FPG1B', 'FPG1B', 'FPG5 ', 'FPG6 ', 'FPG10', 'FPG15', 'FPG24', 'FPG31'/)
            max_order = 7
        elseif (typema(1:3) == "PR6" .or. typema(1:3) == "P15" .or. &
                typema(1:3) == "P18" .or. typema(1:3) == "P21") then
            rules = (/'FPG1 ', 'FPG6B', 'FPG6B', 'FPG8 ', 'FPG21', 'FPG21', 'FPG29', 'XXXXX'/)
            max_order = 6
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call check_order(order, max_order)
!
        rule = rules(order)
!
    end function

!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FE_transfo(coorno, nbnodes, typema, coorref, coorac, jacob)
!
        implicit none
!
        integer, intent(in)                             :: nbnodes
        real(kind=8), dimension(3, nbnodes), intent(in) :: coorno
        character(len=8), intent(in)                    :: typema
        real(kind=8), dimension(3), intent(in)          :: coorref
        real(kind=8), dimension(3), intent(out)         :: coorac
        real(kind=8), intent(out)                       :: jacob
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
        real(kind=8), dimension(27) :: basis
        real(kind=8), dimension(3, 27) :: dbasis
        real(kind=8), dimension(3, 3) :: jaco
        integer :: i, ndim
!
! ----- shape function
!
        call elrfno(typema, ndim=ndim)
        call elrfvf(typema, coorref, basis)
!
! ----- derivative of shape function
!
        call elrfdf(typema, coorref, dbasis)
!
        coorac = 0.d0
!
        do i = 1, nbnodes
            coorac(1:3) = coorac(1:3)+coorno(1:3, i)*basis(i)
        end do
!
! ---  Compute the jacobienne
        jaco = 0.d0
        do i = 1, nbnodes
            jaco(1:3, 1) = jaco(1:3, 1)+coorno(1, i)*dbasis(1:3, i)
            jaco(1:3, 2) = jaco(1:3, 2)+coorno(2, i)*dbasis(1:3, i)
            jaco(1:3, 3) = jaco(1:3, 3)+coorno(3, i)*dbasis(1:3, i)
        end do
!
        if (ndim == 3) then
            jacob = jaco(1, 1)*jaco(2, 2)*jaco(3, 3)+jaco(1, 3)*jaco(2, 1)*jaco(3, 2) &
                    +jaco(3, 1)*jaco(1, 2)*jaco(2, 3)-jaco(3, 1)*jaco(2, 2)*jaco(1, 3) &
                    -jaco(3, 3)*jaco(2, 1)*jaco(1, 2)-jaco(1, 1)*jaco(2, 3)*jaco(3, 2)
        elseif (ndim == 2) then
            jacob = jaco(1, 1)*jaco(2, 2)-jaco(2, 1)*jaco(1, 2)
        elseif (ndim == 1) then
            jacob = jaco(1, 1)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
! --------------------------------------------------------------------------------------------------
!       !!!! BE CAREFULL ALL THE TRANSFORMATIONS ARE LINEAR (PLANAR ELEMENTS) !!!!
! --------------------------------------------------------------------------------------------------
    subroutine FE_rules(this, coorno, nbnodes, typema)
!
        implicit none
!
        real(kind=8), dimension(3, *), intent(in)  :: coorno
        integer, intent(in)                        :: nbnodes
        character(len=8), intent(in)               :: typema
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
        character(len=8) :: rule
        integer :: dimp, nbpg, ipg
        real(kind=8), dimension(MAX_QP) :: poidpg
        real(kind=8), dimension(3*MAX_QP) :: coorpg
        real(kind=8) :: x, y, z, coorac(3), jaco
!
!------ get quadrature points
        rule = FE_order_rule(typema, this%order)
        coorpg = 0.d0
        poidpg = 0.d0
        call elraga(typema, rule, dimp, nbpg, coorpg, poidpg)
!
! ----- fill FEQuad
        ASSERT(nbpg <= MAX_QP)
        this%nbQuadPoints = nbpg
!
        do ipg = 1, nbpg
            x = coorpg(dimp*(ipg-1)+1)
            y = coorpg(dimp*(ipg-1)+2)
            if (dimp == 3) then
                z = coorpg(dimp*(ipg-1)+3)
            else
                z = 0.d0
            end if
            call FE_transfo(coorno, nbnodes, typema, (/x, y, z/), coorac, jaco)
            this%points_param(1:3, ipg) = (/x, y, z/)
            this%weights_param(ipg) = poidpg(ipg)
!
            this%points(1:3, ipg) = coorac
            this%weights(ipg) = abs(jaco)*poidpg(ipg)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FESelectOrder(typema, npg, order)
!
        implicit none
!
        character(len=8), intent(in)  :: typema
        integer, intent(in)           :: npg
        integer, intent(out)          :: order
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Get the order of the quadrature rules from a familly definied in the catalogue of code_aster
!   In typema   : type of Cell or Face
!   In npg      : number of quadrature points
!   Out order   : order of the quadrature
!
! --------------------------------------------------------------------------------------------------
!
        order = 0
!
        if (typema(1:4) == 'HEXA') then
            select case (npg)
            case (1)
                order = 1
            case (8)
                order = 3
            case (27)
                order = 5
            case (64)
                order = 7
            case default
                ASSERT(ASTER_FALSE)
            end select
        elseif (typema(1:5) == 'PENTA') then
            select case (npg)
            case (1)
                order = 0
            case (6)
                order = 2
            case (8)
                order = 3
            case (21)
                order = 5
            case (29)
                order = 6
            case default
                ASSERT(ASTER_FALSE)
            end select
        elseif (typema(1:5) == 'PYRAM') then
            select case (npg)
            case (1)
                order = 1
            case (5)
                order = 2
            case (6)
                order = 3
            case (10)
                order = 4
            case (15)
                order = 5
            case (24)
                order = 6
            case (31)
                order = 7
            case default
                ASSERT(ASTER_FALSE)
            end select
        elseif (typema(1:5) == 'TETRA') then
            select case (npg)
            case (1)
                order = 1
            case (4)
                order = 2
            case (5)
                order = 3
            case (11)
                order = 4
            case (15)
                order = 5
            case (23)
                order = 6
            case default
                ASSERT(ASTER_FALSE)
            end select
        elseif (typema(1:4) == 'QUAD') then
            select case (npg)
            case (1)
                order = 1
            case (4)
                order = 3
            case (9)
                order = 5
            case (16)
                order = 7
            case default
                ASSERT(ASTER_FALSE)
            end select
        elseif (typema(1:4) == 'TRIA') then
            select case (npg)
            case (1)
                order = 1
            case (3)
                order = 2
            case (4)
                order = 3
            case (6)
                order = 4
            case (7)
                order = 5
            case (12)
                order = 6
            case (13)
                order = 7
            case (16)
                order = 8
            case default
                ASSERT(ASTER_FALSE)
            end select
        elseif (typema(1:3) == 'SEG') then
            select case (npg)
            case (1)
                order = 1
            case (2)
                order = 3
            case (3)
                order = 5
            case (4)
                order = 7
            case default
                ASSERT(ASTER_FALSE)
            end select
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine FEinitCellFamiQ(this, FECell, order, fami)
!
        implicit none
!
        type(FE_cell), intent(in)          :: FECell
        integer, intent(in), optional      :: order
        character(len=*), intent(in), optional :: fami
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
        integer :: npg, ipg
!
        if (present(fami)) then
            call elrefe_info(fami=fami, npg=npg)
!
            ASSERT(npg .le. MAX_QP)
            call FESelectOrder(FECell%typema, npg, this%order)
        else
            ASSERT(present(order))
            this%order = order
        end if
!
        call this%FE_rules(FECell%coorno(1:3, 1:FECell%nbnodes), FECell%nbnodes, FECell%typemas)
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
    subroutine FEinitFaceFamiQ(this, FEFace, order, fami)
!
        implicit none
!
        type(FE_Face), intent(in)          :: FEFace
        integer, intent(in), optional      :: order
        character(len=*), intent(in), optional :: fami
        class(FE_quadrature), intent(out)  :: this
!
! --------------------------------------------------------------------------------------------------
!   FE
!
!   Get the quadrature rules from a familly definied in the catalogue of code_aster
!   In FEFace  : FEFace
!   In npg      : number of quadrature points
!   In axis     : axisymetric ? multpiply by r the weith if True
!   Out this    : FE quadrature
!
! --------------------------------------------------------------------------------------------------
!
        integer :: npg, ipg
!
        if (present(fami)) then
            call elrefe_info(fami=fami, npg=npg)
!
            ASSERT(npg .le. MAX_QP)
            call FESelectOrder(FEFace%typema, npg, this%order)
        else
            ASSERT(present(order))
            this%order = order
        end if
!
        call this%FE_rules(FEFace%coorno(1:3, 1:FEFace%nbnodes), FEFace%nbnodes, FEFace%typemas)
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
        integer :: ipg
!
        write (6, *) "QUADRATURE INFORMATIONS"
        write (6, *) "number of qp: ", this%nbQuadPoints
        write (6, *) "order: ", this%order
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
