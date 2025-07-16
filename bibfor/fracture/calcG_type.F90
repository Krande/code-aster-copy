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
module calcG_type
!
    implicit none
!
    private
#include "asterc/getres.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcG_type.h"
#include "asterfort/ccbcop.h"
#include "asterfort/cgComporNodes.h"
#include "asterfort/cgcrio.h"
#include "asterfort/cgReadCompor.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/comp_info.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/fointe.h"
#include "asterfort/gcncon.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/ismali.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "jeveux.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/medomg.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsmena.h"
#include "asterfort/rsrusd.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbajvi.h"
#include "asterfort/tbajvk.h"
#include "asterfort/tbajvr.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xcourb.h"
#include "asterfort/char8_to_int.h"

!
    public :: CalcG_Field, CalcG_Study, CalcG_Theta, CalcG_Table, CalcG_Stat
!
! --------------------------------------------------------------------------------------------------
!
! CALC_G
!
! Define types for datastructures
!
! --------------------------------------------------------------------------------------------------
!
    type CalcG_Field
!
        aster_logical      :: l_debug = ASTER_FALSE
! ----- name of result (in)
        character(len=8)   :: result_in = ' '
! ----- name of temproray result (out)
        character(len=8)   :: result_out = 'resuout'
! ----- type of result (in)
        character(len=16)  :: result_in_type = ' '
! ----- name of table container (out)
        character(len=24)  :: table_out = ' '
! ----- CARTE for behavior
        character(len=24)  :: compor = ' '
! ----- Incremental behavior or not ?
        aster_logical      :: l_incr
!------ temperature in VARC
        aster_logical      :: l_temp = ASTER_FALSE
! ----- topological dimension
        integer(kind=8)            :: ndim = 0
! ----- name of list of NUME
        character(len=24)  :: list_nume_name = ' '
! ----- list of nume
        integer(kind=8), pointer   :: list_nume(:) => null()
! ----- number of order
        integer(kind=8)            :: nb_nume = 0
! ----- number option to compute
        integer(kind=8)            :: nb_option = 0
! ----- list of options
        character(len=8)   :: list_option(NB_MAX_OPT) = ' '
! ----- level information
        integer(kind=8)            :: level_info = 1
    contains
        procedure, pass    :: initialize => initialize_field
        procedure, pass    :: print => print_field
        procedure, pass    :: isModeMeca
        procedure, pass    :: isDynaTrans
        procedure, pass    :: clean => clean_field
!
    end type CalcG_Field
!
!===================================================================================================
!
!===================================================================================================
!
    type CalcG_Study
! ----- name of model
        character(len=8)   :: model = ' '
! ----- name of mesh
        character(len=8)   :: mesh = ' '
! ----- name of material
        character(len=24)  :: material = ' '
! ----- name of coded material
        character(len=24)  :: mateco = ' '
! ----- name of elementary characteristics
        character(len=24)  :: carael = ' '
! ----- name of loading
        character(len=24)  :: loading = ' '
! ----- index order
        integer(kind=8)            :: nume_ordre = -1
! ----- option to compute
        character(len=8)   :: option = ' '
!------ modal analysis ?
        aster_logical      :: l_modal = ASTER_FALSE
!------ axisymetric model
        aster_logical      :: l_axis = ASTER_FALSE
!------ inco model
        aster_logical      :: l_exi_inco = ASTER_FALSE
! ----- displacement field
        character(len=24)  :: depl = ' '
! ----- speed field
        character(len=24)  :: vitesse = ' '
! ----- acceleration field
        character(len=24)  :: acce = ' '
! ----- time
        real(kind=8)       :: time = 0.d0
! ----- pulse
        real(kind=8)       :: pulse = 0.d0
! ----- computed values (G, K1, K2, K3, FIC1, FIC2, FIC3)
        real(kind=8)       :: gth(NB_MAX_TERM) = 0.d0
! ----- member function
    contains
        procedure, pass    :: initialize => initialize_study
        procedure, pass    :: print => print_study
        procedure, pass    :: setOption
        procedure, pass    :: getField
        procedure, pass    :: getParameter
    end type CalcG_Study
!
!===================================================================================================
!
!===================================================================================================
!
    type CalcG_Theta
! ----- name of theta field
        character(len=24)       :: theta_field = ' '
!------ if a input field for theta factors
        aster_logical           :: theta_factors_in = ASTER_FALSE
! ----- name of factors necessary to create theta field in te
        character(len=24)       :: theta_factors = ' '
! ----- name of the matrix from A*G(s)=g(theta)
        character(len=24)       :: matrix = ' '
! ----- number of theta field
        integer(kind=8)                 :: nb_theta_field = 0
! ----- name of crack
        character(len=8)        :: crack = ' '
! ----- type of crack
        character(len=24)       :: crack_type = ' '
! ----- initial configuration of the crack
        character(len=8)        :: config_init = ' '
! ----- name of the mesh (support of crack)
        character(len=8)        :: mesh = ' '
!------ mesh linear or quadratic
        aster_logical           :: milieu = ASTER_FALSE
! ----- name of the nodes of the mesh (support of crack)
        character(len=24)       :: nomNoeud = ' '
! ----- number of nodes in the crack
        integer(kind=8)                 :: nb_fondNoeud = 0
! ----- radius
        real(kind=8)            :: r_inf = 0.d0, r_sup = 0.d0
! ----- radius function of curvilinear abcissa
        character(len=8)        :: r_inf_fo = ' ', r_sup_fo = ' '
! ----- How is theta torus defined ?
        character(len=24)       :: radius_type = ' '
! ----- lenght of crack
        real(kind=8)            :: lonfis = 0.d0
! ----- number of layer
        integer(kind=8)                 :: nb_couche_inf = 0, nb_couche_sup = 0
! ----- type of discretization
        character(len=8)        :: discretization = ' '
! ----- nubmer of nodes if nb_pts_fond is defined (for linear discretization)
        integer(kind=8)                 :: nb_point_fond = 0
! ----- name of fields for nb_point_fond
        character(len=19)       :: npf_baseloc = '&&CALC_G.BASLOC'
        character(len=24)       :: npf_basefond = '&&CALC_G.BASFON'
        character(len=24)       :: npf_fondNoeud = '&&CALC_G.NOEUF'
        character(len=24)       :: npf_absfon = '&&CALC_G.ABSFON'
        character(len=24)       :: npf_abscurv = '&&CALC_G.ABSCUR'
! ----- number of nodes (for linear discretization)
        integer(kind=8)                 :: nnof = 0
! ----- nubmer of nodes (for legendre discretization)
        integer(kind=8)                 :: degree = 0
!-------the crack is symetric ?
        character(len=8)        :: symech = 'NON'
! ----- the crack is closed ?
        aster_logical           :: l_closed = ASTER_FALSE
! ----- name of the curvature
        character(len=24)       :: curvature = ' '
! ----- name of the curvilinear abscissa of crack nodes
        character(len=24)       :: absfond = ' '
! ----- id of crack nodes
        character(len=24)       :: fondNoeudNume = ' '
! ----- name of coordinates of nodes in the crack
        character(len=24)       :: fondNoeudCoor = '&&CALC_G.COORN'
! ----- circular or elliptical crack
        character(len=24)       :: form_fiss = ''
! ----- parameters of curve crack
        real(kind=8)            :: rayon = 0.
        real(kind=8)            :: demi_grand_axe = 0.
        real(kind=8)            :: demi_petit_axe = 0.
! ----- member function
    contains
        procedure, pass    :: initialize => initialize_theta
        procedure, pass    :: print => print_theta
        procedure, pass, private :: create_npf
        procedure, pass    :: compute_curvature
        procedure, pass    :: getCoorNodes
        procedure, pass    :: getAbscurv
        procedure, pass    :: getAbsfon
        procedure, pass    :: getAbsfonName
        procedure, pass    :: getBaseLoc
        procedure, pass    :: getFondTailleR
        procedure, pass    :: getFondNoeud
        procedure, pass    :: getFondNoeudCoor
        procedure, pass    :: getFondNoeudNume
    end type CalcG_Theta
!
!=================================================================================================
!
!=================================================================================================
!
    type CalcG_Table
! ----- name of table G (out)
        character(len=24)  :: table_g = ' '
! ----- list of parameters
        character(len=24)  :: list_name_para(NB_MAX_PARA) = ' '
        character(len=8)   :: list_type_para(NB_MAX_PARA) = ' '
        integer(kind=8)            :: nb_para = 0
! ----- grandeur
        real(kind=8), pointer :: v_G(:) => null()
        real(kind=8), pointer :: v_K1(:) => null()
        real(kind=8), pointer :: v_K2(:) => null()
        real(kind=8), pointer :: v_K3(:) => null()
        real(kind=8), pointer :: v_G_IRWIN(:) => null()
        real(kind=8), pointer :: v_G_EPSI(:) => null()
        real(kind=8), pointer :: v_KJ(:) => null()
        real(kind=8), pointer :: v_KJ_EPSI(:) => null()
        real(kind=8), pointer :: v_TEMP(:) => null()
        character(len=8), pointer :: v_COMPOR(:) => null()
! ----- nb point
        integer(kind=8) :: nb_point = 1
!
! ----- member function
    contains
        procedure, pass    :: initialize => initialize_table
        procedure, pass    :: addValues
        procedure, pass    :: addPara
        procedure, pass    :: save => save_table
        procedure, pass    :: clean => clean_table
!
    end type CalcG_Table
!
!===================================================================================================
!
!===================================================================================================
!
    type CalcG_Stat
! ----- level information
        integer(kind=8)            :: level_info = 1
        real(kind=8)       :: time
! ----- cgField
        real(kind=8)       :: init_cgField = 0.d0
        real(kind=8)       :: clean_cgField = 0.d0
! ----- cgTheta
        real(kind=8)       :: init_cgTheta = 0.d0
        real(kind=8)       :: npf_cgTheta = 0.d0
        integer(kind=8)            :: nb_npf_cgTheta = 0
! ----- cgStudy
        real(kind=8)       :: init_cgStudy = 0.d0
        integer(kind=8)            :: nb_init_cgStudy = 0
! ----- cgTable
        real(kind=8)       :: init_cgTable = 0.d0
        real(kind=8)       :: clean_cgTable = 0.d0
        real(kind=8)       :: save_cgTable = 0.d0
        integer(kind=8)            :: nb_save_cgTable = 0
! ----- cgComputeGtheta
        real(kind=8)       :: cgCmpGtheta = 0.d0
        integer(kind=8)            :: nb_cgCmpGtheta = 0
        real(kind=8)       :: cgCmpGtheta_te = 0.d0
        integer(kind=8)            :: nb_cgCmpGtheta_te = 0
        real(kind=8)       :: cgCmpGtheta_sys = 0.d0
        integer(kind=8)            :: nb_cgCmpGtheta_sys = 0
        real(kind=8)       :: cgCmpGtheta_disc = 0.d0
        integer(kind=8)            :: nb_cgCmpGtheta_disc = 0
        real(kind=8)       :: cgCmpGtheta_mes = 0.d0
        integer(kind=8)            :: nb_cgCmpGtheta_mes = 0
! ----- routines
        real(kind=8)       :: cgVerif = 0.d0
        real(kind=8)       :: cgThetaFact = 0.d0
        real(kind=8)       :: cgCmpMat = 0.d0
        real(kind=8)       :: cgExpTabl = 0.d0

! ----- member function
    contains
        procedure, pass    :: initialize => initialize_stat
        procedure, pass    :: print => print_stat
        procedure, pass    :: finish => finish_stat
    end type CalcG_Stat
!
contains
!
!---------------------------------------------------------------------------------------------------
! -- member functions
!---------------------------------------------------------------------------------------------------
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_stat(this)
!
        implicit none
!
        class(CalcG_Stat), intent(inout)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Stat type
!   In this     : calG Stat
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ifm
!
! --- Level of information
!
        call infniv(ifm, this%level_info)
        call cpu_time(this%time)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_field(this, cgStat)
!
        implicit none
!
        class(CalcG_Field), intent(inout)  :: this
        type(CalcG_Stat), intent(inout)   :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Field type
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ier, ifm, nbropt
        character(len=3) :: repk
        character(len=16) :: k16bid
        character(len=19) :: lisopt
        character(len=8)  :: modele, mater
        character(len=8)  :: templist(NB_MAX_OPT)
        integer(kind=8), pointer :: v_nume(:) => null()
        real(kind=8) :: start, finish
!
        call cpu_time(start)
        templist = ' '
!
        call jemarq()
!
! --- Concept de sortie (table container)
!
        call getres(this%table_out, k16bid, k16bid)
!
! --- Get name and type of result (in)
!
        call getvid(' ', 'RESULTAT', scal=this%result_in, nbret=ier)
        ASSERT(ier == 1)
        call gettco(this%result_in, this%result_in_type, ASTER_TRUE)
        call dismoi('MODELE', this%result_in, 'RESULTAT', repk=modele)
        call dismoi('DIM_GEOM', modele, 'MODELE', repi=this%ndim)
!
        call dismoi('CHAM_MATER', this%result_in, 'RESULTAT', repk=mater)
        call dismoi('EXI_VARC', mater, 'CHAM_MATER', repk=repk)
        if (repk == "OUI") then
            this%l_temp = ASTER_TRUE
        end if
!
! --- Get name of option
!
        call getvtx(' ', 'OPTION', nbret=ier)
        if (ier == 1) then
            this%nb_option = ier
            call getvtx(' ', 'OPTION', scal=this%list_option(1))
        else
            this%nb_option = -ier
            ASSERT(this%nb_option <= NB_MAX_OPT)
            call getvtx(' ', 'OPTION', nbval=this%nb_option, vect=this%list_option)
!
! --- Remove option G* if option KJ* exists
!
            if (any('G' == this%list_option) .and. any('KJ' == this%list_option)) then
                templist = ' '
                this%nb_option = this%nb_option-1
                templist = pack(this%list_option, this%list_option .ne. 'G')
                this%list_option = templist
            end if
            if (any('G_EPSI' == this%list_option) .and. any('KJ_EPSI' == this%list_option)) then
                templist = ' '
                this%nb_option = this%nb_option-1
                templist = pack(this%list_option, this%list_option .ne. 'G_EPSI')
                this%list_option = templist
            end if
        end if
!
! --- Level of information
!
        call infniv(ifm, this%level_info)
!
! --- List of nume
!
        this%list_nume_name = '&&OP0027.VECTORDR'
        call cgcrio(this%result_in, this%list_nume_name, this%nb_nume)
        ASSERT(this%nb_nume > 0)
        call jeveuo(this%list_nume_name, 'L', vi=v_nume)
        AS_ALLOCATE(vi=this%list_nume, size=this%nb_nume)
        this%list_nume(1:this%nb_nume) = v_nume(1:this%nb_nume)
!
! --- Read <CARTE> COMPORTEMENT
!
        call cgReadCompor(this%result_in, this%compor, this%list_nume(1), this%l_incr)
!
        if (this%level_info > 1) then
            call comp_info(modele, this%compor)
        end if
!
! --- if ELAS_INCR
!
        if (this%l_incr) then
            lisopt = ' '
            nbropt = 0
            call ccbcop(this%result_in, this%result_out, this%list_nume_name, &
                        this%nb_nume, lisopt, nbropt)
        end if
!
        call jedema()
!
        call cpu_time(finish)
        cgStat%init_cgField = cgStat%init_cgField+finish-start
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine print_field(this)
!
        implicit none
!
        class(CalcG_Field), intent(in)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Field type
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i
!
        print *, "----------------------------------------------------------------------"
        print *, "Informations about CalcG_Field"
        print *, "Level of informations: ", this%level_info
        print *, "Dimension of the problem: ", this%ndim
        print *, "Result: ", this%result_in, " of type ", this%result_in_type
        print *, "Output: ", this%table_out
        print *, "Number of option to compute: ", this%nb_option
        do i = 1, this%nb_option
            print *, "** option ", i, ": ", this%list_option(i)
        end do
        print *, "Number of step/mode to compute: ", this%nb_nume
        print *, "----------------------------------------------------------------------"
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine print_stat(this)
!
        implicit none
!
        class(CalcG_Stat), intent(in)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Stat type
!   In this     : calG stat
! --------------------------------------------------------------------------------------------------
!
        print *, "---------------------------------------------------------------------------------"
        print *, "Informations about CalcG_Stat"
        print *, "Level of informations: ", this%level_info
        print *, "Total time of OP: ", this%time, "s"
        print *, "*CalcG_Field: "
        print *, "**initialize(): ", this%init_cgField, "s with ", 1, "calls"
        print *, "**clean():      ", this%clean_cgField, "s with ", 1, "calls"
        print *, "*CalcG_Theta: "
        print *, "**initialize(): ", this%init_cgTheta, "s with ", 1, "calls"
        print *, "**create_npf(): ", this%npf_cgTheta, "s with ", this%nb_npf_cgTheta, "calls"
        print *, "*CalcG_Study: "
        print *, "**initialize(): ", this%init_cgStudy, "s with ", this%nb_init_cgStudy, "calls"
        print *, "*CalcG_Table: "
        print *, "**initialize(): ", this%init_cgTable, "s with ", 1, "calls"
        print *, "**save():       ", this%save_cgTable, "s with ", this%nb_save_cgTable, "calls"
        print *, "**clean():      ", this%clean_cgTable, "s with ", 1, "calls"
        print *, "*Others: "
        print *, "**cgVerification():       ", this%cgVerif, "s with ", 1, "calls"
        print *, "**cgComputeThetaFactor(): ", this%cgThetaFact, "s with ", 1, "calls"
        print *, "**cgComputeMatrix():      ", this%cgCmpMat, "s with ", 1, "calls"
        print *, "**cgExportTableG():       ", this%cgExpTabl, "s with ", 1, "calls"
        print *, "**cgComputeGtheta():      ", this%cgCmpGtheta, "s with ", &
            this%nb_cgCmpGtheta, "calls"
        print *, "***calcul():       ", this%cgCmpGtheta_te, "s with ", &
            this%nb_cgCmpGtheta_te, "calls"
        print *, "***gsyste():       ", this%cgCmpGtheta_sys, "s with ", &
            this%nb_cgCmpGtheta_sys, "calls"
        print *, "***cgDiscrField(): ", this%cgCmpGtheta_disc, "s with ", &
            this%nb_cgCmpGtheta_disc, "calls"
        print *, "***mesomm():       ", this%cgCmpGtheta_mes, "s with ", &
            this%nb_cgCmpGtheta_mes, "calls"
        print *, "---------------------------------------------------------------------------------"
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine finish_stat(this)
!
        implicit none
!
        class(CalcG_Stat), intent(inout)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Stat type
!   In this     : calG Stat
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: finish
!
        call cpu_time(finish)
        this%time = finish-this%time
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function isModeMeca(this) result(lmode)
!
        implicit none
!
        class(CalcG_Field), intent(in)  :: this
        aster_logical                   :: lmode
!
! --------------------------------------------------------------------------------------------------
!
!   The result is a MODE_MECA ?
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        if (this%result_in_type == "MODE_MECA") then
            lmode = ASTER_TRUE
        else
            lmode = ASTER_FALSE
        end if
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function isDynaTrans(this) result(lmode)
!
        implicit none
!
        class(CalcG_Field), intent(in)  :: this
        aster_logical                   :: lmode
!
! --------------------------------------------------------------------------------------------------
!
!   The result is a DYNA_TRANS ?
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        if (this%result_in_type == "DYNA_TRANS") then
            lmode = ASTER_TRUE
        else
            lmode = ASTER_FALSE
        end if
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine clean_field(this, cgStat)
!
        implicit none
!
        class(CalcG_Field), intent(inout)  :: this
        type(CalcG_stat), intent(inout) :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!   Clean objects
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret
        integer(kind=8), pointer :: ordr(:) => null()
        real(kind=8) :: start, finish
!
        call cpu_time(start)
!
        call jemarq()
!
        if (this%l_incr) then
!
            call jeexin(this%result_out//'           .ORDR', iret)
            if (iret .ne. 0) then
                call jeveuo(this%result_out//'           .ORDR', 'L', vi=ordr)
                call rsrusd(this%result_out, ordr(1))
                call detrsd('RESULTAT', this%result_out)
            end if
!
            call jedetr(this%list_nume_name)
            call jedetr(this%result_out)
            call rsmena(this%result_in)
        end if
!
        AS_DEALLOCATE(vi=this%list_nume)
!
        call jedema()
!
        call cpu_time(finish)
        cgStat%clean_cgField = cgStat%clean_cgField+finish-start
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine clean_table(this, cgStat)
!
        implicit none
!
        class(CalcG_Table), intent(inout)  :: this
        type(CalcG_stat), intent(inout) :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!   Clean objects
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: start, finish
!
        call cpu_time(start)
!
        call jemarq()
!
        AS_DEALLOCATE(vr=this%v_G)
        AS_DEALLOCATE(vr=this%v_K1)
        AS_DEALLOCATE(vr=this%v_K2)
        AS_DEALLOCATE(vr=this%v_K3)
        AS_DEALLOCATE(vr=this%v_G_IRWIN)
        AS_DEALLOCATE(vr=this%v_G_EPSI)
        AS_DEALLOCATE(vr=this%v_KJ)
        AS_DEALLOCATE(vr=this%v_KJ_EPSI)
        AS_DEALLOCATE(vk8=this%v_COMPOR)
        AS_DEALLOCATE(vr=this%v_TEMP)
!
        call jedema()
!
        call cpu_time(finish)
        cgStat%clean_cgTable = cgStat%clean_cgTable+finish-start
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_study(this, result_in, nume_index, cgStat)
!
        implicit none
!
        class(CalcG_Study), intent(inout)  :: this
        character(len=8), intent(in)       :: result_in
        integer(kind=8), intent(in)                :: nume_index
        type(CalcG_Stat), intent(inout)    :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Study type
!   In this     : study type
!   In nume_index : index nume
! --------------------------------------------------------------------------------------------------
!
        real(kind=8) :: start, finish
        character(len=3) :: exiInco, exiAxis
!
        call cpu_time(start)
!
        this%nume_ordre = nume_index
        this%loading = "&&STUDY.CHARGES"
        call medomg(result_in, this%nume_ordre, this%model, this%material, &
                    this%mateco, this%loading)
        call dismoi('CARA_ELEM', result_in, 'RESULTAT', repk=this%carael)
        call dismoi('NOM_MAILLA', this%model, 'MODELE', repk=this%mesh)
!       Cas axisymétrique
        call dismoi('EXI_AXIS', this%model, 'MODELE', repk=exiAxis)
        if (exiAxis(1:3) .eq. 'OUI') then
            this%l_axis = ASTER_TRUE
        else
            this%l_axis = ASTER_FALSE
        end if
!
        call dismoi('EXI_INCO', this%model, 'MODELE', repk=exiInco)
        if (exiInco(1:3) .eq. 'OUI') then
            this%l_exi_inco = ASTER_TRUE
        else
            this%l_exi_inco = ASTER_FALSE
        end if
!
        call cpu_time(finish)
        cgStat%init_cgStudy = cgStat%init_cgStudy+finish-start
        cgStat%nb_init_cgStudy = cgStat%nb_init_cgStudy+1
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine print_study(this)
!
        implicit none
!
        class(CalcG_Study), intent(in)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   print informations of a CalcG_Study type
!   In this     : study type
! --------------------------------------------------------------------------------------------------
!
        print *, "----------------------------------------------------------------------"
        print *, "Informations about CalcG_Study"
        print *, "Option: ", this%option
        print *, "Model: ", this%model
        print *, "Mesh: ", this%mesh
        print *, "Material: ", this%material
        print *, "Coded material: ", this%mateco
        print *, "Loading: ", this%loading
        print *, "Nume index: ", this%nume_ordre
        print *, "----------------------------------------------------------------------"
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine setOption(this, option, isModeMeca)
!
        implicit none
!
        class(CalcG_Study), intent(inout)  :: this
        character(len=8), intent(in)       :: option
        aster_logical, intent(in)          :: isModeMeca
!
! --------------------------------------------------------------------------------------------------
!   print informations of a CalcG_Study type
!   In this     : study type
!   In option   : name of option
! --------------------------------------------------------------------------------------------------
!
        this%option = option
!
        if (isModeMeca) then
            if (this%option .eq. 'K') then
                this%l_modal = ASTER_TRUE
            else
                call utmess('F', 'RUPTURE0_27')
            end if
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getField(this, result_in)
!
        implicit none
!
        class(CalcG_Study), intent(inout)  :: this
        character(len=8), intent(in)       :: result_in
!
! --------------------------------------------------------------------------------------------------
!   print informations of a CalcG_Study type
!   In this     : study type
!   In result_in   : name of result field
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: iret
!
        call rsexch('F', result_in, 'DEPL', this%nume_ordre, this%depl, iret)
        call rsexch(' ', result_in, 'VITE', this%nume_ordre, this%vitesse, iret)
        if (iret .ne. 0) then
            this%vitesse = ' '
            this%acce = ' '
        else
            call rsexch(' ', result_in, 'ACCE', this%nume_ordre, this%acce, iret)
        end if
!
    end subroutine
!
    !
!===================================================================================================
!
!===================================================================================================
!
    subroutine getParameter(this, result_in)
!
        implicit none
!
        class(CalcG_Study), intent(inout)  :: this
        character(len=8), intent(in)       :: result_in
!
! --------------------------------------------------------------------------------------------------
!   print informations of a CalcG_Study type
!   In this     : study type
!   In result_in   : name of result field
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ipuls, jinst
        character(len=8)  :: k8bid

!
        if (this%l_modal) then
            call rsadpa(result_in, 'L', 1, 'OMEGA2', this%nume_ordre, &
                        0, sjv=ipuls, styp=k8bid)
            this%pulse = sqrt(zr(ipuls))
            this%time = 0.d0
        else
            call rsadpa(result_in, 'L', 1, 'INST', this%nume_ordre, &
                        0, sjv=jinst, styp=k8bid)
            this%pulse = 0.d0
            this%time = zr(jinst)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_theta(this, cgStat)
!
        implicit none
!
        class(CalcG_Theta), intent(inout)  :: this
        type(CalcG_Stat), intent(inout)   :: cgStat
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Theta type
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ier, i, ndim, jma, nume
        character(len=8) :: typfon, nompar(1), typma
        real(kind=8) :: maxtai, mintai, valpar(1), valres_i, valres_s
        aster_logical :: l_disc
        real(kind=8), pointer :: fondTailleR(:) => null()
        real(kind=8), pointer :: absfon(:) => null()
        integer(kind=8), pointer :: typmail(:) => null()
        real(kind=8) :: start, finish
        character(len=8) :: thetafactorsin
!
        call cpu_time(start)
!
        call jemarq()

! --- get automatic name
        call gcncon("_", this%theta_field)
        call gcncon("_", this%matrix)
        call gcncon("_", this%fondNoeudNume)
!
! --- get informations about the crack
!
        call getvtx('THETA', 'FISSURE', iocc=1, scal=this%crack, nbret=ier)
        ASSERT(ier == 1)
        call gettco(this%crack, this%crack_type, ASTER_TRUE)
!
        call dismoi('NOM_MAILLA', this%crack, 'FOND_FISS', arret='F', repk=this%mesh)
        call dismoi('TYPE_FOND', this%crack, 'FOND_FISS', arret='F', repk=typfon)
        call dismoi('CONFIG_INIT', this%crack, 'FOND_FISS', repk=this%config_init)
        call dismoi('SYME', this%crack, 'FOND_FISS', repk=this%symech)
        call dismoi('DIM_GEOM', this%mesh, 'MAILLAGE', repi=ndim)
!
        ASSERT(.not. isParallelMesh(this%mesh))
!
        this%nomNoeud = this%mesh//'.NOMNOE'
        this%absfond = this%crack//'.ABSFON'
!
!       Maillage linéaire ou quadratique
        call jeveuo(this%crack//'.LEVRESUP.MAIL', 'L', jma)
        call jeveuo(this%mesh//'.TYPMAIL', 'L', vi=typmail)
!
        nume = char8_to_int(zk8(jma))
!
        call jenuno(jexnum('&CATA.TM.NOMTM', typmail(nume)), typma)
!
        if (.not. ismali(typma)) then
            this%milieu = ASTER_TRUE
        end if
!
! --- the crack is closed ?
        if (typfon .eq. 'FERME') then
            this%l_closed = ASTER_TRUE
        else
            this%l_closed = ASTER_FALSE
        end if
! --- number of nodes in the crack
        call jelira(this%crack//'.FOND.NOEU', 'LONMAX', this%nb_fondNoeud)
!
! --- get informations about theta discretization
!
        call getvtx('THETA', 'DISCRETISATION', iocc=1, scal=this%discretization, nbret=ier)
        l_disc = (this%discretization == "LINEAIRE") .or. (this%discretization == "LEGENDRE")
        ASSERT(l_disc)
!
        if (this%discretization == "LINEAIRE") then
            call getvis('THETA', 'NB_POINT_FOND', iocc=1, scal=this%nb_point_fond, nbret=ier)
            ASSERT(this%nb_point_fond >= 0)
            ASSERT(this%degree == 0)
            if (this%nb_point_fond .ne. 0) then
                this%milieu = ASTER_FALSE
                this%nnof = this%nb_point_fond
            else
                this%nnof = this%nb_fondNoeud
            end if
        end if
!
        if (this%discretization == "LEGENDRE") then
            call getvis('THETA', 'DEGRE', iocc=1, scal=this%degree, nbret=ier)
            this%nnof = this%nb_fondNoeud
            ASSERT(this%nb_point_fond == 0)
            ASSERT(this%degree >= 0)
            if (this%l_closed) then
                call utmess('F', 'RUPTURE0_90')
            end if
        end if
!
! --- Find nb_theta_field
        if (ndim == 2) then
            this%nb_theta_field = 1
        else
            if (this%discretization .eq. 'LINEAIRE') then
                this%nb_theta_field = this%nnof
            elseif (this%discretization .eq. 'LEGENDRE') then
                this%nb_theta_field = this%degree+1
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
!
! ---- Create pointers for NB_POINT_FOND
        call this%create_npf(cgStat)
!
! ---- Theta torus definition
!
        call getvis('THETA', 'NB_COUCHE_INF', iocc=1, scal=this%nb_couche_inf, nbret=ier)
        call getvis('THETA', 'NB_COUCHE_SUP', iocc=1, scal=this%nb_couche_sup, nbret=ier)
        if (ier .ne. 0) then
            this%radius_type = 'NB_COUCHE'
            !       Blindage pour interdire NB_COUCHE avec HEXA27 et PENTA18
            !       temporaire à supprimer quand issue 33988 sera corrigée
            call jeveuo(this%crack//'.LEVRESUP.MAIL', 'L', jma)
            call jeveuo(this%mesh//'.TYPMAIL', 'L', vi=typmail)
            nume = char8_to_int(zk8(jma))
            call jenuno(jexnum('&CATA.TM.NOMTM', typmail(nume)), typma)
!
            if ((typma .eq. 'QUAD9') .or. (typma .eq. 'TRIA7')) then
                call utmess('F', 'RUPTURE1_89')
            end if
!
        end if
!
        if (ier == 1 .and. ((this%nb_couche_inf < 0) .or. &
                            (this%nb_couche_inf >= this%nb_couche_sup))) then
            call utmess('F', 'RUPTURE3_4', ni=2, vali=[this%nb_couche_inf, this%nb_couche_sup])
        end if

        call this%getAbsfon(absfon)
        this%lonfis = absfon(this%nnof)

        ier = 0
        call getvtx('THETA', 'CHAM_THETA', iocc=1, scal=this%theta_factors, nbret=ier)
        if (ier > 0) then
            call exisd('CHAM_NO', this%theta_factors, ier)
        end if

        if (ier == 0) then
! --- get automatic name
            call gcncon("_", this%theta_factors)
!
! --- Get RINF and RSUP from command file or from SD FOND_FISSURE

            call getvr8('THETA', 'R_INF', iocc=1, scal=this%r_inf, nbret=ier)
            call getvr8('THETA', 'R_SUP', iocc=1, scal=this%r_sup, nbret=ier)
            if (ier .ne. 0) then
                this%radius_type = 'R'
            end if
!
            call getvid('THETA', 'R_INF_FO', iocc=1, scal=this%r_inf_fo, nbret=ier)
            call getvid('THETA', 'R_SUP_FO', iocc=1, scal=this%r_sup_fo, nbret=ier)
            if (ier .ne. 0) then
                this%radius_type = 'R_FO'
            end if
!
!       Verifications
            if (this%radius_type == 'R' .and. &
                ((this%r_inf < 0.d0) .or. (this%r_inf >= this%r_sup))) then
                call utmess('F', 'RUPTURE3_3', nr=2, valr=[this%r_inf, this%r_sup])
            end if
!
            if (this%radius_type == 'R_FO') then
                do i = 1, this%nnof
                    nompar(1) = 'ABSC'
                    valpar(1) = absfon(i)
                    call fointe('FM', this%r_inf_fo, 1, nompar, valpar, valres_i, ier)
                    call fointe('FM', this%r_sup_fo, 1, nompar, valpar, valres_s, ier)
                    if (valres_s .le. valres_i) then
                        call utmess('F', 'RUPTURE1_6')
                    end if
                end do
            end if
!
!       if no radius is defined
            if (this%radius_type .ne. 'R' .and. this%radius_type .ne. 'R_FO' &
                .and. this%radius_type .ne. 'NB_COUCHE') then
                this%radius_type = 'R'
                if (this%config_init == 'DECOLLEE') then
!               A vérifier si toujours impossible de calculer r dans ce cas
!               (fond_taille_r disponible dans la sd_fond_fissure même si
!               DECOLLEE?)
                    call utmess('F', 'RUPTURE1_7')
                end if
                call this%getFondTailleR(fondTailleR)
                maxtai = fondTailleR(1)
                mintai = fondTailleR(1)
                do i = 1, this%nb_fondNoeud
                    maxtai = max(maxtai, fondTailleR(i))
                    mintai = min(mintai, fondTailleR(i))
                end do
                this%r_inf = 2*maxtai
                this%r_sup = 4*maxtai
                call utmess('I', 'RUPTURE1_5', nr=2, valr=[this%r_inf, this%r_sup])
                if (maxtai .gt. 2*mintai) then
                    call utmess('A', 'RUPTURE1_16', nr=2, valr=[mintai, maxtai])
                end if
            end if
        else
            this%theta_factors_in = ASTER_TRUE
!           pour COPIER ET STOCKER LE INPUT CHAMP_NO THETA_FACTORS
            thetafactorsin = this%theta_factors(1:8)
            this%theta_factors = this%theta_factors(1:8)//'_CHAM_THETA_FACT'
            call copisd('CHAMP_GD', 'G', thetafactorsin, this%theta_factors)
        end if
!
! --- Obtain crack parameters to take into account the second order due to crack curvature
        call getvtx('', 'FORM_FISS', iocc=1, scal=this%form_fiss, nbret=ier)
        if (this%form_fiss .eq. 'CERCLE') then
            call getvr8('', 'RAYON', iocc=1, scal=this%rayon, nbret=ier)
        else
            call getvr8('', 'DEMI_GRAND_AXE', iocc=1, scal=this%demi_grand_axe, nbret=ier)
            call getvr8('', 'DEMI_PETIT_AXE', iocc=1, scal=this%demi_petit_axe, nbret=ier)
        end if

        call jedema()
!
        call cpu_time(finish)
        cgStat%init_cgTheta = cgStat%init_cgTheta+finish-start
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine compute_curvature(this, model)
!
        implicit none
!
        class(CalcG_Theta), intent(inout)  :: this
        character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
!   Compute the curvature in 3D
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        character(len=24) :: baseloc
!
        baseloc = this%crack//'.BASLOC'
        this%curvature = '&&cgtheta.COURB'
        call xcourb(baseloc, this%mesh, model, this%curvature)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getCoorNodes(this, v_coor)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        real(kind=8), pointer :: v_coor(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on coordinates of nodes
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        call jeveuo(this%mesh//'.COORDO    .VALE', 'L', vr=v_coor)
        call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getFondNoeudCoor(this, v_coor)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        real(kind=8), pointer :: v_coor(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on coordinates of nodes
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        call jeveuo(this%fondNoeudCoor, 'L', vr=v_coor)
        call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getAbscurv(this, v_abs)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        real(kind=8), pointer :: v_abs(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on abscisse curviligne
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        if (this%discretization == "LINEAIRE" .and. this%nb_point_fond > 0) then
            call jeveuo(this%npf_abscurv, 'L', vr=v_abs)
        else
            call jeveuo(this%crack//'.ABSCUR', 'L', vr=v_abs)
        end if
        call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getAbsfon(this, v_absfon)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        real(kind=8), pointer :: v_absfon(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on abscisse curviligne
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        if (this%discretization == "LINEAIRE" .and. this%nb_point_fond > 0) then
            call jeveuo(this%npf_absfon, 'L', vr=v_absfon)
        else
            call jeveuo(this%crack//'.ABSFON', 'L', vr=v_absfon)
        end if
        call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getAbsfonName(this, absfon)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        character(len=24), intent(out) :: absfon
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on abscisse curviligne
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        if (this%discretization == "LINEAIRE" .and. this%nb_point_fond > 0) then
            absfon = this%npf_absfon
        else
            absfon = this%crack//'.ABSFON'
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getFondTailleR(this, v_taille)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        real(kind=8), pointer :: v_taille(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        call jeveuo(this%crack//'.FOND.TAILLE_R', 'L', vr=v_taille)
        call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine getBaseLoc(this, v_base)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        real(kind=8), pointer :: v_base(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on baseloc
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        if (this%discretization == "LINEAIRE" .and. this%nb_point_fond > 0) then
            call jeveuo(this%npf_baseloc//'.VALE', 'L', vr=v_base)
        else
            call jeveuo(this%crack//'.BASLOC    .VALE', 'L', vr=v_base)
        end if
        call jedema()
!
    end subroutine
!===================================================================================================
!
    subroutine getFondNoeud(this, v_fondnoeu)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        character(len=8), pointer :: v_fondnoeu(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on baseloc
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        if (this%discretization == "LINEAIRE" .and. this%nb_point_fond > 0) then
            call jeveuo(this%npf_fondNoeud, 'L', vk8=v_fondnoeu)
        else
            call jeveuo(this%crack//'.FOND.NOEU', 'L', vk8=v_fondnoeu)
        end if
        call jedema()
!
    end subroutine
!===================================================================================================
!
    subroutine getFondNoeudNume(this, v_fondnoeuNume)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        integer(kind=8), pointer :: v_fondnoeuNume(:)
!
! --------------------------------------------------------------------------------------------------
!
!   Get pointer on baseloc
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        call jemarq()
        call jeveuo(this%fondNoeudNume, 'L', vi=v_fondnoeuNume)
        call jedema()
!
    end subroutine
!===================================================================================================
!
!===================================================================================================
!
    subroutine print_theta(this)
!
        implicit none
!
        class(CalcG_Theta), intent(in)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   print informations of a CalcG_Theta type
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        print *, "----------------------------------------------------------------------"
        print *, "Informations about CalcG_Theta"
        print *, "Field theta: ", this%theta_field
        print *, "theta_factors: ", this%theta_factors
        print *, "Crack: ", this%crack, " of type ", this%crack_type
        print *, "Number of nodes in the crack: ", this%nb_fondNoeud
        print *, "Mesh support: ", this%mesh
        print *, "linear or quadratic", this%milieu
        print *, "Initial configuration: ", this%config_init
        print *, "the crack is symetric: ", this%symech
        print *, "The crack is closed ?: ", this%l_closed
!        print*, "Nombre de champs THETA: ", this%nb_theta_field
        print *, "Discretization : ", this%discretization, " with number/degree ", &
            this%nnof, this%degree
        print *, "Radius:"
        print *, "*** Inferior: ", this%r_inf
        print *, "*** Superior: ", this%r_sup
        print *, "Number of cell layers:"
        print *, "*** Inferior: ", this%nb_couche_inf
        print *, "*** Superior: ", this%nb_couche_sup
        print *, "----------------------------------------------------------------------"
!
    end subroutine
!
!=======================================================================================
!
!=======================================================================================
!
    subroutine addPara(this, name, type)
!
        implicit none
!
        class(CalcG_Table), intent(inout)  :: this
        character(len=*), intent(in) :: name, type
!
! ---------------------------------------------------------------------------------------
!
!   add a parameter to the table
!   In this     : calcG Table
! ---------------------------------------------------------------------------------------
        this%nb_para = this%nb_para+1
        ASSERT(this%nb_para .le. NB_MAX_PARA)
!
        this%list_name_para(this%nb_para) = name
        this%list_type_para(this%nb_para) = type

    end subroutine
!
!=======================================================================================
!
!=======================================================================================
!
    subroutine initialize_table(this, cgField, cgTheta, cgStat)
!
        implicit none
!
        class(CalcG_Table), intent(inout)  :: this
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
        type(CalcG_Stat), intent(inout)   :: cgStat
!
! --------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Table type
!   In this     : calcG Table
! --------------------------------------------------------------------------------------
        integer(kind=8) :: iopt, nbValues
        character(len=8) :: option
        integer(kind=8), pointer :: fondNoeudNume(:) => null()
        real(kind=8) :: start, finish
!
        call cpu_time(start)
        call jemarq()
!
! --- Table pour les valeurs (table)
!
        call gcncon("_", this%table_g)
        call tbcrsd(this%table_g, 'G')
!
        this%nb_point = cgTheta%nnof
        this%nb_para = 0
!
! --- INST or FREQ
        if (cgField%isModeMeca()) then
            call this%addPara('NUME_MODE', 'I')
        else
            call this%addPara('NUME_ORDRE', 'I')
            call this%addPara('INST', 'R')
        end if
! --- Node name
        call this%addPara('NOEUD', 'K8')
        call this%addPara('NUM_PT', 'I')
! --- Coordinates of nodes
        call this%addPara('COOR_X', 'R')
        call this%addPara('COOR_Y', 'R')
        if (cgField%ndim .eq. 3) then
            call this%addPara('COOR_Z', 'R')
            call this%addPara('ABSC_CURV', 'R')
            call this%addPara('ABSC_CURV_NORM', 'R')
        end if
! --- Tempature
        if (cgField%l_temp) then
!            call this%addPara('TEMP', 'R')
            AS_ALLOCATE(vr=this%v_TEMP, size=this%nb_point)
        end if
! --- Behavior
        call this%addPara('COMPORTEMENT', 'K8')
        call cgTheta%getFondNoeudNume(fondNoeudNume)
        AS_ALLOCATE(vk8=this%v_COMPOR, size=this%nb_point)
        call cgComporNodes(cgField%result_in, cgField%list_nume(1), this%nb_point, &
                           fondNoeudNume, this%v_COMPOR)
! --- Option
        nbValues = this%nb_point
        do iopt = 1, cgField%nb_option
            option = cgField%list_option(iopt)

            if (option == "G") then
                call this%addPara('G', 'R')
            elseif (option == "K") then
                call this%addPara('K1', 'R')
                call this%addPara('K2', 'R')
                if (cgField%ndim .eq. 3) then
                    call this%addPara('K3', 'R')
                end if
                call this%addPara('G_IRWIN', 'R')
            elseif (option == "KJ") then
                call this%addPara('KJ', 'R')
                call this%addPara('G', 'R')
            elseif (option == "KJ_EPSI") then
                call this%addPara('KJ_EPSI', 'R')
                call this%addPara('G_EPSI', 'R')
            elseif (option == "G_EPSI") then
                call this%addPara('G_EPSI', 'R')
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
!
        AS_ALLOCATE(vr=this%v_G, size=nbValues)
        this%v_G(:) = r8vide()
        AS_ALLOCATE(vr=this%v_K1, size=nbValues)
        this%v_K1(:) = r8vide()
        AS_ALLOCATE(vr=this%v_K2, size=nbValues)
        this%v_K2(:) = r8vide()
        AS_ALLOCATE(vr=this%v_K3, size=nbValues)
        this%v_K3(:) = r8vide()
        AS_ALLOCATE(vr=this%v_G_IRWIN, size=nbValues)
        this%v_G_IRWIN(:) = r8vide()
        AS_ALLOCATE(vr=this%v_G_EPSI, size=nbValues)
        this%v_G_EPSI(:) = r8vide()
        AS_ALLOCATE(vr=this%v_KJ, size=nbValues)
        this%v_KJ(:) = r8vide()
        AS_ALLOCATE(vr=this%v_KJ_EPSI, size=nbValues)
        this%v_KJ_EPSI(:) = r8vide()
!
! --- create table
        call tbajpa(this%table_g, this%nb_para, this%list_name_para, this%list_type_para)
        call jedema()
!
        call cpu_time(finish)
        cgStat%init_cgTable = cgStat%init_cgTable+finish-start
!
    end subroutine
!
!==========================================================================================
!
!==========================================================================================
!
    subroutine addValues(this, cgField, cgStudy, node_id)
!
        implicit none
!
        class(CalcG_Table), intent(inout)  :: this
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_study), intent(in) :: cgStudy
        integer(kind=8), intent(in) :: node_id
!
! ----------------------------------------------------------------------------------------
!
!   add Values in the table
!   In this     : calcG Table
! ----------------------------------------------------------------------------------------
!
        if (cgStudy%option == "G") then
            this%v_G(node_id) = cgStudy%gth(1)
        elseif (cgStudy%option == "K") then
            this%v_K1(node_id) = cgStudy%gth(2)
            this%v_K2(node_id) = cgStudy%gth(3)
            if (cgField%ndim == 3) then
                this%v_K3(node_id) = cgStudy%gth(4)
            end if
            this%v_G_IRWIN(node_id) = cgStudy%gth(5)**2+cgStudy%gth(6)**2+cgStudy%gth(7)**2
        elseif (cgStudy%option == "KJ") then
            if (cgStudy%gth(2) .ge. 0) then
                this%v_KJ(node_id) = sqrt(cgStudy%gth(2))
            else
                this%v_KJ(node_id) = 0.d0
            end if
            this%v_G(node_id) = cgStudy%gth(1)
        elseif (cgStudy%option == "KJ_EPSI") then
            if (cgStudy%gth(2) .ge. 0) then
                this%v_KJ_EPSI(node_id) = sqrt(cgStudy%gth(2))
            else
                this%v_KJ_EPSI(node_id) = 0.d0
            end if
            this%v_G_EPSI(node_id) = cgStudy%gth(1)
        elseif (cgStudy%option == "G_EPSI") then
            this%v_G_EPSI(node_id) = cgStudy%gth(1)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine addValues
!
!
!=======================================================================================
!
!=======================================================================================
!
    subroutine save_table(this, cgField, cgTheta, cgStudy, cgStat)
!
        implicit none
!
        class(CalcG_Table), intent(in)  :: this
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
        type(CalcG_study), intent(in) :: cgStudy
        type(CalcG_stat), intent(inout) :: cgStat
!
! --------------------------------------------------------------------------------------
!
!   save values in the table
!   In this     : calcG Table
! --------------------------------------------------------------------------------------
!
        integer(kind=8)            :: livi(NB_MAX_PARA)
        real(kind=8)       :: livr(NB_MAX_PARA)
        complex(kind=8)    :: livc(NB_MAX_PARA)
        character(len=24)  :: livk(NB_MAX_PARA)
        real(kind=8) :: coor(3), start, finish
        integer(kind=8) :: i_node, iopt
        character(len=8) :: option
        real(kind=8), pointer   :: fondNoeudCoor(:) => null()
        real(kind=8), pointer   :: absfon(:) => null()
        character(len=8), pointer :: fondNoeud(:) => null()
!
        call cpu_time(start)
!
        livi(:) = ismaem()
        livr(:) = r8vide()
        livc(:) = cmplx(r8vide(), r8vide())
        livk(:) = ' '
!
        if (cgField%isModeMeca()) then
            call tbajvi(this%table_g, this%nb_para, 'NUME_MODE', cgStudy%nume_ordre, livi)
        else
            call tbajvi(this%table_g, this%nb_para, 'NUME_ORDRE', cgStudy%nume_ordre, livi)
            call tbajvr(this%table_g, this%nb_para, 'INST', cgStudy%time, livr)
        end if
!
        call cgTheta%getAbsfon(absfon)
        call cgTheta%getFondNoeud(fondNoeud)
        call cgTheta%getFondNoeudCoor(fondNoeudCoor)
!
        do i_node = 1, this%nb_point
            call tbajvk(this%table_g, this%nb_para, 'NOEUD', fondNoeud(i_node), livk)
            call tbajvi(this%table_g, this%nb_para, 'NUM_PT', i_node, livi)
!
            coor = fondNoeudCoor((i_node-1)*3+1:(i_node-1)*3+3)
            call tbajvr(this%table_g, this%nb_para, 'COOR_X', coor(1), livr)
            call tbajvr(this%table_g, this%nb_para, 'COOR_Y', coor(2), livr)
            if (cgField%ndim .eq. 3) then
                call tbajvr(this%table_g, this%nb_para, 'COOR_Z', coor(3), livr)
                call tbajvr(this%table_g, this%nb_para, 'ABSC_CURV', absfon(i_node), livr)
                call tbajvr(this%table_g, this%nb_para, 'ABSC_CURV_NORM', &
                            absfon(i_node)/cgTheta%lonfis, livr)
            end if
!
            if (cgField%l_temp) then
!                call tbajvr(this%table_g, this%nb_para, 'TEMP', this%v_TEMP(i_node), livr)
            end if
            call tbajvk(this%table_g, this%nb_para, 'COMPORTEMENT', this%v_COMPOR(i_node), livk)
!
            do iopt = 1, cgField%nb_option
                option = cgField%list_option(iopt)
                if (option == "G") then
                    call tbajvr(this%table_g, this%nb_para, 'G', this%v_G(i_node), livr)
                elseif (option == "K") then
                    call tbajvr(this%table_g, this%nb_para, 'K1', this%v_K1(i_node), livr)
                    call tbajvr(this%table_g, this%nb_para, 'K2', this%v_K2(i_node), livr)
                    if (cgField%ndim .eq. 3) then
                        call tbajvr(this%table_g, this%nb_para, 'K3', this%v_K3(i_node), livr)
                    end if
                    call tbajvr(this%table_g, this%nb_para, 'G_IRWIN', this%v_G_IRWIN(i_node), livr)
                elseif (option == "KJ") then
                    call tbajvr(this%table_g, this%nb_para, 'KJ', this%v_KJ(i_node), livr)
                    call tbajvr(this%table_g, this%nb_para, 'G', this%v_G(i_node), livr)
                elseif (option == "KJ_EPSI") then
                    call tbajvr(this%table_g, this%nb_para, 'KJ_EPSI', this%v_KJ_EPSI(i_node), livr)
                    call tbajvr(this%table_g, this%nb_para, 'G_EPSI', this%v_G_EPSI(i_node), livr)
                elseif (option == "G_EPSI") then
                    call tbajvr(this%table_g, this%nb_para, 'G_EPSI', this%v_G_EPSI(i_node), livr)
                else
                    ASSERT(ASTER_FALSE)
                end if
            end do
!
            call tbajli(this%table_g, this%nb_para, this%list_name_para, livi, livr, livc, livk, 0)
!
        end do
!
        call cpu_time(finish)
        cgStat%save_cgTable = cgStat%save_cgTable+finish-start
        cgStat%nb_save_cgTable = cgStat%nb_save_cgTable+1
!
    end subroutine save_table
!
!
!=======================================================================================
!
!=======================================================================================
!
    subroutine create_npf(this, cgStat)
!
        implicit none
!
        class(CalcG_theta), intent(in) :: this
        type(CalcG_Stat), intent(inout)   :: cgStat
!
! --------------------------------------------------------------------------------------
!
!   Create new datastructures for nb_point_fond
!   In this     : calcG Theta
! --------------------------------------------------------------------------------------
!
        real(kind=8), pointer   :: v_basfon(:) => null()
        real(kind=8), pointer   :: basfon(:) => null()
        real(kind=8), pointer   :: coorNoeud(:) => null()
        real(kind=8), pointer   :: fondNoeudCoor(:) => null()
        real(kind=8), pointer   :: absfon(:) => null()
        real(kind=8), pointer   :: v_absfon(:) => null()
        real(kind=8), pointer   :: v_abscur(:) => null()
        real(kind=8), pointer :: gsv(:) => null()
        aster_logical, pointer :: gsl(:) => null()
        integer(kind=8), pointer   :: fondNoeudNume(:) => null()
        character(len=8), pointer   :: v_noeuf(:) => null()
        character(len=8), pointer   :: fondNoeud(:) => null()
        character(len=19) :: cnsbas
        integer(kind=8) :: i_node, nume, k, node_nume, ndim, nbno, ibid
        integer(kind=8) :: ina, inb, nseg, iseg, indica, indicb
        real(kind=8) :: smax, s, s1, s2, coor1(3), coor2(3), start, finish
        real(kind=8) :: base1(6), base2(6), eps, dmin, sn, n(3)
        real(kind=8) :: xm, ym, zm, xa, ya, za, xb, yb, zb, d
        real(kind=8) :: xab, yab, zab, xam, yam, zam, norm2, xnm, ynm, znm
        real(kind=8) :: vnora(3), vnorb(3), vdira(3), vdirb(3), vnorn(3), vdirn(3)
        character(len=8), parameter :: licmp(9) = ['X1', 'X2', 'X3', 'X4', 'X5', &
                                                   'X6', 'X7', 'X8', 'X9']
!
        call jemarq()
        call jeveuo(this%crack//'.FOND.NOEU', 'L', vk8=fondNoeud)
        call this%getCoorNodes(coorNoeud)
!
        call wkvect(this%fondNoeudNume, 'V V I', this%nnof, vi=fondNoeudNume)
        call wkvect(this%fondNoeudCoor, 'V V R', 3*this%nnof, vr=fondNoeudCoor)
!
        if (this%nb_point_fond > 0) then
            call cpu_time(start)
            call wkvect(this%npf_basefond, 'V V R', 6*this%nb_point_fond, vr=v_basfon)
            call wkvect(this%npf_fondNoeud, 'V V K8', this%nb_point_fond, vk8=v_noeuf)
            call wkvect(this%npf_absfon, 'V V R', this%nb_point_fond, vr=v_absfon)
!
            call jeveuo(this%crack//'.ABSFON', 'L', vr=absfon)
            call jeveuo(this%crack//'.BASNOF', 'L', vr=basfon)
!
! --- Create virtual nodes (id and name)
            do i_node = 2, this%nb_point_fond-1
                fondNoeudNume(i_node) = -1
                v_noeuf(i_node) = "XXXXXXXX"
            end do
!
! --- First and last nodes
!
            nume = char8_to_int(fondNoeud(1))
            fondNoeudNume(1) = nume
            v_noeuf(1) = fondNoeud(1)
            fondNoeudCoor(1:3) = coorNoeud(3*(nume-1)+1:3*(nume-1)+3)
            v_absfon(1) = 0.d0
            v_basfon(1:6) = basfon(1:6)
!
            nume = char8_to_int(fondNoeud(this%nb_fondNoeud))
            fondNoeudNume(this%nb_point_fond) = nume
            v_noeuf(this%nb_point_fond) = fondNoeud(this%nb_fondNoeud)
            fondNoeudCoor(3*(this%nb_point_fond-1)+1:3*(this%nb_point_fond-1)+3) = &
                coorNoeud(3*(nume-1)+1:3*(nume-1)+3)
            v_absfon(this%nb_point_fond) = absfon(this%nb_fondNoeud)
            v_basfon(6*(this%nb_point_fond-1)+1:6*(this%nb_point_fond-1)+6) = &
                basfon(6*(this%nb_fondNoeud-1)+1:6*(this%nb_fondNoeud-1)+6)
!
! --- Create virtual nodes (coordinates)
            smax = absfon(this%nb_fondNoeud)
            do i_node = 2, this%nb_point_fond-1
                s = (i_node-1)*smax/(this%nb_point_fond-1)
                do k = 1, this%nb_fondNoeud
                    if (absfon(k) .gt. s) then
                        node_nume = k
                        exit
                    end if
                end do
!         ON INTERPOLE LES COORD ENTRE CELLES DU SEGMENT [K-1,K]
                s1 = absfon(node_nume-1)
                nume = char8_to_int(fondNoeud(node_nume-1))
                coor1(1:3) = coorNoeud(3*(nume-1)+1:3*(nume-1)+3)
                base1(1:6) = basfon(6*(node_nume-1-1)+1:6*(node_nume-1-1)+6)
                s2 = absfon(node_nume)
                nume = char8_to_int(fondNoeud(node_nume))
                coor2(1:3) = coorNoeud(3*(nume-1)+1:3*(nume-1)+3)
                base2(1:6) = basfon(6*(node_nume-1)+1:6*(node_nume-1)+6)
!
                fondNoeudCoor(3*(i_node-1)+1:3*(i_node-1)+3) = &
                    coor1(1:3)+(coor2(1:3)-coor1(1:3))*(s-s1)/(s2-s1)
                v_basfon(6*(i_node-1)+1:6*(i_node-1)+6) = &
                    base1(1:6)+(base2(1:6)-base1(1:6))*(s-s1)/(s2-s1)
                v_absfon(i_node) = s
            end do

            do i_node = 1, this%nb_point_fond
                nume = char8_to_int(fondNoeud(i_node))
            end do
!
! --- Create Local Basis
            call dismoi('DIM_GEOM', this%mesh, 'MAILLAGE', repi=ndim)
            ASSERT(ndim == 3)
            call dismoi('NB_NO_MAILLA', this%mesh, 'MAILLAGE', repi=nbno)
            nseg = this%nb_point_fond-1
            call wkvect(this%npf_abscurv, 'V V R', nbno, vr=v_abscur)
!
!     INITIALISATION DU CHAMP SIMPLE DE LA BASE LOCALE
            cnsbas = '&&CALC_G.CNSBAS'
            call cnscre(this%mesh, 'NEUT_R', ndim*3, licmp, 'V', &
                        cnsbas)
            call jeveuo(cnsbas//'.CNSV', 'E', vr=gsv)
            call jeveuo(cnsbas//'.CNSL', 'E', vl=gsl)
!
            eps = 1.d-12
            do i_node = 1, nbno
!
!       COORD DU NOEUD M DU MAILLAGE
                xm = coorNoeud((i_node-1)*3+1)
                ym = coorNoeud((i_node-1)*3+2)
                zm = coorNoeud((i_node-1)*3+3)
!
!         RECHERCHE DU PROJETE DE INO SUR LE FOND DE FISSURE
!         --------------------------------------------------
                dmin = r8maem()
!
!         BOUCLE SUR LES "SEGMENTS" DU FOND DE FISSURE
                do iseg = 1, nseg
                    ina = iseg
                    inb = iseg+1
!
!           COORD DES POINTS A ET B, EXTREMITES DU SEGMENT ISEG
                    xa = fondNoeudCoor((ina-1)*3+1)
                    ya = fondNoeudCoor((ina-1)*3+2)
                    za = fondNoeudCoor((ina-1)*3+3)
                    xb = fondNoeudCoor((inb-1)*3+1)
                    yb = fondNoeudCoor((inb-1)*3+2)
                    zb = fondNoeudCoor((inb-1)*3+3)
!
!           VECTEUR AB ET AM
                    xab = xb-xa
                    yab = yb-ya
                    zab = zb-za
                    xam = xm-xa
                    yam = ym-ya
                    zam = zm-za
!
!           PARAM S (PRODUIT SCALAIRE...)
                    s = xab*xam+yab*yam+zab*zam
                    norm2 = xab*xab+yab*yab+zab*zab
                    s = s/norm2
!
!           SI N EN DEHORS DU SEGMENT AB
                    if ((s-1) .ge. eps) s = 1.d0
                    if (s .le. eps) s = 0.d0
!
!           COORD DU PROJETE DE M SUR ISEG: N
                    xnm = xm-(s*xab+xa)
                    ynm = ym-(s*yab+ya)
                    znm = zm-(s*zab+za)
!
!           DISTANCE MN
                    d = sqrt(xnm*xnm+ynm*ynm+znm*znm)
!
                    if (d .lt. (dmin*(1-abs(r8prem())*1.d04))) then
                        dmin = d
                        sn = s
                        indica = ina
                        indicb = inb
!
                        n(1) = s*xab+xa
                        n(2) = s*yab+ya
                        n(3) = s*zab+za
                    end if
                end do
!           ABSCISSE CURVILIGNE DU NOEUD N SUR LE FRONT DE FISSURE
                v_abscur(i_node) = (1-sn)*v_absfon(indica)+sn*v_absfon(indicb)
!
!         CALCUL DES VECTEURS DE LA BASE LOCALE AU POINT PROJETE
!
                do k = 1, ndim
                    vnora(k) = v_basfon(6*(indica-1)+k)
                    vdira(k) = v_basfon(6*(indica-1)+k+ndim)
                    vnorb(k) = v_basfon(6*(indicb-1)+k)
                    vdirb(k) = v_basfon(6*(indicb-1)+k+ndim)
                    vnorn(k) = sn*vnorb(k)+(1-sn)*vnora(k)
                    vdirn(k) = sn*vdirb(k)+(1-sn)*vdira(k)
!
!               STOCKAGE COORDONNEES DU PROJETE
                    gsv(3*ndim*(i_node-1)+k) = n(k)
                    gsl(3*ndim*(i_node-1)+k) = ASTER_TRUE
!               STOCKAGE VECTEUR DIRECTION
                    gsv(3*ndim*(i_node-1)+k+3) = vdirn(k)
                    gsl(3*ndim*(i_node-1)+k+3) = ASTER_TRUE
!               STOCKAGE VECTEUR NORMAL
                    gsv(3*ndim*(i_node-1)+k+6) = vnorn(k)
                    gsl(3*ndim*(i_node-1)+k+6) = ASTER_TRUE

                end do
            end do
            call cnscno(cnsbas, ' ', 'NON', 'G', this%npf_baseloc, &
                        'F', ibid)
            call cpu_time(finish)
            cgStat%npf_cgTheta = cgStat%npf_cgTheta+finish-start
        else
            do i_node = 1, this%nb_fondNoeud
!               Récupération du numéro de noeud
                nume = char8_to_int(fondNoeud(i_node))
                fondNoeudNume(i_node) = nume
                fondNoeudCoor(3*(i_node-1)+1:3*(i_node-1)+3) = coorNoeud(3*(nume-1)+1:3*(nume-1)+3)
            end do
        end if
        call jedema()
!
    end subroutine
!
end module
