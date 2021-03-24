! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcG_type.h"
#include "asterfort/ccbcop.h"
#include "asterfort/cgcrio.h"
#include "asterfort/cgReadCompor.h"
#include "asterfort/cgComporNodes.h"
#include "asterfort/comp_info.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcncon.h"
#include "asterfort/gettco.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/ismali.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
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
#include "jeveux.h"
!
public :: CalcG_Field, CalcG_Study, CalcG_Theta, CalcG_Table
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
        integer            :: ndim      = 0
! ----- name of list of NUME
        character(len=24)  :: list_nume_name = ' '
! ----- list of nume
        integer, pointer   :: list_nume(:) => null()
! ----- number of order
        integer            :: nb_nume = 0
! ----- number option to compute
        integer            :: nb_option = 0
! ----- list of options
        character(len=8)   :: list_option(NB_MAX_OPT) = ' '
! ----- level information
        integer            :: level_info = 1
        contains
        procedure, pass    :: initialize => initialize_field
        procedure, pass    :: print      => print_field
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
        character(len=8)   :: model     = ' '
! ----- name of mesh
        character(len=8)   :: mesh      = ' '
! ----- name of material
        character(len=24)  :: material  = ' '
! ----- name of coded material
        character(len=24)  :: mateco    = ' '
! ----- name of elementary characteristics
        character(len=24)  :: carael    = ' '
! ----- name of loading
        character(len=24)  :: loading   = ' '
! ----- index order
        integer            :: nume_ordre = -1
! ----- option to compute
        character(len=8)   :: option    = ' '
!------ linear or quadratic
        aster_logical      :: milieu = ASTER_FALSE
!------ modal analysis ?
        aster_logical      :: l_modal = ASTER_FALSE
!------ axisymetric model
        aster_logical      :: l_axis = ASTER_FALSE
! ----- displacement field
        character(len=24)  :: depl   = ' '
! ----- speed field
        character(len=24)  :: vitesse   = ' '
! ----- acceleration field
        character(len=24)  :: acce   = ' '
! ----- time
        real(kind=8)       :: time   = 0.d0
! ----- pulse
        real(kind=8)       :: pulse   = 0.d0
! ----- computed values (G, K1, K2, K3, FIC1, FIC2, FIC3)
        real(kind=8)       :: gth(NB_MAX_TERM)   = 0.d0
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
! ----- name of factors necessary to create theta field in te
        character(len=24)       :: theta_factors = ' '
! ----- name of the matrix from A*G(s)=g(theta)
        character(len=24)       :: matrix = ' '
! ----- number of theta field
        integer                 :: nb_theta_field = 0
! ----- name of crack
        character(len=8)        :: crack = ' '
! ----- type of crack
        character(len=24)       :: crack_type = ' '
! ----- initial configuration of the crack
        character(len=8)        :: config_init = ' '
! ----- name of the mesh (support of crack)
        character(len=8)        :: mesh = ' '
! ----- name of the nodes of the mesh (support of crack)
        character(len=24)       :: nomNoeud = ' '
! ----- number of nodes in the crack
        integer                 :: nb_fondNoeud = 0
! ----- rayon
        real(kind=8)            :: r_inf = 0.d0, r_sup = 0.d0
! ----- lenght of crack
        real(kind=8)            :: lonfis = 0.d0
! ----- number of layer
        integer                 :: nb_couche_inf = 0, nb_couche_sup = 0
! ----- type of discretization
        character(len=8)        :: discretization = ' '
! ----- nubmer of nodes if nb_pts_fond is defined (for linear discretization)
        integer                 :: nb_point_fond = 0
! ----- number of nodes (for linear discretization)
        integer                 :: nnof = 0
! ----- nubmer of nodes (for legendre discretization)
        integer                 :: degree = 0
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
! ----- member function
        contains
        procedure, pass    :: initialize => initialize_theta
        procedure, pass    :: print => print_theta
        procedure, pass    :: compute_curvature
        procedure, pass    :: getCoorNodes
        procedure, pass    :: getAbscurv
        procedure, pass    :: getAbsfon
        procedure, pass    :: getBaseLoc
        procedure, pass    :: getFondTailleR
        procedure, pass    :: getFondNoeu
        procedure, pass    :: getFondNoeuNume
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
        integer            :: nb_para = 0
! ----- grandeur
        real(kind=8), pointer :: v_G(:) => null()
        real(kind=8), pointer :: v_K1(:) => null()
        real(kind=8), pointer :: v_K2(:) => null()
        real(kind=8), pointer :: v_K3(:) => null()
        real(kind=8), pointer :: v_G_IRWIN(:) => null()
        real(kind=8), pointer :: v_G_EPSI(:) => null()
        real(kind=8), pointer :: v_TEMP(:) => null()
        character(len=8), pointer :: v_COMPOR(:) => null()
! ----- nb point
        integer :: nb_point = 1
!
! ----- member function
        contains
        procedure, pass    :: initialize => initialize_table
        procedure, pass    :: addValues
        procedure, pass    :: addPara
        procedure, pass    :: save => save_table
!
    end type CalcG_Table
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
    subroutine initialize_field(this)
!
    implicit none
!
        class(CalcG_Field), intent(inout)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Field type
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
        integer :: ier, ifm, jopt, nbropt
        character(len=3) :: repk
        character(len=16) :: k16bid
        character(len=19) :: lisopt
        character(len=8)  :: modele, mater
        integer, pointer :: v_nume(:) => null()
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
        ASSERT(ier==1)
        call gettco(this%result_in, this%result_in_type, ASTER_TRUE)
        call dismoi('MODELE', this%result_in, 'RESULTAT', repk=modele)
        call dismoi('DIM_GEOM', modele, 'MODELE', repi=this%ndim)
!
        call dismoi('CHAM_MATER', this%result_in, 'RESULTAT', repk=mater)
        call dismoi('EXI_VARC', mater, 'CHAM_MATER', repk=repk)
        if( repk == "OUI" ) then
            this%l_temp = ASTER_TRUE
        end if
!
! --- Get name of option
!
        call getvtx(' ', 'OPTION', nbret=ier)
        if(ier == 1) then
            this%nb_option = ier
            call getvtx(' ', 'OPTION', scal=this%list_option(1))
        else
            this%nb_option = -ier
            ASSERT(this%nb_option <= NB_MAX_OPT)
            call getvtx(' ', 'OPTION', nbval=this%nb_option, vect=this%list_option)
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
        if(this%level_info > 1) then
            call comp_info(modele, this%compor)
        end if
!
! --- if ELAS_INCR
!
        if(this%l_incr) then
            lisopt = '&&OP0027.LISOPT'
            nbropt = 2
            call wkvect(lisopt, 'V V K16', nbropt, jopt)
            zk16(jopt) = 'VARI_ELNO'
            zk16(jopt+1) = 'EPSP_ELNO'
            call ccbcop(this%result_in, this%result_out, this%list_nume_name,&
                        this%nb_nume, lisopt, nbropt)
        end if
!
        call jedema()
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
        integer :: i
!
        print*, "----------------------------------------------------------------------"
        print*, "Informations about CalcG_Field"
        print*, "Level of informations: ", this%level_info
        print*, "Dimension of the problem: ", this%ndim
        print*, "Result: ", this%result_in, " of type ", this%result_in_type
        print*, "Output: ", this%table_out
        print*, "Number of option to compute: ", this%nb_option
        do i =1, this%nb_option
            print*, "** option ", i, ": ", this%list_option(i)
        end do
        print*, "Number of step/mode to compute: ", this%nb_nume
        print*, "----------------------------------------------------------------------"
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
        if(this%result_in_type == "MODE_MECA") then
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
        if(this%result_in_type == "DYNA_TRANS") then
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
    subroutine clean_field(this)
!
    implicit none
!
        class(CalcG_Field), intent(inout)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   Clean objects
!   In this     : calG field
! --------------------------------------------------------------------------------------------------
!
    integer :: iret
    integer, pointer :: ordr(:) => null()
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
        endif
!
        call jedetr(this%list_nume_name)
        call jedetr(this%result_out)
        call rsmena(this%result_in)
    endif
!
    AS_DEALLOCATE(vi=this%list_nume)
!
    call jedema()
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_study(this, result_in, nume_index)
!
    implicit none
!
        class(CalcG_Study), intent(inout)  :: this
        character(len=8), intent(in)       :: result_in
        integer, intent(in)                :: nume_index
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Study type
!   In this     : study type
!   In nume_index : index nume
! --------------------------------------------------------------------------------------------------
!
        character(len=8) :: typma, typmo
        integer :: jma

        call jemarq()

        this%nume_ordre = nume_index
        this%loading    = "&&STUDY.CHARGES"
        call medomg(result_in, this%nume_ordre, this%model, this%material, &
                    this%mateco, this%loading)
        call dismoi('CARA_ELEM', result_in, 'RESULTAT', repk=this%carael)
        call dismoi('NOM_MAILLA', this%model,'MODELE', repk=this%mesh)
        call dismoi('MODELISATION', this%model, 'MODELE', repk=typmo)
        if( typmo(1:4) == "AXIS") then
            this%l_axis = ASTER_TRUE
        else
            this%l_axis = ASTER_FALSE
        end if
!
!       Maillage linéaire ou quadratique
        call jeveuo(this%mesh//'.TYPMAIL', 'L', jma)
        call jenuno(jexnum('&CATA.TM.NOMTM', zi(jma)), typma)
        if (.not. ismali(typma)) then
            this%milieu = ASTER_TRUE
        endif

        call jedema()

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
        print*, "----------------------------------------------------------------------"
        print*, "Informations about CalcG_Study"
        print*, "Option: ", this%option
        print*, "Model: ", this%model
        print*, "Mesh: ", this%mesh
        print*, "Material: ", this%material
        print*, "Coded material: ", this%mateco
        print*, "Loading: ", this%loading
        print*, "Nume index: ", this%nume_ordre
        PRINT*, "linear or quadratic", this%milieu
        print*, "----------------------------------------------------------------------"
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
            endif
        endif
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
        integer :: iret
!
        call rsexch('F', result_in, 'DEPL', this%nume_ordre, this%depl, iret)
        call rsexch(' ', result_in, 'VITE', this%nume_ordre, this%vitesse, iret)
        if (iret .ne. 0) then
            this%vitesse = ' '
            this%acce = ' '
        else
            call rsexch(' ', result_in, 'ACCE', this%nume_ordre, this%acce, iret)
        endif
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
        integer :: ipuls, jinst
        character(len=8)  :: k8bid

!
        if (this%l_modal) then
            call rsadpa(result_in, 'L', 1, 'OMEGA2', this%nume_ordre,&
                        0, sjv=ipuls, styp=k8bid)
            this%pulse = sqrt(zr(ipuls))
            this%time = 0.d0
        else
            call rsadpa(result_in, 'L', 1, 'INST', this%nume_ordre,&
                        0, sjv=jinst, styp=k8bid)
            this%pulse = 0.d0
            this%time = zr(jinst)
        endif
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine initialize_theta(this)
!
    implicit none
!
        class(CalcG_Theta), intent(inout)  :: this
!
! --------------------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Theta type
!   In this     : theta type
! --------------------------------------------------------------------------------------------------
!
        integer :: ier, j, ibasf, i, num, inume
        character(len=8) :: typfon
        real(kind=8) :: maxtai, mintai
        aster_logical :: l_disc
        real(kind=8), pointer :: fondTailleR(:) => null()
        real(kind=8), pointer :: absfon(:)  => null()
        character(len=8), pointer :: fondNoeud(:)  => null()
!
        call jemarq()
! --- get automatic name
        call gcncon("_", this%theta_field)
        call gcncon("_", this%theta_factors)
        call gcncon("_", this%matrix)
        call gcncon("_", this%absfond)
        call gcncon("_", this%fondNoeudNume)
!
! --- get informations about the crack
!
        call getvtx('THETA', 'FISSURE', iocc=1, scal=this%crack, nbret=ier)
        ASSERT(ier==1)
        call gettco(this%crack, this%crack_type, ASTER_TRUE)
!
        call dismoi('NOM_MAILLA',this%crack,'FOND_FISS', arret='F', repk=this%mesh)
        call dismoi('TYPE_FOND', this%crack,'FOND_FISS', arret='F', repk=typfon)
        call dismoi('CONFIG_INIT', this%crack, 'FOND_FISS', repk=this%config_init)
        call dismoi('SYME', this%crack, 'FOND_FISS', repk=this%symech)
! --- the crack is closed ?
        if (typfon .eq. 'FERME') then
            this%l_closed = .true.
        else
            this%l_closed = .false.
        endif
! --- number of nodes in the crack
        call jelira(this%crack//'.FOND.NOEU', 'LONMAX', this%nb_fondNoeud)
!
! --- get informations about theta discretization
!
        call getvtx('THETA', 'DISCRETISATION', iocc=1, scal=this%discretization, nbret=ier)
        l_disc = (this%discretization == "LINEAIRE") .or. (this%discretization == "LEGENDRE")
        ASSERT(l_disc)
!
        call getvis('THETA', 'DEGRE', iocc=1, scal=this%degree, nbret=ier)
        call jelira(this%crack//'.FOND.NOEU', 'LONMAX', this%nnof)

!
!        if ( this%discretization == "LINEAIRE") then
!            this%nb_theta_field = this%nnof
!        else
!           this%nb_theta_field = this%degree
!        endif

!        if(this%discretization == "LINEAIRE") then
!            ASSERT(this%nb_point_fond >= 0)
!            ASSERT(this%degree == 0)
!        end if
!
        if(this%discretization == "LEGENDRE") then
!            ASSERT(this%nb_point_fond == 0)
!            ASSERT(this%degree >= 0)
            if(this%l_closed) then
                call utmess('F', 'RUPTURE0_90')
            end if
        end if
!
        call getvis('THETA', 'NB_COUCHE_INF', iocc=1, scal=this%nb_couche_inf, nbret=ier)
        call getvis('THETA', 'NB_COUCHE_SUP', iocc=1, scal=this%nb_couche_sup, nbret=ier)
!
        if(ier == 1 .and. ((this%nb_couche_inf < 0) .or. &
            (this%nb_couche_inf >= this%nb_couche_sup))) then
            call utmess('F', 'RUPTURE3_4', ni=2, vali=[this%nb_couche_inf, this%nb_couche_sup])
        end if
!
        this%nomNoeud = this%mesh//'.NOMNOE'

! --- Get RINF and DE RSUP from command file or from SD FOND_FISSURE

        call getvr8('THETA', 'R_INF', iocc=1, scal=this%r_inf, nbret=ier)
        call getvr8('THETA', 'R_SUP', iocc=1, scal=this%r_sup, nbret=ier)
!
        if(ier == 1 .and. ((this%r_inf < 0.d0) .or. (this%r_inf >= this%r_sup))) then
            call utmess('F', 'RUPTURE3_3', nr=2, valr=[this%r_inf, this%r_sup])
        end if

        if (ier .eq. 0) then
            if (this%config_init .eq. 'DECOLLEE') then
                call utmess('F', 'RUPTURE1_7')
            endif
            call this%getFondTailleR(fondTailleR)
            maxtai = fondTailleR(1)
            mintai = fondTailleR(1)
            do j = 1, this%nb_fondNoeud
                maxtai = max(maxtai,fondTailleR(j))
                mintai = min(mintai,fondTailleR(j))
            end do
            this%r_inf = 2*maxtai
            this%r_sup = 4*maxtai
            call utmess('I', 'RUPTURE1_5', nr=2, valr=[this%r_inf, this%r_sup])
            if (maxtai .gt. 2*mintai) then
                call utmess('A', 'RUPTURE1_16', nr=2, valr=[mintai, maxtai])
            endif
        endif
!
!       Extraction à modifier lors de la résolution de issue30288 (NB_POINT_FOND)
        call wkvect(this%absfond, 'V V R8', this%nb_fondNoeud, ibasf)
        call wkvect(this%fondNoeudNume, 'V V I', this%nb_fondNoeud, inume)
        call this%getAbscurv(absfon)
        call this%getFondNoeu(fondNoeud)
!
        do i = 1, this%nb_fondNoeud
!
!           Récupération du numéro de noeud
            call jenonu(jexnom(this%nomNoeud, fondNoeud(i)), num)
!
            zr(ibasf-1+i) = absfon(i)
            zi(inume-1+i) = num
        enddo
!
        call jedema()
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
        call jeveuo(this%crack//'.ABSCUR', 'L', vr=v_abs)
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
        call jeveuo(this%crack//'.ABSFON', 'L', vr=v_absfon)
        call jedema()
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
        call jeveuo(this%crack//'.BASLOC    .VALE', 'L', vr=v_base)
        call jedema()
!
    end subroutine
!===================================================================================================
!
    subroutine getFondNoeu(this, v_fondnoeu)
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
        call jeveuo(this%crack//'.FOND.NOEU', 'L', vk8=v_fondnoeu)
        call jedema()
!
    end subroutine
!===================================================================================================
!
    subroutine getFondNoeuNume(this, v_fondnoeuNume)
!
    implicit none
!
        class(CalcG_Theta), intent(in)  :: this
        integer, pointer :: v_fondnoeuNume(:)
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
        print*, "----------------------------------------------------------------------"
        print*, "Informations about CalcG_Theta"
        print*, "Field theta: ", this%theta_field
        print*, "theta_factors: ", this%theta_factors
        print*, "Crack: ", this%crack, " of type ", this%crack_type
        print*, "Number of nodes in the crack: ", this%nb_fondNoeud
        print*, "Mesh support: ", this%mesh
        print*, "Initial configuration: ", this%config_init
        print*, "the crack is symetric: ", this%symech
        print*, "The crack is closed ?: ", this%l_closed
!        print*, "Nombre de champs THETA: ", this%nb_theta_field
        print*, "Discretization : ", this%discretization,  " with number/degree ", &
                this%nnof, this%degree
        print*, "Radius:"
        print*, "*** Inferior: ", this%r_inf
        print*, "*** Superior: ", this%r_sup
        print*, "Number of cell layers:"
        print*, "*** Inferior: ", this%nb_couche_inf
        print*, "*** Superior: ", this%nb_couche_sup
        print*, "----------------------------------------------------------------------"
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
        this%nb_para = this%nb_para + 1
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
    subroutine initialize_table(this, cgField, cgTheta)
!
    implicit none
!
        class(CalcG_Table), intent(inout)  :: this
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
!
! --------------------------------------------------------------------------------------
!
!   initialization of a CalcG_Table type
!   In this     : calcG Table
! --------------------------------------------------------------------------------------
        integer :: iopt, nbValues
        character(len=8) :: option
        integer, pointer :: fondNoeudNume(:) => null()

!
! --- Table pour les valeurs (table)
!
        call gcncon("_", this%table_g)
        call tbcrsd(this%table_g, 'G')
!
        this%nb_point = cgTheta%nb_fondNoeud
        this%nb_para = 0
!
! --- INST or FREQ
        if (cgField%isModeMeca()) then
            call this%addPara('NUME_MODE', 'I')
        else
            call this%addPara('NUME_ORDRE', 'I')
            call this%addPara('INST', 'R')
        endif
! --- Node name
        call this%addPara('NOEUD', 'K8')
        call this%addPara('NUM_PT', 'I')
! --- Coordinates of nodes
        call this%addPara('COOR_X', 'R')
        call this%addPara('COOR_Y', 'R')
        if (cgField%ndim.eq.3) then
            call this%addPara('COOR_Z', 'R')
            call this%addPara('ABSC_CURV', 'R')
            call this%addPara('ABSC_CURV_NORM', 'R')
        endif
! --- Tempature
        if( cgField%l_temp ) then
            call this%addPara('TEMP', 'R')
            call wkvect("&&TABLEG.TEMP", 'V V R', this%nb_point, vr=this%v_TEMP)
        end if
! --- Behavior
        call this%addPara('COMPORTEMENT', 'K8')
        call cgTheta%getFondNoeuNume(fondNoeudNume)
        call cgComporNodes(cgField%result_in, cgField%list_nume(1), this%nb_point, &
                            fondNoeudNume, this%v_COMPOR)
! --- Option
        nbValues = this%nb_point
        do iopt = 1, cgField%nb_option
            option = cgField%list_option(iopt)

            if (option == "G" ) then
                call this%addPara('G', 'R')
                call wkvect("&&TABLEG.G", 'V V R', nbValues, vr=this%v_G)
            elseif (option == "K" ) then
                call this%addPara('K1', 'R')
                call wkvect("&&TABLEG.K1", 'V V R', nbValues, vr=this%v_K1)
                call this%addPara('K2', 'R')
                call wkvect("&&TABLEG.K2", 'V V R', nbValues, vr=this%v_K2)
                if (cgField%ndim.eq.3) then
                    call this%addPara('K3', 'R')
                    call wkvect("&&TABLEG.K3", 'V V R', nbValues, vr=this%v_K3)
                endif
                call this%addPara('G_IRWIN', 'R')
                call wkvect("&&TABLEG.GIR", 'V V R', nbValues, vr=this%v_G_IRWIN)
            elseif (option == "G_EPSI" ) then
                call this%addPara('G_EPSI', 'R')
                call wkvect("&&TABLEG.GEP", 'V V R', nbValues, vr=this%v_G_EPSI)
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
!
! --- create table
        call tbajpa(this%table_g, this%nb_para, this%list_name_para, this%list_type_para)
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
        integer, intent(in) :: node_id
!
! ----------------------------------------------------------------------------------------
!
!   add Values in the table
!   In this     : calcG Table
! ----------------------------------------------------------------------------------------
!
        if(cgStudy%option == "G") then
            this%v_G(node_id) = cgStudy%gth(1)
        elseif(cgStudy%option == "K") then
            this%v_K1(node_id) = cgStudy%gth(2)
            this%v_K2(node_id) = cgStudy%gth(3)
            if(cgField%ndim == 3) then
                this%v_K3(node_id) = cgStudy%gth(4)
            end if
            this%v_G_IRWIN(node_id) = cgStudy%gth(5)**2 + cgStudy%gth(6)**2 + cgStudy%gth(7)**2
        elseif(cgStudy%option == "G_EPSI") then
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
    subroutine save_table(this, cgField, cgTheta, cgStudy)
!
    implicit none
!
        class(CalcG_Table), intent(in)  :: this
        type(CalcG_field), intent(in) :: cgField
        type(CalcG_theta), intent(in) :: cgTheta
        type(CalcG_study), intent(in) :: cgStudy
!
! --------------------------------------------------------------------------------------
!
!   save values in the table
!   In this     : calcG Table
! --------------------------------------------------------------------------------------
!
        integer            :: livi(NB_MAX_PARA)
        real(kind=8)       :: livr(NB_MAX_PARA)
        complex(kind=8)    :: livc(NB_MAX_PARA)
        character(len=24)  :: livk(NB_MAX_PARA)
        real(kind=8) :: coor(3)
        integer :: i_node, iopt, node_id
        character(len=8) :: option
        real(kind=8), pointer   :: coorNoeud(:) => null()
        real(kind=8), pointer   :: abscur(:) => null()
        integer, pointer   :: fondNoeudNume(:) => null()
        character(len=8), pointer :: fondNoeud(:) => null()
!
        if (cgField%isModeMeca()) then
            call tbajvi(this%table_g, this%nb_para, 'NUME_MODE', cgStudy%nume_ordre, livi)
        else
            call tbajvi(this%table_g, this%nb_para, 'NUME_ORDRE', cgStudy%nume_ordre, livi)
            call tbajvr(this%table_g, this%nb_para, 'INST', cgStudy%time, livr)
        endif
!
        call cgTheta%getCoorNodes(coorNoeud)
        call cgTheta%getAbscurv(abscur)
        call cgTheta%getFondNoeuNume(fondNoeudNume)
        call cgTheta%getFondNoeu(fondNoeud)
!
        do i_node = 1, this%nb_point
            node_id = fondNoeudNume(i_node)
            call tbajvk(this%table_g, this%nb_para, 'NOEUD', fondNoeud(i_node), livk)
            call tbajvi(this%table_g, this%nb_para, 'NUM_PT', i_node, livi)
!
            coor = coorNoeud((node_id-1)*3+1:(node_id-1)*3+3)
            call tbajvr(this%table_g, this%nb_para, 'COOR_X', coor(1), livr)
            call tbajvr(this%table_g, this%nb_para, 'COOR_Y', coor(2), livr)
            if (cgField%ndim.eq.3) then
                call tbajvr(this%table_g, this%nb_para, 'COOR_Z', coor(3), livr)
                call tbajvr(this%table_g, this%nb_para, 'ABSC_CURV', abscur(node_id), livr)
                call tbajvr(this%table_g, this%nb_para, 'ABSC_CURV_NORM', &
                            abscur(node_id)/cgTheta%lonfis, livr)
            endif
!
            if( cgField%l_temp ) then
                call tbajvr(this%table_g, this%nb_para, 'TEMP', this%v_TEMP(i_node), livr)
            end if
            call tbajvk(this%table_g, this%nb_para, 'COMPORTEMENT', this%v_COMPOR(i_node), livk)
!
            do iopt = 1, cgField%nb_option
                option = cgField%list_option(iopt)
                if(option == "G") then
                    call tbajvr(this%table_g, this%nb_para, 'G', this%v_G(i_node), livr)
                elseif (option == "K") then
                    call tbajvr(this%table_g, this%nb_para, 'K1', this%v_K1(i_node), livr)
                    call tbajvr(this%table_g, this%nb_para, 'K2', this%v_K2(i_node), livr)
                    if (cgField%ndim .eq. 3) then
                        call tbajvr(this%table_g, this%nb_para, 'K3', this%v_K3(i_node), livr)
                    endif
                    call tbajvr(this%table_g, this%nb_para, 'G_IRWIN', this%v_G_IRWIN(i_node), livr)
                elseif(option == "G_EPSI") then
                    call tbajvr(this%table_g, this%nb_para, 'G_EPSI', this%v_G_EPSI(i_node), livr)
                else
                    ASSERT(ASTER_TRUE)
                end if
            end do
!
            call tbajli(this%table_g, this%nb_para, this%list_name_para, livi, livr, livc, livk, 0)
!
        end do
!
    end subroutine save_table
!
end module
