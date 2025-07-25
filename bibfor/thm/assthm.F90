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
! aslint: disable=W1504,W1306
!
subroutine assthm(ds_thm, option, j_mater, &
                  lMatr, lSigm, lVect, &
                  lVari, lMatrPred, l_axi, &
                  typmod, inte_type, angl_naut, &
                  ndim, nbvari, nno, nnos, &
                  npg, npi, &
                  nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                  dimdef, dimcon, dimuel, &
                  mecani, press1, press2, tempe, second, &
                  compor, carcri, &
                  jv_poids, jv_poids2, &
                  jv_func, jv_func2, &
                  jv_dfunc, jv_dfunc2, &
                  elem_coor, &
                  dispm, dispp, &
                  congem, congep, &
                  vintm, vintp, &
                  time_prev, time_curr, &
                  matuu, vectu, codret)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cabthm.h"
#include "asterfort/equthm.h"
#include "asterfort/pmathm.h"
#include "asterfort/thmGetBehaviour.h"
#include "asterfort/thmGetBehaviourVari.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/thmGetParaInit.h"
#include "asterfort/thmSelectMatrix.h"
#include "asterfort/thmGetBehaviourChck.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    integer(kind=8), parameter :: dimmat = 120
    character(len=16), intent(in) :: option
    aster_logical, intent(in) :: lMatr, lSigm, lVari, lMatrPred, lVect
    integer(kind=8), intent(in) :: j_mater
    aster_logical, intent(in)  :: l_axi
    character(len=8), intent(in) :: typmod(2)
    character(len=3), intent(in) :: inte_type
    real(kind=8), intent(in)  :: angl_naut(3)
    integer(kind=8), intent(in) :: nbvari, ndim
    integer(kind=8), intent(in) :: nno, nnos
    integer(kind=8), intent(in) :: npg, npi
    integer(kind=8), intent(in) :: nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd
    integer(kind=8), intent(in) :: dimuel, dimdef, dimcon
    integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
    character(len=16), intent(in)  :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    integer(kind=8), intent(in) :: jv_poids, jv_poids2
    integer(kind=8), intent(in) :: jv_func, jv_func2
    integer(kind=8), intent(in) :: jv_dfunc, jv_dfunc2
    real(kind=8), intent(in) :: elem_coor(ndim, nno)
    real(kind=8), intent(in) :: dispm(dimuel), dispp(dimuel)
    real(kind=8), intent(inout) :: congem(dimcon*npi)
    real(kind=8), intent(inout) :: congep(dimcon*npi)
    real(kind=8), intent(in) :: vintm(nbvari*npi)
    real(kind=8), intent(inout) :: vintp(nbvari*npi)
    real(kind=8), intent(in) :: time_prev, time_curr
    real(kind=8), intent(inout) :: matuu(dimuel*dimuel)
    real(kind=8), intent(inout) :: vectu(dimuel)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! THM - Compute
!
! Compute non-linear options - General assembling for all physics
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! In  option           : name of option to compute
! In  j_mater          : coded material address
! In  l_axi            : flag is axisymmetric model
! In  l_steady         : flag for no-transient problem
! In  typmod           : type of modelization (TYPMOD2)
! In  inte_type        : type of integration - classical, lumped (D), reduced (R)
! In  angl_naut        : nautical angles
! In  ndim             : dimension of space (2 or 3)
! In  nbvari           : total number of internal state variables
! In  nno              : total number of nodes
! In  nnos             : number of nodes (not middle ones)
! In  npg              : number of Gauss points
! In  npi              : number of Gauss points for linear
! In  nddls            : number of dof at nodes (not middle ones)
! In  nddlm            : number of dof at nodes (middle ones)
! In  nddl_meca        : number of dof for mechanical quantity
! In  nddl_p1          : number of dof for first hydraulic quantity
! In  nddl_p2          : number of dof for second hydraulic quantity
! In  nddl_2nd         : number of dof for second gradient
! In  dimdef           : dimension of generalized strains vector
! In  dimcon           : dimension of generalized stresses vector
! In  dimuel           : number of dof for element
! In  mecani           : parameters for mechanic
! In  press1           : parameters for hydraulic (capillary pressure)
! In  press2           : parameters for hydraulic (gaz pressure)
! In  tempe            : parameters for thermic
! In  second           : parameters for second gradient
! In  compor           : behaviour
! In  carcri           : parameters for comportment
! In  jv_poids         : JEVEUX adress for weight of Gauss points (linear)
! In  jv_poids2        : JEVEUX adress for weight of Gauss points (quadratic)
! In  jv_func          : JEVEUX adress for shape functions (linear)
! In  jv_func2         : JEVEUX adress for shape functions (quadratic)
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions (linear)
! In  jv_dfunc2        : JEVEUX adress for derivative of shape functions (quadratic)
! In  elem_coor        : coordinates of node for current element
! In  dispm            : displacements - At begin of current step
! In  dispp            : displacements - At end of current step
! IO  congem           : generalized stresses - At begin of current step
! IO  congep           : generalized stresses - At end of current step
! In  vintm            : internal state variables - At begin of current step
! IO  vintp            : internal state variables - At end of current step
! In  time_prev        : time at beginning of step
! In  time_curr        : time at end of step
! IO  matuu            : tangent matrix
! IO  vectu            : non-linear forces
! Out codret           : return code for error
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: time_incr, parm_theta
    integer(kind=8) :: kpi, ipi
    integer(kind=8) :: i, j, n, k, kji
    integer(kind=8) :: nb_vari_meca
    real(kind=8) :: a(2), as(2), ak(2), poids, poids2
    real(kind=8) :: c(dimdef), ck(dimdef), cs(dimdef)
    integer(kind=8) :: addeme, addep1, addep2, addete, adde2nd, ii, jj
    real(kind=8) :: defgep(dimdef), defgem(dimdef)
    real(kind=8) :: dfdi(nno, 3), dfdi2(nnos, 3), b(dimdef, dimuel)
    real(kind=8) :: drds(dimdef+1, dimcon), drdsr(dimdef, dimcon), dsde(dimcon, dimdef)
    real(kind=8) :: r(dimdef+1), sigbar(dimdef)
    real(kind=8) :: work1(dimcon, dimuel), work2(dimdef, dimuel)
    real(kind=8), allocatable, dimension(:, :) :: matri
!
! --------------------------------------------------------------------------------------------------
!
    allocate (matri(dimmat, dimmat))
    ASSERT(nddls*nnos .le. dimmat)
    ASSERT(dimuel .le. dimmat)
    codret = 0
    defgep = 0.d0
    defgem = 0.d0
    dfdi = 0.d0
    dfdi2 = 0.d0
    b = 0.d0
    drds = 0.d0
    drdsr = 0.d0
    dsde = 0.d0
    r = 0.d0
    sigbar = 0.d0
    work1 = 0.d0
    work2 = 0.d0
    addeme = mecani(2)
    addep1 = press1(3)
    addep2 = press2(3)
    addete = tempe(2)
    adde2nd = second(2)

! - Get parameters for behaviour
    call thmGetBehaviour(compor, ds_thm)

! - Get parameters for internal variables
    call thmGetBehaviourVari(ds_thm)

! - Some checks between behaviour and model
    call thmGetBehaviourChck(ds_thm)

! - Get storage parameters for behaviours
    nb_vari_meca = ds_thm%ds_behaviour%nb_vari_meca

! - Get initial parameters (THM_INIT)
    call thmGetParaInit(j_mater, ds_thm, l_check_=ASTER_TRUE)

! - Time parameters
    time_incr = time_curr-time_prev
    parm_theta = carcri(PARM_THETA_THM)

! - Create matrix for selection of dof for reduced integration
    call thmSelectMatrix(ds_thm, &
                         ndim, dimdef, inte_type, &
                         addeme, addete, addep1, addep2, adde2nd, &
                         a, as, &
                         c, cs)

! - Initialization of output fields
    if (lVect) then
        vectu(1:dimuel) = 0.d0
    end if
    if (lMatr) then
        matuu(1:dimuel*dimuel) = 0.d0
        matri(:, :) = 0.d0
    end if

! - Loop on integration points
    do ipi = 1, npi
        kpi = ipi

! ----- Compute [B] matrix for generalized strains
        call cabthm(ds_thm, l_axi, ndim, &
                    nddls, nddlm, &
                    nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                    nno, nnos, &
                    dimuel, dimdef, kpi, &
                    addeme, addete, addep1, addep2, adde2nd, &
                    elem_coor, &
                    jv_poids, jv_poids2, &
                    jv_func, jv_func2, &
                    jv_dfunc, jv_dfunc2, &
                    dfdi, dfdi2, &
                    poids, poids2, &
                    b)

! ----- Compute generalized strains
        do i = 1, dimdef
            defgem(i) = 0.d0
            defgep(i) = 0.d0
            do n = 1, dimuel
                defgem(i) = defgem(i)+b(i, n)*dispm(n)
                defgep(i) = defgep(i)+b(i, n)*dispp(n)
            end do
        end do

! ----- Compute generalized stresses and derivatives at current Gauss point
        call equthm(ds_thm, option, j_mater, &
                    lMatr, lSigm, &
                    lVari, lMatrPred, &
                    typmod, angl_naut, parm_theta, &
                    ndim, nbvari, &
                    kpi, npg, &
                    dimdef, dimcon, &
                    mecani, press1, press2, tempe, second, &
                    carcri, &
                    defgem, defgep, &
                    congem((kpi-1)*dimcon+1), congep((kpi-1)*dimcon+1), &
                    vintm((kpi-1)*nbvari+1), vintp((kpi-1)*nbvari+1), &
                    time_prev, time_curr, time_incr, &
                    r, drds, dsde, codret)

! ----- For selective integrations => move Gauss points to nodes
        if (ds_thm%ds_elem%l_dof_meca) then
            if (kpi .gt. npg) then
                if (lSigm) then
                    do i = 1, 6
                        congep((kpi-1)*dimcon+i) = congep((kpi-npg-1)*dimcon+i)
                    end do
                end if
                if (lVari) then
                    do i = 1, nb_vari_meca
                        vintp((kpi-1)*nbvari+i) = vintp((kpi-npg-1)*nbvari+i)
                    end do
                end if
            end if
        end if

        if (codret .ne. 0) then
            goto 99
        end if
! ======================================================================
! --- CONTRIBUTION DU POINT D'INTEGRATION KPI A LA MATRICE TANGENTE ET -
! --- AU RESIDU --------------------------------------------------------
! ----------------------------------------------------------------------
! --- MATRICE TANGENTE : REMPLISSAGE EN NON SYMETRIQUE -----------------
! ======================================================================
! --- CHOIX DU JEU DE MATRICES ADAPTE AU POINT D'INTEGRATION -----------
! --- SI KPI<NPG ALORS ON EST SUR UN POINT DE GAUSS: CK = C  -----------
! --- SINON ON EST SUR UN SOMMET                   : CK = CS -----------
! ======================================================================
        if (kpi .le. npg) then
            ck(1:dimdef) = c(1:dimdef)
            ak(1:2) = a(1:2)
        else
            ck(1:dimdef) = cs(1:dimdef)
            ak(1:2) = as(1:2)
        end if
! ======================================================================
! --- CALCUL DE MATUU (MATRI) ------------------------------------------
! --- ON MODIFIE LA 7EME LIGNE (TERME EN TEMPERATURE) DE LA MATRICE ----
! --- DE MANIERE A L'ADAPTER AU PI -------------------------------------
! ======================================================================
        if (lMatr) then
            do i = 1, dimdef
                do j = 1, dimcon
                    drdsr(i, j) = drds(i, j)
                end do
            end do
!
            if (ds_thm%ds_elem%l_dof_ther) then
                do i = 1, dimcon
                    drdsr(addete, i) = ak(1)*drds(addete, i)+ak(2)*drds(dimdef+1, i)
                end do
            end if
! ======================================================================
! --- ON ASSEMBLE: DF=BT.CK.DRDSR.DK.DSDE.FK.B.POIDS -------------------
! ======================================================================
            call pmathm(dimmat, dimdef, dimcon, dimuel, dsde, &
                        drdsr, ck, b, poids, work1, work2, matri)
!
        end if
! ======================================================================
! --- CALCUL DE VECTUU -------------------------------------------------
! ======================================================================
! --- ON SELECTIONNE LES COMPOSANTES UTILES DE R POUR CE PI ------------
! ======================================================================
        if (lVect) then
            do i = 1, dimdef
                sigbar(i) = ck(i)*r(i)
            end do
! ======================================================================
! --- ON SELECTIONNE LA BONNE COMPOSANTE 7 POUR CE PI ------------------
! ======================================================================
            if (ds_thm%ds_elem%l_dof_ther) then
                sigbar(addete) = ak(1)*r(addete)+ak(2)*r(dimdef+1)
            end if
! ======================================================================
! --- ON ASSEMBLE R=BT.SIGBAR.POIDS ------------------------------------
! ======================================================================
            do i = 1, dimuel
                do k = 1, dimdef
                    vectu(i) = vectu(i)+b(k, i)*sigbar(k)*poids
                end do
            end do
        end if
    end do
! ======================================================================
! --- SORTIE DE BOUCLE SUR LES POINTS D'INTEGRATION --------------------
! ======================================================================
    if (lMatr) then
        kji = 1
        do ii = 1, dimuel
            do jj = 1, dimuel
                matuu(kji) = matri(ii, jj)
                kji = kji+1
            end do
        end do
    end if
! ======================================================================
99  continue
    deallocate (matri)
! ======================================================================
end subroutine
