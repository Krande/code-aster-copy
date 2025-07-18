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
subroutine nmrelp(model, nume_dof, ds_material, cara_elem, ds_system, &
                  ds_constitutive, list_load, list_func_acti, iter_newt, ds_measure, &
                  sdnume, ds_algopara, ds_contact, valinc, solalg, &
                  veelem, veasse, ds_conv, ldccvg, sddyna_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmcha0.h"
#include "asterfort/nmchai.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmchso.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nonlinRForceCompute.h"
#include "asterfort/nonlinIntForce.h"
#include "asterfort/nmmaji.h"
#include "asterfort/nmrebo.h"
#include "asterfort/nmrech.h"
#include "asterfort/nmrecz.h"
#include "asterfort/vlaxpy.h"
#include "asterfort/vtcreb.h"
#include "asterfort/vtzero.h"
#include "asterfort/zbinit.h"
#include "blas/daxpy.h"
!
    integer(kind=8) :: list_func_acti(*)
    integer(kind=8) :: iter_newt, ldccvg
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19) :: list_load, sdnume
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24) :: model, nume_dof, cara_elem
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19) :: veelem(*), veasse(*)
    character(len=19) :: solalg(*), valinc(*)
    type(NL_DS_Conv), intent(inout) :: ds_conv
    character(len=19), intent(in), optional :: sddyna_
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! RECHERCHE LINEAIRE DANS LA DIRECTION DE DESCENTE
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  cara_elem        : name of elementary characteristics (field)
! In  ds_material      : datastructure for material parameters
! In  list_load        : name of datastructure for list of loads
! In  nume_dof         : name of numbering object (NUME_DDL)
! In  ds_constitutive  : datastructure for constitutive laws management
! IO  ds_measure       : datastructure for measure and statistics management
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IN  ITERAT : NUMERO D'ITERATION DE NEWTON
! IN  SDNUME : SD NUMEROTATION
! In  ds_contact       : datastructure for contact management
! In  ds_system        : datastructure for non-linear system management
! In  ds_algopara      : datastructure for algorithm parameters
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! OUT LDCCVG : CODE RETOUR DE L'INTEGRATION DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
! IO  ds_conv          : datastructure for convergence management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: phaseType = CORR_NEWTON
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: zsolal = 17, zvalin = 28
    character(len=19) :: solalt(zsolal), valint(zvalin, 2)
    integer(kind=8) :: itrlmx, iterho, neq, act, opt, ldcopt
    integer(kind=8) :: dimmem, nmax
    real(kind=8) :: rhomin, rhomax, rhoexm, rhoexp
    real(kind=8) :: rhom, rhoopt, rho
    real(kind=8) :: f0, fm, f, fopt, fcvg
    real(kind=8) :: parmul, relirl, sens
    real(kind=8) :: mem(2, 10)
    aster_logical :: stite, lnkry
    aster_logical :: lgrot, lendo
    character(len=19) :: cnfint2(2), cndiri2(2)
    character(len=19) :: cndiri, cnfext, cnsstr, k19bla
    character(len=19) :: depplu, sigplu, varplu, complu
    character(len=19) :: sigplt, varplt, depplt
    character(len=19) :: vediri
    character(len=19) :: depdet, ddepla, depdel, sddyna
    aster_logical :: echec
    real(kind=8), pointer :: vale(:) => null()
    type(NL_DS_System) :: ds_system2
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... RECHERCHE LINEAIRE'
    end if
!
! --- FONCTIONNALITES ACTIVEES
!
    sddyna = ' '
    if (present(sddyna_)) then
        sddyna = sddyna_
    end if
    lgrot = isfonc(list_func_acti, 'GD_ROTA')
    lendo = isfonc(list_func_acti, 'ENDO_NO')
    lnkry = isfonc(list_func_acti, 'NEWTON_KRYLOV')
!
! --- INITIALISATIONS
!
    opt = 1
    parmul = 3.d0
    fopt = r8maem()
    k19bla = ' '
    ldccvg = -1
    call nmchai('VALINC', 'LONMAX', nmax)
    ASSERT(nmax .eq. zvalin)
    call nmchai('SOLALG', 'LONMAX', nmax)
    ASSERT(nmax .eq. zsolal)
!
! --- PARAMETRES RECHERCHE LINEAIRE
!
    itrlmx = ds_algopara%line_search%iter_maxi
    rhomin = ds_algopara%line_search%rho_mini
    rhomax = ds_algopara%line_search%rho_maxi
    rhoexm = -ds_algopara%line_search%rho_excl
    rhoexp = ds_algopara%line_search%rho_excl
    relirl = ds_algopara%line_search%resi_rela
    ASSERT(itrlmx .le. 1000)
    dimmem = 10
!
! --- DECOMPACTION VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPPLU', depplu)
    call nmchex(valinc, 'VALINC', 'SIGPLU', sigplu)
    call nmchex(valinc, 'VALINC', 'VARPLU', varplu)
    call nmchex(valinc, 'VALINC', 'COMPLU', complu)
    call nmchex(veasse, 'VEASSE', 'CNDIRI', cndiri)
    call nmchex(veasse, 'VEASSE', 'CNSSTR', cnsstr)
    call nmchex(veasse, 'VEASSE', 'CNFEXT', cnfext)
    call nmchex(veelem, 'VEELEM', 'CNDIRI', vediri)
    call nmchex(solalg, 'SOLALG', 'DDEPLA', ddepla)
    call nmchex(solalg, 'SOLALG', 'DEPDEL', depdel)
!
! - Copy datastructure for solving system
!
    ds_system2 = ds_system
!
! --- ACCES VARIABLES
!
    call jeveuo(ddepla(1:19)//'.VALE', 'E', vr=vale)
!
! --- PREPARATION DES ZONES TEMPORAIRES POUR ITERATION COURANTE
!
    cnfint2(1) = ds_system%cnfint
    cnfint2(2) = '&&NMRECH.RESI'
    cndiri2(1) = cndiri
    cndiri2(2) = '&&NMRECH.DIRI'
    depdet = '&&CNPART.CHP1'
    depplt = '&&CNPART.CHP2'
    sigplt = '&&NMRECH.SIGP'
    varplt = '&&NMRECH.VARP'
    call vtzero(depdet)
    call vtzero(depplt)
    call copisd('CHAMP_GD', 'V', varplu, varplt)
    call copisd('CHAMP_GD', 'V', sigplu, sigplt)
    call vtcreb('&&NMRECH.RESI', 'V', 'R', nume_ddlz=nume_dof, nb_equa_outz=neq)
    call vtcreb('&&NMRECH.DIRI', 'V', 'R', nume_ddlz=nume_dof, nb_equa_outz=neq)
!
! --- CONSTRUCTION DES VARIABLES CHAPEAUX
!
    call nmcha0('VALINC', 'ALLINI', ' ', valint(1, 1))
    call nmchso(valinc, 'VALINC', '      ', k19bla, valint(1, 1))
    call nmchso(valint(1, 1), 'VALINC', 'DEPPLU', depplt, valint(1, 1))
    call nmcha0('VALINC', 'ALLINI', ' ', valint(1, 2))
    call nmchso(valinc, 'VALINC', '      ', k19bla, valint(1, 2))
    call nmchso(valint(1, 2), 'VALINC', 'DEPPLU', depplt, valint(1, 2))
    call nmchso(valint(1, 2), 'VALINC', 'SIGPLU', sigplt, valint(1, 2))
    call nmchso(valint(1, 2), 'VALINC', 'VARPLU', varplt, valint(1, 2))
    call nmchso(solalg, 'SOLALG', 'DEPDEL', depdet, solalt)
!
! --- CALCUL DE F(RHO=0)
!
    call nmrecz(nume_dof, ds_contact, list_func_acti, cndiri, ds_system%cnfint, &
                cnfext, cnsstr, ddepla, f0)
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... FONCTIONNELLE INITIALE: ', f0
    end if
!
! --- VALEUR DE CONVERGENCE
!
    fcvg = abs(relirl*f0)
!
! --- INITIALISATION ET DIRECTION DE DESCENTE
!
    if (ds_algopara%line_search%method .eq. 'CORDE') then
        sens = 1.d0
        rhom = 0.d0
        fm = f0
        rhoopt = 1.d0
    else if (ds_algopara%line_search%method .eq. 'MIXTE') then
        if (f0 .le. 0.d0) then
            sens = 1.d0
        else
            sens = -1.d0
        end if
        call zbinit(sens*f0, parmul, dimmem, mem)
        rhoopt = 1.d0
    else
        ASSERT(.false.)
    end if
!
! --- BOUCLE DE RECHERCHE LINEAIRE
!
    rho = sens
    act = 1
!
    do iterho = 0, itrlmx
!
! ----- CALCUL DE L'INCREMENT DE DEPLACEMENT TEMPORAIRE
!
        call nmmaji(nume_dof, lgrot, lendo, sdnume, rho, &
                    depdel, ddepla, depdet, 0)
        call nmmaji(nume_dof, lgrot, lendo, sdnume, rho, &
                    depplu, ddepla, depplt, 1)
        if (lnkry) then
            call vlaxpy(1.d0-rho, ddepla, depdet)
            call vlaxpy(1.d0-rho, ddepla, depplt)
        end if
! ----- Print
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE> ...... ITERATION <', iterho, '>'
            write (ifm, *) '<MECANONLINE> ...... RHO COURANT = ', rho
            write (ifm, *) '<MECANONLINE> ...... INCREMENT DEPL.'
            call nmdebg('VECT', depplt, 6)
            write (ifm, *) '<MECANONLINE> ...... INCREMENT DEPL. TOTAL'
            call nmdebg('VECT', depdet, 6)
        end if
! ----- Update internal forces
        ds_system2%cnfint = cnfint2(act)
        ds_system2%veinte = ds_system%veinte
        call nonlinIntForce(phaseType, model, cara_elem, list_func_acti, iter_newt, &
                            sdnume, ds_material, ds_constitutive, ds_system2, ds_measure, &
                            valint(1, act), solalt, ldccvg)
! ----- Update force for Dirichlet boundary conditions (dualized) - BT.LAMBDA
        call nonlinRForceCompute(model, ds_material, cara_elem, list_load, nume_dof, &
                                 ds_measure, depplt, veelem, cndiri_=cndiri2(act))
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE> ...... FORCES INTERNES'
            call nmdebg('VECT', cnfint2(act), 6)
            write (ifm, *) '<MECANONLINE> ...... REACTIONS D''APPUI'
            call nmdebg('VECT', cndiri2(act), 6)
        end if
!
! ----- ON A NECESSAIREMENT INTEGRE LA LOI DE COMPORTEMENT
!
        ASSERT(ldccvg .ne. -1)
!
! ----- ECHEC A L'INTEGRATION DE LA LOI DE COMPORTEMENT
!
        if (ldccvg .ne. 0) then
!
! ------- S'IL EXISTE DEJA UN RHO OPTIMAL, ON LE CONSERVE
!
            if (iterho .gt. 0) then
                goto 100
            else
                goto 999
            end if
        end if
!
! ----- CALCUL DE F(RHO)
!
        call nmrecz(nume_dof, ds_contact, list_func_acti, cndiri2(act), cnfint2(act), &
                    cnfext, cnsstr, ddepla, f)
!
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE> ... FONCTIONNELLE COURANTE: ', f
        end if
!
! ----- CALCUL DU RHO OPTIMAL
!
        if (ds_algopara%line_search%method .eq. 'CORDE') then
            call nmrech(fm, f, fopt, fcvg, rhomin, &
                        rhomax, rhoexm, rhoexp, rhom, rho, &
                        rhoopt, ldcopt, ldccvg, opt, act, &
                        stite)
!
        else if (ds_algopara%line_search%method .eq. 'MIXTE') then
            call nmrebo(f, mem, sens, rho, rhoopt, &
                        ldcopt, ldccvg, fopt, fcvg, opt, &
                        act, rhomin, rhomax, rhoexm, rhoexp, &
                        stite, echec)
            if (echec) then
                goto 100
            end if
        else
            ASSERT(.false.)
        end if
        if (stite) then
            goto 100
        end if
    end do
    iterho = itrlmx
!
! --- STOCKAGE DU RHO OPTIMAL ET DES CHAMPS CORRESPONDANTS
!
100 continue
!
! --- AJUSTEMENT DE LA DIRECTION DE DESCENTE
!
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call daxpy(b_n, rhoopt-1.d0, vale, b_incx, vale, &
               b_incy)
!
! --- RECUPERATION DES VARIABLES EN T+ SI NECESSAIRE
!
    if (opt .ne. 1) then
        call copisd('CHAMP_GD', 'V', sigplt, sigplu)
        call copisd('CHAMP_GD', 'V', varplt, varplu)
        call copisd('CHAMP_GD', 'V', cnfint2(opt), ds_system%cnfint)
        call copisd('CHAMP_GD', 'V', cndiri2(opt), cndiri)
    end if
!
! --- INFORMATIONS SUR LA RECHERCHE LINEAIRE
!
    ldccvg = ldcopt
!
999 continue
!
! - Save results of line search
!
    ds_conv%line_sear_coef = rhoopt
    ds_conv%line_sear_iter = iterho
!
    call detrsd('CHAMP', '&&NMRECH.RESI')
    call detrsd('CHAMP', '&&NMRECH.DIRI')
!
end subroutine
