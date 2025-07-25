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
! aslint: disable=W1504
!
subroutine nmnewt(mesh, model, numins, numedd, numfix, &
                  ds_material, cara_elem, ds_constitutive, list_load, ds_system, &
                  hhoField, &
                  sddyna, nlDynaDamping, &
                  ds_algopara, fonact, ds_measure, sderro, ds_print, &
                  sdnume, sddisc, sdcrit, sdsuiv, &
                  sdpilo, ds_conv, solveu, maprec, matass, &
                  ds_inout, valinc, solalg, meelem, measse, &
                  veelem, veasse, ds_contact, ds_algorom, eta, &
                  nbiter)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use Rom_Datastructure_type
    use HHO_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/isfonc.h"
#include "asterfort/nmactf.h"
#include "asterfort/nmactn.h"
#include "asterfort/nmaffi.h"
#include "asterfort/nmconv.h"
#include "asterfort/nmcrel.h"
#include "asterfort/nmcvgf.h"
#include "asterfort/nmcvgn.h"
#include "asterfort/nmdepl.h"
#include "asterfort/nmdesc.h"
#include "asterfort/nmeceb.h"
#include "asterfort/nmeraz.h"
#include "asterfort/nmevdt.h"
#include "asterfort/nmevr0.h"
#include "asterfort/nmfcon.h"
#include "asterfort/nmfcor.h"
#include "asterfort/nmible.h"
#include "asterfort/nmimcr.h"
#include "asterfort/nmimr0.h"
#include "asterfort/nmleeb.h"
#include "asterfort/nmnble.h"
#include "asterfort/nmnpas.h"
#include "asterfort/nmpred.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmrini.h"
#include "asterfort/nmsuiv.h"
#include "asterfort/nmtble.h"
#include "asterfort/nmtime.h"
#include "asterfort/nmtimr.h"
#include "asterfort/nmforc_step.h"
#include "asterfort/infdbg.h"
#include "asterfort/nonlinDSPrintSepLine.h"
!
    character(len=8), intent(in) :: mesh
    character(len=24), intent(in) :: model
    integer(kind=8) :: numins
    character(len=24) :: numedd
    character(len=24) :: numfix
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=24), intent(in) :: cara_elem
    type(NL_DS_Constitutive), intent(inout) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
    type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
    integer(kind=8) :: fonact(*)
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(HHO_Field), intent(in) :: hhoField
    character(len=24) :: sderro
    type(NL_DS_Print), intent(inout) :: ds_print
    character(len=19) :: sdnume
    character(len=19), intent(in) :: list_load, sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    character(len=19) :: sddisc
    character(len=19) :: sdcrit
    character(len=24) :: sdsuiv
    character(len=19) :: sdpilo
    type(NL_DS_Conv), intent(inout) :: ds_conv
    character(len=19) :: solveu
    character(len=19) :: maprec
    character(len=19) :: matass
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19) :: valinc(*)
    character(len=19) :: solalg(*)
    character(len=19) :: meelem(*)
    character(len=19) :: measse(*)
    character(len=19) :: veelem(*)
    character(len=19) :: veasse(*)
    type(NL_DS_Contact), intent(inout) :: ds_contact
    type(ROM_DS_AlgoPara), intent(inout) :: ds_algorom
    real(kind=8) :: eta
    integer(kind=8) :: nbiter
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Newton
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  ds_material      : datastructure for material parameters
! In  cara_elem        : name of elementary characteristics (field)
! In  list_load        : name of datastructure for list of loads
! IO  ds_algopara      : datastructure for algorithm parameters
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IO  ds_constitutive  : datastructure for constitutive laws management
! In  solver           : name of datastructure for solver
! In  ds_system        : datastructure for non-linear system management
! In  sd_suiv          : datastructure for dof monitoring parameters
! In  sd_obsv          : datastructure for observation parameters
! IO  ds_inout         : datastructure for input/output management
! IO  ds_energy        : datastructure for energy management
! IO  ds_conv          : datastructure for convergence management
! IO  ds_contact       : datastructure for contact management
! IO  ds_measure       : datastructure for measure and statistics management
! IO  ds_algorom       : datastructure for ROM parameters
! IO  ds_print         : datastructure for printing parameters
! IO  ds_algorom       : datastructure for ROM parameters
! I/O ETA    : PARAMETRE DE PILOTAGE
! OUT NBITER : NOMBRE D'ITERATIONS DE NEWTON
! OUT ETATIN : ETAT DE LA CONVERGENCE DU PAS DE TEMPS
!     0 - CVG  -> LE PAS DE TEMPS A CONVERGE
!     1 - NOOK -> UN EVENEMENT DURANT LE PAS DE TEMPS
!     2 - NCVG -> LE PAS DE TEMPS N'A PAS CONVERGE
!     3 - STOP -> ERREUR FATALE - ARRET DU CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: niveau, iterat
    aster_logical :: lerrit
    aster_logical :: l_loop_exte, l_cont_disc, l_cont, l_hrom_corref
    character(len=4) :: etnewt, etfixe
    real(kind=8) :: time
!
! --------------------------------------------------------------------------------------------------
!
    iterat = 0
    niveau = 0
    nbiter = 0
    lerrit = ASTER_FALSE
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call nonlinDSPrintSepLine()
    end if
!
! - Active functionnalities
!
    l_loop_exte = isfonc(fonact, 'BOUCLE_EXTERNE')
    l_cont_disc = isfonc(fonact, 'CONT_DISCRET')
    l_cont = isfonc(fonact, 'CONTACT')
    l_hrom_corref = isfonc(fonact, 'HROM_CORR_EF')
!
! - Reset events
!
    call nmeraz(sderro, 'TOUS')
    call nmevr0(sddisc)
!
! - Reset values in convergence table for Newton loop
!
    call nmimr0(ds_print, 'NEWT')
!
! - Supplementary loops
!
    if (l_loop_exte) then
        if (l_cont) then
            niveau = 3
        elseif (l_hrom_corref) then
            niveau = 10
        else
            ASSERT(.false.)
        end if
    end if
!
! - Updates for new time step
!
    call nmnpas(mesh, model, cara_elem, &
                fonact, list_load, &
                ds_material, ds_constitutive, &
                ds_measure, ds_print, &
                sddisc, numins, &
                sdsuiv, sddyna, &
                ds_contact, ds_conv, &
                sdnume, numedd, solveu, &
                valinc, solalg)
!
! - Compute forces for second member when constant in time step
!
    call nmforc_step(fonact, &
                     model, cara_elem, numedd, &
                     list_load, sddyna, &
                     ds_material, ds_constitutive, &
                     ds_measure, ds_inout, &
                     sddisc, numins, &
                     valinc, solalg, hhoField, &
                     veelem, veasse)
!
! ======================================================================
!     BOUCLE POINTS FIXES
! ======================================================================
!
100 continue
!
    iterat = 0
    nbiter = nbiter+1
    ds_contact%iteration_newton = iterat
    ds_contact%it_adapt_maxi = ds_conv%iter_glob_maxi
!
! - External loop management - BEGIN
!
    call nmible(niveau, model, ds_contact, &
                fonact, ds_measure, ds_print, ds_algorom)
!
! - External loop management - Initializations for new loop
!
    call nmnble(mesh, model, fonact, sddisc, numins, &
                sddyna, sdnume, numedd, ds_measure, ds_contact, &
                valinc, solalg)
!
! ======================================================================
!     PREDICTION
! ======================================================================
!
!
! - Launch timer
!
    call nmtime(ds_measure, 'Launch', 'Newt_Iter')
!
! --- PREDICTION D'UNE DIRECTION DE DESCENTE
!
    call nmpred(model, numedd, numfix, ds_material, cara_elem, &
                ds_constitutive, list_load, ds_algopara, solveu, ds_system, &
                fonact, ds_print, ds_measure, ds_algorom, sddisc, &
                sdnume, sderro, numins, valinc, solalg, &
                matass, maprec, ds_contact, &
                sddyna, nlDynaDamping, &
                meelem, measse, veelem, veasse, lerrit)
!
    if (lerrit) goto 315
!
! ======================================================================
!     BOUCLE SUR LES ITERATIONS DE NEWTON
! ======================================================================
!
300 continue
!
! - Launch timer
!
    if (iterat .ne. 0) then
        call nmtime(ds_measure, 'Launch', 'Newt_Iter')
    end if
!
! --- CALCUL PROPREMENT DIT DE L'INCREMENT DE DEPLACEMENT
! --- EN CORRIGEANT LA (LES) DIRECTIONS DE DESCENTE
! --- SI CONTACT OU PILOTAGE OU RECHERCHE LINEAIRE
!
    call nmdepl(model, numedd, ds_material, cara_elem, &
                ds_constitutive, list_load, fonact, ds_measure, ds_algopara, &
                mesh, numins, iterat, solveu, matass, &
                sddyna, nlDynaDamping, &
                sddisc, sdnume, sdpilo, sderro, &
                ds_contact, valinc, solalg, veelem, veasse, &
                eta, ds_conv, ds_system, lerrit)
!
    if (lerrit) goto 315
!
! --- CALCUL DES FORCES APRES CORRECTION
!
    call nmfcor(model, numedd, ds_material, cara_elem, ds_system, &
                ds_constitutive, list_load, fonact, ds_algopara, numins, &
                iterat, ds_measure, sddisc, &
                sddyna, nlDynaDamping, &
                sdnume, sderro, ds_contact, &
                valinc, solalg, &
                veelem, veasse, measse, matass, lerrit)
!
    if (lerrit) goto 315
!
! - DOF monitoring
!
    call nmsuiv(mesh, sdsuiv, ds_print, cara_elem, model, &
                ds_material, ds_constitutive, valinc, sddisc, numins)
!
! --- ESTIMATION DE LA CONVERGENCE
!
315 continue
    call nmconv(mesh, model, ds_material, numedd, sdnume, fonact, &
                sddyna, nlDynaDamping, &
                ds_conv, ds_print, ds_measure, &
                sddisc, sdcrit, sderro, ds_algopara, ds_algorom, &
                ds_inout, matass, solveu, ds_system, numins, &
                iterat, eta, ds_contact, valinc, &
                solalg, measse, veasse)
!
! --- MISE A JOUR DES EFFORTS DE CONTACT
!
    call nmfcon(model, numedd, ds_material, fonact, ds_contact, &
                ds_measure, valinc, solalg, ds_constitutive)
!
! - Evaluate events at current Newton iteration
!
    call nmcvgn(sddisc, sderro, valinc, ds_contact)
!
! - Print during Newton loop
!
    call nmaffi(fonact, ds_conv, ds_print, sderro, sddisc, &
                'NEWT')
!
! - Stop Newton iterations
!
    call nmleeb(sderro, 'NEWT', etnewt)

    if (etnewt .ne. 'CONT') then
        goto 330
    end if
!
! --- ON CONTINUE LES ITERATIONS DE NEWTON : CALCUL DE LA DESCENTE
!
320 continue
!
    call nmdesc(model, numedd, &
                numfix, ds_material, cara_elem, &
                ds_constitutive, list_load, ds_contact, &
                ds_algopara, ds_system, solveu, &
                fonact, numins, iterat, &
                sddisc, ds_print, ds_measure, &
                ds_algorom, sddyna, nlDynaDamping, sdnume, &
                sderro, matass, maprec, &
                valinc, solalg, meelem, &
                measse, veasse, lerrit)
!
    if (lerrit) goto 315
!
! --- ON CONTINUE NEWTON
!
    iterat = iterat+1
    nbiter = nbiter+1
    ds_contact%iteration_newton = iterat
!
! --- CAS DU CONTACT DISCRET
!
    call nmleeb(sderro, 'NEWT', etnewt)
    if (etnewt .eq. 'CTCD') then
        call nmeceb(sderro, 'NEWT', 'CONT')
        call nmtime(ds_measure, 'Stop', 'Newt_Iter')
        goto 300
    end if
!
330 continue
!
! - Timer for current Newton iteration (not for prediction)
!
    call nmtime(ds_measure, 'Stop', 'Newt_Iter')
    call nmrinc(ds_measure, 'Newt_Iter')
    call nmtimr(ds_measure, 'Newt_Iter', 'N', time)
    call nmimcr(ds_print, 'ITER_TIME', time, .true._1)
!
! --- VERIFICATION DU DECLENCHEMENT DES ERREURS FATALES
!
    call nmevdt(ds_measure, sderro, 'ITE')
!
! - Reset times and counters
!
    call nmrini(ds_measure, 'N')
!
! --- ON CONTINUE NEWTON ?
!
    call nmleeb(sderro, 'NEWT', etnewt)
    if (etnewt .eq. 'CONT') goto 300
!
! ======================================================================
!     FIN BOUCLE SUR LES ITERATIONS DE NEWTON
! ======================================================================
!
!
!
! --- GESTION DES ACTIONS A LA FIN DE LA BOUCLE DE NEWTON
!
    call nmactn(ds_print, sddisc, sderro, ds_contact, &
                ds_conv, iterat, numins)
!
! --- ON FAIT DES ITERATIONS SUPPLEMENTAIRES ?
!
    call nmleeb(sderro, 'NEWT', etnewt)
    if (etnewt .eq. 'CONT') then
        call nmtime(ds_measure, 'Launch', 'Newt_Iter')
        call nmcrel(sderro, 'ITER_MAXI', .false._1)
        goto 320
    end if
!
! - External loop management - END
!
    call nmtble(niveau, model, mesh, ds_contact, &
                fonact, ds_print, &
                sderro, ds_conv, sddisc, numins, valinc, &
                solalg, ds_algorom)
!
! --- ETAT DE LA CONVERGENCE POINT FIXE
!
    call nmcvgf(sddisc, sderro, valinc, ds_contact)
!
! --- GESTION DES ACTIONS A LA FIN D'UNE BOUCLE DE POINT FIXE
!
    call nmactf(ds_print, sddisc, sderro, ds_contact, &
                ds_conv, iterat, numins)
!
! --- POUR LA CONTINUATION DU POINT FIXE: GLUTE DUE AU CONTACT DISCRET
!
    call nmleeb(sderro, 'FIXE', etfixe)
    if (etfixe .eq. 'CONT') then
        if (l_cont_disc) then
            if (.not. ds_conv%l_stop) then
                call nmeceb(sderro, 'FIXE', 'CONV')
            else
                call nmeceb(sderro, 'NEWT', 'CTCD')
                call nmtime(ds_measure, 'Launch', 'Newt_Iter')
                goto 320
            end if
        else if (l_loop_exte) then
            goto 100
        else
            call nmeceb(sderro, 'FIXE', 'CONV')
        end if
    end if
!
! ======================================================================
!     FIN BOUCLE POINTS FIXES
! ======================================================================
!
    if (niv .ge. 2) then
        call nonlinDSPrintSepLine()
    end if
!
end subroutine
