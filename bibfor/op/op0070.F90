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
subroutine op0070()
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use Rom_Datastructure_type
    use NonLinearDyna_module
    use HHO_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/infmaj.h"
#include "asterfort/inidbg.h"
#include "asterfort/jerecu.h"
#include "asterfort/ndexpl.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmactp.h"
#include "asterfort/nmaffi.h"
#include "asterfort/nmarch.h"
#include "asterfort/nmcvgc.h"
#include "asterfort/nmcvgp.h"
#include "asterfort/nmdata.h"
#include "asterfort/nmeceb.h"
#include "asterfort/nmerro.h"
#include "asterfort/nmevdt.h"
#include "asterfort/nmfpas.h"
#include "asterfort/nmini0.h"
#include "asterfort/nminit.h"
#include "asterfort/nmleeb.h"
#include "asterfort/nmlost.h"
#include "asterfort/nmmeng.h"
#include "asterfort/nmnewt.h"
#include "asterfort/nmpost.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmstat.h"
#include "asterfort/nmtime.h"
#include "asterfort/onerrf.h"
#include "asterfort/titre.h"
#include "asterfort/setTimeListProgressBar.h"
!
! --------------------------------------------------------------------------------------------------
!
! STAT_NON_LINE
! DYNA_NON_LINE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: fonact(100)
    integer(kind=8), parameter :: zmeelm = 8
    integer(kind=8), parameter :: zmeass = 4
    integer(kind=8), parameter :: zveelm = 12
    integer(kind=8), parameter :: zveass = 19
    integer(kind=8), parameter :: zsolal = 17
    integer(kind=8), parameter :: zvalin = 28
!
! --- GESTION BOUCLES
!
    integer(kind=8) :: numins, nbiter
    character(len=4) :: etfixe, etinst, etcalc
!
! --- GESTION ERREUR
!
    integer(kind=8) :: lenout
    character(len=16) :: compex
!
    integer(kind=8) :: ibid
!
    real(kind=8) :: eta
!
    character(len=8) :: mesh
!
    character(len=16) :: k16bid
    character(len=19) :: list_load
    character(len=19) :: maprec
    character(len=19) :: matass, solver
    character(len=24) :: model, mateco, cara_elem, mater
    character(len=24) :: numedd, numfix
!
! --- FONCTIONNALITES ACTIVEES
!
    aster_logical :: lexpl, limpl, lstat
!
! --- STRUCTURES DE DONNEES
!
    character(len=24), parameter :: sderro = '&&OP0070.ERRE.'
    character(len=24) :: sd_suiv
    character(len=19), parameter :: sdpilo = '&&OP0070.PILO.', sdnume = '&&OP0070.NUME.ROTAT'
    character(len=19), parameter :: sddisc = '&&OP0070.DISC.', sdcrit = '&&OP0070.CRIT.'
    character(len=19) :: sddyna
    character(len=19) :: sd_obsv
    type(NL_DS_Print)        :: ds_print
    type(NL_DS_Conv)         :: ds_conv
    type(NL_DS_AlgoPara)     :: ds_algopara
    type(NL_DS_InOut)        :: ds_inout
    type(NL_DS_Contact)      :: ds_contact
    type(NL_DS_Measure)      :: ds_measure
    type(NL_DS_Energy)       :: ds_energy
    type(ROM_DS_AlgoPara)    :: ds_algorom
    type(NL_DS_Constitutive) :: ds_constitutive
    type(NL_DS_PostTimeStep) :: ds_posttimestep
    type(NL_DS_Material)     :: ds_material
    type(NL_DS_ErrorIndic)   :: ds_errorindic
    type(NL_DS_System)       :: ds_system
    type(HHO_Field)          :: hhoField
    type(NLDYNA_DAMPING)     :: nlDynaDamping
!
! --- VARIABLES CHAPEAUX
!
    character(len=19) :: valinc(zvalin), solalg(zsolal)
!
! --- MATR_ELEM, VECT_ELEM ET MATR_ASSE
!
    character(len=19) :: meelem(zmeelm), veelem(zveelm)
    character(len=19) :: measse(zmeass), veasse(zveass)
!
! ----------------------------------------------------------------------
!
    call titre()
    call infmaj()
    call inidbg()
!
    fonact(:) = 0
    solver = '&&OP0070.SOLVEUR'
    list_load = '&&OP0070.LISCHA'
    maprec = '&&OP0070.MAPREC'
!
! ======================================================================
!     RECUPERATION DES OPERANDES ET INITIALISATION
! ======================================================================
!
! --- ON STOCKE LE COMPORTEMENT EN CAS D'ERREUR AVANT MNL : COMPEX
! --- PUIS ON PASSE DANS LE MODE "VALIDATION DU CONCEPT EN CAS D'ERREUR"
!
    call onerrf(' ', compex, lenout)
    call onerrf('EXCEPTION+VALID', k16bid, ibid)
!
! - Creation of datastructures
!
    call nmini0(eta, numins, matass, &
                zmeelm, zmeass, zveelm, &
                zveass, zsolal, zvalin, &
                ds_print, ds_conv, ds_algopara, &
                ds_inout, ds_contact, ds_measure, &
                ds_energy, ds_material, sderro)

! - Read parameters
    call nmdata(model, mesh, mater, mateco, cara_elem, ds_constitutive, &
                list_load, solver, ds_conv, sddyna, ds_posttimestep, &
                ds_energy, ds_errorindic, ds_print, ds_algopara, &
                ds_inout, ds_contact, ds_measure, ds_algorom, &
                nlDynaDamping)

! - Initializations of datastructures
    call nminit(mesh, model, mater, mateco, cara_elem, list_load, &
                numedd, numfix, ds_algopara, ds_constitutive, maprec, &
                solver, numins, sddisc, sdnume, sdcrit, &
                ds_material, fonact, sdpilo, ds_print, &
                sddyna, nlDynaDamping, &
                sd_suiv, sd_obsv, sderro, ds_posttimestep, ds_inout, &
                ds_energy, ds_conv, ds_errorindic, valinc, solalg, &
                measse, veelem, meelem, veasse, ds_contact, &
                ds_measure, ds_algorom, ds_system, hhoField)
!
! - Launch timer for total time
!
    call nmtime(ds_measure, 'Launch', 'Compute')
!
! --- PREMIER INSTANT
!
    numins = 1
!
! --- QUELQUES FONCTIONNALITES ACTIVEES
!
    limpl = ndynlo(sddyna, 'IMPLICITE')
    lexpl = ndynlo(sddyna, 'EXPLICITE')
    lstat = ndynlo(sddyna, 'STATIQUE')
!
! ======================================================================
!  DEBUT DU PAS DE TEMPS
! ======================================================================
!
200 continue
!
! --- AUCUNE BOUCLE N'EST CONVERGE
!
    call nmeceb(sderro, 'RESI', 'CONT')
    call nmeceb(sderro, 'NEWT', 'CONT')
    call nmeceb(sderro, 'FIXE', 'CONT')
    call nmeceb(sderro, 'INST', 'CONT')
!
    call jerecu('V')
!
! - Launch timer for current step time
!
    call nmtime(ds_measure, 'Launch', 'Time_Step')
!
! --- REALISATION DU PAS DE TEMPS
!
    if (lexpl) then
        call ndexpl(model, numedd, ds_material, cara_elem, &
                    ds_constitutive, list_load, ds_algopara, fonact, ds_system, &
                    ds_print, ds_measure, sdnume, &
                    sddyna, nlDynaDamping, &
                    sddisc, sderro, valinc, numins, solalg, solver, &
                    matass, maprec, ds_inout, meelem, measse, &
                    veelem, veasse, nbiter)
    else if (lstat .or. limpl) then
        call nmnewt(mesh, model, numins, numedd, numfix, &
                    ds_material, cara_elem, ds_constitutive, list_load, ds_system, &
                    hhoField, &
                    sddyna, nlDynaDamping, &
                    ds_algopara, fonact, ds_measure, sderro, ds_print, &
                    sdnume, sddisc, sdcrit, sd_suiv, &
                    sdpilo, ds_conv, solver, maprec, matass, &
                    ds_inout, valinc, solalg, meelem, measse, &
                    veelem, veasse, ds_contact, ds_algorom, eta, &
                    nbiter)
    else
        ASSERT(ASTER_FALSE)
    end if
!
! - End of timer for current step time
!
    call nmtime(ds_measure, 'Stop', 'Time_Step')
    call nmrinc(ds_measure, 'Time_Step')
!
! ======================================================================
!  FIN DU PAS DE TEMPS
! ======================================================================
!
!
! - Time lost (step was cut)
!
    call nmleeb(sderro, 'FIXE', etfixe)
    if (etfixe .eq. 'ERRE') then
        call nmlost(ds_measure)
    end if
!
! - Post-treatment
!
    call nmpost(model, mesh, cara_elem, list_load, &
                numedd, numfix, ds_system, &
                ds_constitutive, ds_material, &
                ds_contact, ds_algopara, fonact, &
                ds_measure, sddisc, numins, eta, &
                sd_obsv, sderro, &
                sddyna, nlDynaDamping, &
                valinc, solalg, &
                meelem, measse, veasse, &
                ds_energy, ds_errorindic, &
                ds_posttimestep)
!
! --- ETAT DE LA CONVERGENCE DU PAS DE TEMPS
!
    call nmcvgp(sddisc, numins, sderro, valinc, fonact, &
                ds_contact)
!
! --- AFFICHAGES PENDANT LA BOUCLE DES PAS DE TEMPS
!
    call nmaffi(fonact, ds_conv, ds_print, sderro, sddisc, &
                'INST')
!
! --- STATISTIQUES SUR PAS DE TEMPS
!
    if (.not. lexpl) then
        call nmstat('P', ds_measure, ds_print, sddisc, numins, sderro)
    end if
!
! --- GESTION DES ACTIONS A LA FIN D'UN PAS DE TEMPS
!
    call nmactp(ds_print, sddisc, sderro, ds_contact, &
                ds_conv, nbiter, numins)
!
! --- INSTANT SUIVANT
!
    call nmleeb(sderro, 'INST', etinst)
    if (etinst .eq. 'ERRE') then
        goto 200
    else if (etinst .eq. 'STOP') then
        goto 800
    end if
!
! --- VERIFICATION DU DECLENCHEMENT DES ERREURS FATALES
!
    call nmevdt(ds_measure, sderro, 'PAS')
!
! --- EVALUATION DE LA CONVERGENCE DU CALCUL
!
    call nmcvgc(sddisc, sderro, numins, fonact)
!
! --- ARCHIVAGE DES RESULTATS
!
    call onerrf(compex, k16bid, ibid)
    call nmarch(numins, model, ds_material, cara_elem, fonact, &
                ds_print, sddisc, sdcrit, &
                ds_measure, sderro, sddyna, sdpilo, ds_energy, &
                ds_inout, ds_errorindic, ds_algorom)
    call onerrf('EXCEPTION+VALID', k16bid, ibid)
!
! --- ETAT DU CALCUL
!
    call nmleeb(sderro, 'CALC', etcalc)
    if ((etcalc .eq. 'ERRE') .or. (etcalc .eq. 'STOP')) then
        goto 800
    else if (etcalc .eq. 'CONV') then
        goto 900
    end if
!
! --- MISE A JOUR DES INFORMATIONS POUR UN NOUVEAU PAS DE TEMPS
!
    ASSERT(etcalc .eq. 'CONT')
    call nmfpas(fonact, sddyna, sdpilo, sddisc, nbiter, &
                numins, eta, valinc, solalg, veasse, ds_system, &
                ds_contact)
    numins = numins+1
!
    goto 200
!
! ======================================================================
!     GESTION DES ERREURS
! ======================================================================
!
800 continue
!
! --- ON COMMENCE PAR ARCHIVER LE PAS DE TEMPS PRECEDENT
!
    if (numins .ne. 1) then
        call nmarch(numins-1, model, ds_material, cara_elem, fonact, &
                    ds_print, sddisc, sdcrit, &
                    ds_measure, sderro, sddyna, sdpilo, ds_energy, &
                    ds_inout, ds_errorindic, ds_algorom)
    end if
!
! - Write messages for errors
!
    call nmerro(sderro, ds_measure, numins)
!
! ======================================================================
!     SORTIE
! ======================================================================
!
900 continue
!
! - Progress bar
!
    call setTimeListProgressBar(sddisc, numins, final_=ASTER_TRUE)
!
! - End of timer for total time
!
    call nmtime(ds_measure, 'Stop', 'Compute')
!
! --- IMPRESSION STATISTIQUES FINALES
!
    if (.not. lexpl) then
        call nmstat('T', ds_measure, ds_print, sddisc, numins, sderro)
    end if
!
! --- ON REMET LE MECANISME D'EXCEPTION A SA VALEUR INITIALE
!
    call onerrf(compex, k16bid, ibid)
!
! - Cleaning datastructures
!
    call nmmeng(fonact, &
                ds_algorom, ds_print, ds_measure, ds_material, &
                ds_energy, ds_inout, ds_posttimestep, hhoField)
!
end subroutine
