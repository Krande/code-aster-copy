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
! aslint: disable=W1504
!
subroutine nmrepl(modele, numedd, ds_material, carele, ds_system, &
                  ds_constitutive, lischa, ds_algopara, fonact, iterat, &
                  ds_measure, sdpilo, sdnume, ds_contact, &
                  deltat, valinc, solalg, veelem, veasse, &
                  sddisc, etan, ds_conv, eta, offset, &
                  ldccvg, pilcvg, matass)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/ismaem.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/infdbg.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmceta.h"
#include "asterfort/nmcha0.h"
#include "asterfort/nmchai.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmchso.h"
#include "asterfort/nmfext.h"
#include "asterfort/nmpilo.h"
#include "asterfort/nmpilr.h"
#include "asterfort/nmrelp.h"
#include "asterfort/nmrep2.h"
!
    integer(kind=8) :: fonact(*)
    integer(kind=8) :: iterat
    real(kind=8) :: deltat, eta, etan, offset
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    character(len=19) :: lischa, sdnume, sdpilo, sddisc, matass
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Contact), intent(in) :: ds_contact
    type(NL_DS_Measure), intent(inout) :: ds_measure
    type(NL_DS_System), intent(in) :: ds_system
    character(len=24) :: modele, numedd, carele
    character(len=19) :: veelem(*), veasse(*)
    character(len=19) :: solalg(*), valinc(*)
    type(NL_DS_Conv), intent(inout) :: ds_conv
    integer(kind=8) :: pilcvg, ldccvg
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE)
!
! CHOIX DU ETA DE PILOTAGE AVEC RECHERCHE LINEAIRE
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  NUMEDD : NUME_DDL
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  LISCHA : LISTE DES CHARGES
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IN  ITERAT : NUMERO D'ITERATION DE NEWTON
! IN  SDPILO : SD PILOTAGE
! IN  SDNUME : SD NUMEROTATION
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algopara      : datastructure for algorithm parameters
! In  ds_contact       : datastructure for contact management
! IO  ds_system        : datastructure for non-linear system management
! IN  DELTAT : INCREMENT DE TEMPS
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  ETAN   : ETA_PILOTAGE AU DEBUT DE L'ITERATION
! IN  SDDISC : SD DISCRETISATION
! IO  ds_conv          : datastructure for convergence management
! OUT ETA    : PARAMETRE DE PILOTAGE
! OUT RHO    : PARAMETRE DE RECHERCHE_LINEAIRE
! OUT OFFSET : DECALAGE DE ETA_PILOTAGE EN FONCTION DE RHO
! OUT PILCVG : CODE DE CONVERGENCE POUR LE PILOTAGE
!                -1 : PAS DE CALCUL DU PILOTAGE
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : PAS DE SOLUTION
!                 2 : BORNE ATTEINTE -> FIN DU CALCUL
! OUT LDCCVG : CODE RETOUR DE L'INTEGRATION DU COMPORTEMENT
!                -1 : PAS D'INTEGRATION DU COMPORTEMENT
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : ECHEC DE L'INTEGRATION DE LA LDC
!                 3 : SIZZ PAS NUL POUR C_PLAN DEBORST
! IN  MATASS : SD MATRICE ASSEMBLEE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8), parameter :: zsolal = 17, zvalin = 28
    character(len=19) :: solalt(zsolal), valint(zvalin, 2)
    integer(kind=8), parameter :: zveass = 19
    character(len=19) :: veasst(zveass)
    aster_logical :: exopt, mieux, irecli
    integer(kind=8) :: itrlmx, iterho, act, opt
    integer(kind=8) :: pilopt
    integer(kind=8) :: nbeffe
    integer(kind=8) :: nr, pos, nbsto, n, nbatte, nmax
    real(kind=8) :: rhomin, rhomax, rhoexm, rhoexp
    real(kind=8) :: f0, fopt, fcvg
    real(kind=8) :: rho, rhoopt, proeta(2)
    real(kind=8) :: r(1002), g(1002), memfg(1002), relirl
    real(kind=8) :: fgmax, fgmin, amelio, residu, etaopt
    character(len=19) :: cnfint2(2), cndiri2(2)
    character(len=19) :: cndiri, cnfext, cnsstr, k19bla
    character(len=19) :: depplu, sigplu, varplu, complu
    character(len=19) :: veinte2(2), depdet
    character(len=19) :: sigplt, varplt, depplt
    character(len=24) :: typilo
    character(len=24), pointer :: pltk(:) => null()
    type(NL_DS_System) :: ds_system2
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('PILOTAGE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<PILOTAGE> ... PILOTAGE AVEC RECH_LINE'
    end if
!
! --- INITIALISATIONS
!
    fopt = r8maem()
    pilopt = ismaem()
    nbsto = 0
    exopt = .false.
    irecli = .true.
    pilcvg = -1
    ldccvg = -1
    k19bla = ' '
    call nmchai('VEASSE', 'LONMAX', nmax)
    ASSERT(nmax .eq. zveass)
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
!
! --- DECOMPACTION VARIABLES CHAPEAUX
!
    call nmchex(valinc, 'VALINC', 'DEPPLU', depplu)
    call nmchex(valinc, 'VALINC', 'SIGPLU', sigplu)
    call nmchex(valinc, 'VALINC', 'VARPLU', varplu)
    call nmchex(valinc, 'VALINC', 'COMPLU', complu)
    call nmchex(veasse, 'VEASSE', 'CNDIRI', cndiri)
    call nmchex(veasse, 'VEASSE', 'CNSSTR', cnsstr)
!
! - Copy datastructure for solving system
!
    ds_system2 = ds_system
!
! --- LECTURE DONNEES PILOTAGE
!
    call jeveuo(sdpilo(1:19)//'.PLTK', 'L', vk24=pltk)
    typilo = pltk(1)
!
! --- FONCTIONS DE PILOTAGE LINEAIRES : RECHERCHE LINEAIRE STANDARD
!
    if (typilo .eq. 'DDL_IMPO') then
        call nmrelp(modele, numedd, ds_material, carele, ds_system, &
                    ds_constitutive, lischa, fonact, iterat, ds_measure, &
                    sdnume, ds_algopara, ds_contact, valinc, &
                    solalg, veelem, veasse, ds_conv, ldccvg)
        goto 999
    end if
!
! --- PREPARATION DES ZONES TEMPORAIRES POUR ITERATION COURANTE
!
    cnfint2(1) = ds_system%cnfint
    cnfint2(2) = '&&CNREPL.CHP1'
    cndiri2(1) = cndiri
    cndiri2(2) = '&&CNREPL.CHP2'
    veinte2(1) = ds_system%veinte
    veinte2(2) = '&&CNREPL.CHPX'
    depdet = '&&CNREPL.CHP3'
    depplt = '&&CNREPL.CHP4'
    sigplt = '&&NMREPL.SIGPLU'
    varplt = '&&NMREPL.VARPLU'
    call copisd('CHAMP_GD', 'V', varplu, varplt)
    call copisd('CHAMP_GD', 'V', sigplu, sigplt)
    call copisd('CHAMP_GD', 'V', depplu, depplt)
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
    call nmchso(veasse, 'VEASSE', 'CNDIRI', cndiri2(1), veasst)
!
! --- CALCUL DE F(RHO=0)
!
    call nmpilr(fonact, numedd, matass, veasse, ds_contact, ds_system%cnfint, &
                etan, f0)
    fcvg = abs(relirl*f0)
!
! --- INITIALISATION ET DIRECTION DE DESCENTE
!
    nr = 2
    r(1) = 0.d0
    r(2) = 1.d0
    g(1) = f0
    pos = 2
    nbatte = 2
!
! --- BOUCLE DE RECHERCHE LINEAIRE
!
    rho = 1.d0
    act = 1
!
    do iterho = 0, itrlmx
!
! ----- RESOLUTION DE L'EQUATION DE PILOTAGE: NVELLE DIRECT. DE DESCENTE
!
        call nmpilo(sdpilo, deltat, rho, solalg, veasse, &
                    modele, ds_material, ds_constitutive, valinc, &
                    nbatte, numedd, nbeffe, proeta, pilcvg, &
                    carele)
        if (pilcvg .eq. 1) goto 999
!
! ----- DECALAGE DU ETA_PILOTAGE
!
        offset = etan*(1-rho)
        do n = 1, nbeffe
            proeta(n) = proeta(n)+offset
        end do
!
! ----- Get right vectors
!
        call nmchso(veasse, 'VEASSE', 'CNDIRI', cndiri2(act), veasst)
        ds_system2%veinte = veinte2(act)
        ds_system2%cnfint = cnfint2(act)
!
! ----- Select ETA
!
        call nmceta(modele, numedd, ds_material, carele, &
                    ds_constitutive, ds_contact, lischa, fonact, ds_measure, &
                    sdpilo, iterat, sdnume, valint(1, act), solalg, &
                    veelem, veasst, sddisc, nbeffe, irecli, &
                    proeta, offset, rho, eta, ldccvg, &
                    pilcvg, residu, matass, ds_system2)
!
! ----- PB CVG: S'IL EXISTE DEJA UN RHO OPTIMAL, ON LE CONSERVE ET ON SORT
!
        if (ldccvg .gt. 0) then
            if (exopt) goto 100
            goto 999
        end if
!
! ---    SI ON A PAS ENCORE CONVERGE LE PILO :
! ---      * ON PREND UN PILO CONVERGE QQ SOIT LE RESIDU
! ---    SINON :
! ---      * ON CHERCHE A BAISSER LE RESIDU AVEC UN PILO CONVERGE
!
        if (pilopt .gt. 0) then
            mieux = ((pilcvg .eq. 0) .or. (pilcvg .eq. 2) .or. (residu .lt. fopt))
        else
            mieux = (((pilcvg .eq. 0) .or. (pilcvg .eq. 2)) .and. (residu .lt. fopt))
        end if
!
        if (mieux) then
            exopt = .true.
            rhoopt = rho
            etaopt = eta
            pilopt = pilcvg
            fopt = residu
            opt = act
            act = 3-act
        end if
!
! ---   MEMOIRE DES RESIDUS ATTEINTS
!
        nbsto = nbsto+1
        memfg(nbsto) = residu
!
! ---   ARRET SI SATISFACTION DU CRITERE
!
        if (residu .lt. fcvg) goto 100
!
! ---   ARRET SI IL N'Y A PLUS D'AMELIORATIONS SIGNIFICATIVES
!
        if (nbsto .ge. 3) then
            fgmax = max(memfg(nbsto), memfg(nbsto-1), memfg(nbsto-2))
            fgmin = min(memfg(nbsto), memfg(nbsto-1), memfg(nbsto-2))
            amelio = fgmin/fgmax
            if (amelio .gt. 0.95d0) goto 100
        end if
!
! ---   CALCUL DE RHO(N+1) PAR INTERPOLATION QUADRATIQUE AVEC BORNES
!
        g(pos) = residu
        call nmrep2(nr, r, g, fcvg, rhomin, &
                    rhomax, rhoexm, rhoexp, pos)
        rho = r(pos)
    end do
    iterho = itrlmx
!
! --- STOCKAGE DU RHO OPTIMAL ET DES CHAMPS CORRESPONDANTS
!
100 continue
!
! --- CALCUL DE ETA_PILOTAGE
!
    eta = etaopt
    rho = rhoopt
!
! --- REACTUALISATION DES EFFORTS EXTERIEURS (AVEC ETA)
!
    call nmchex(veasst, 'VEASSE', 'CNFEXT', cnfext)
    call nmfext(eta, fonact, veasst, cnfext, ds_contact)
!
! --- RECUPERATION DES VARIABLES EN T+ (PAS DE RECALCUL)
!
    if (opt .ne. 1) then
        call copisd('CHAMP_GD', 'V', sigplt, sigplu)
        call copisd('CHAMP_GD', 'V', varplt, varplu)
        call copisd('CHAMP_GD', 'V', cnfint2(opt), ds_system%cnfint)
        call copisd('CHAMP_GD', 'V', cndiri2(opt), cndiri)
    end if
!
! - Save results of line search
!
    ds_conv%line_sear_coef = rhoopt
    ds_conv%line_sear_iter = iterho
    pilcvg = pilopt
999 continue
!
! --- LE CALCUL DE PILOTAGE A FORCEMENT ETE REALISE
!
    ASSERT(pilcvg .ge. 0)
!
end subroutine
