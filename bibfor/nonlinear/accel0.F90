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
subroutine accel0(model, numeDof, listFuncActi, listLoad, &
                  ds_contact, maprec, solveu, &
                  sddyna, nlDynaDamping, &
                  ds_measure, ds_system, &
                  hval_meelem, hval_measse, &
                  hval_veelem, hval_veasse, &
                  hval_incr, hval_algo)
!
    use NonLin_Datastructure_type
    use NonLinearDyna_type
!
    implicit none
!
#include "asterfort/copisd.h"
#include "asterfort/detlsp.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/lspini.h"
#include "asterfort/nmassi.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmprac.h"
#include "asterfort/nmreso.h"
#include "asterfort/utmess.h"
#include "asterfort/vtzero.h"
#include "asterfort/nonlinLoadDirichletCompute.h"
!
    character(len=19), intent(in) :: solveu, maprec, listLoad
    character(len=19), intent(in) :: sddyna
    type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=24), intent(in) :: numeDof, model
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: hval_meelem(*), hval_measse(*)
    character(len=19), intent(in) :: hval_veasse(*), hval_veelem(*)
    character(len=19), intent(in) :: hval_algo(*), hval_incr(*)
    integer(kind=8), intent(in) :: listFuncActi(*)
    type(NL_DS_System), intent(in) :: ds_system
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (DYNAMIQUE)
!
! CALCUL DE L'ACCELERATION INITIALE
!
! --------------------------------------------------------------------------------------------------
!
!     ==> ON SUPPOSE QUE LA VITESSE INITIALE EST NULLE
!                    QUE LES DEPLACEMENTS IMPOSES SONT NULS
!     ==> ON NE PREND EN COMPTE QUE LES CHARGES DYNAMIQUES, CAR LES
!         CHARGES STATIQUES SONT EQUILIBREES PAR LES FORCES INTERNES
!
!
! IN  MODELE : NOM DU MODELE
! IN  NUMEDD : NUME_DDL (VARIABLE AU COURS DU CALCUL)
! IN  LISCHA : LISTE DES CHARGES
! In  ds_contact       : datastructure for contact management
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_system        : datastructure for non-linear system management
! In  sddyna           : name of datastructure for dynamic parameters
! In  nlDynaDamping    : damping parameters
! IN  VEELEM : VARIABLE CHAPEAU POUR NOM DES VECT_ELEM
! IN  VEASSE : VARIABLE CHAPEAU POUR NOM DES VECT_ASSE
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: neq
    integer(kind=8) :: faccvg, rescvg
    character(len=19) :: depso1, depso2
    character(len=19) :: cncine, cncinx, cndonn, k19bla
    character(len=19) :: disp_prev, acce_prev
    character(len=19) :: matrAsse
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_26')
    end if
    call utmess('I', 'MECANONLINE_24')

! - Initializations
    k19bla = ' '
    cndonn = '&&CNCHAR.DONN'
    cncinx = '&&CNCHAR.CINE'
    matrAsse = '&&ACCEL0.MATASS'
    call dismoi('NB_EQUA', numeDof, 'NUME_DDL', repi=neq)

! - DECOMPACTION VARIABLES CHAPEAUX
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', disp_prev)
    call nmchex(hval_incr, 'VALINC', 'ACCMOI', acce_prev)
    call nmchex(hval_veasse, 'VEASSE', 'CNCINE', cncine)
    call nmchex(hval_algo, 'SOLALG', 'DEPSO1', depso1)
    call nmchex(hval_algo, 'SOLALG', 'DEPSO2', depso2)

! - ASSEMBLAGE ET FACTORISATION DE LA MATRICE
    call nmprac(listFuncActi, listLoad, numeDof, solveu, &
                sddyna, nlDynaDamping, &
                ds_measure, ds_contact, &
                hval_meelem, hval_measse, &
                maprec, matrAsse, &
                faccvg)
    if (faccvg .eq. 1) then
        call utmess('A', 'MECANONLINE_78')
    end if
    if (faccvg .eq. 2) then
        call vtzero(acce_prev)
        call utmess('A', 'MECANONLINE_69')
        goto 999
    end if

! - Compute values of Dirichlet conditions
    call nonlinLoadDirichletCompute(listLoad, model, numeDof, &
                                    ds_measure, matrAsse, disp_prev, &
                                    hval_veelem, hval_veasse)

! - Evaluate second member for initial acceleration
    call nmassi(listFuncActi, sddyna, nlDynaDamping, ds_system, hval_incr, hval_veasse, cndonn)

! --- POUR LE CALCUL DE DDEPLA, IL FAUT METTRE CNCINE A ZERO

    call copisd('CHAMP_GD', 'V', cncine, cncinx)
    call vtzero(cncinx)
!
! --- RESOLUTION DIRECTE
!
    call nmreso(listFuncActi, cndonn, k19bla, cncinx, solveu, &
                maprec, matrAsse, depso1, depso2, rescvg)
    if (rescvg .eq. 1) then
        call vtzero(acce_prev)
        call utmess('A', 'MECANONLINE_70')
        goto 999
    end if
!
! --- DEPENDAMMENT DU SOLVEUR, TRAITEMENT PARTICULIER
!
    call lspini(solveu)
!
! --- RECOPIE SOLUTION
!
    call copisd('CHAMP_GD', 'V', depso1, acce_prev)
!
999 continue
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ...... ACCMOI : '
        call nmdebg('VECT', acce_prev, ifm)
    end if

!
! --- MENAGE
!
! --- EN PREMIER DE L'EVENTUELLE INSTANCE MUMPS (SI PRE_COND 'LDLT_SP')
    call detlsp(matrAsse, solveu)
! --- EN SECOND DE LA MATRICE ASSEMBLEE
    call detrsd('MATR_ASSE', matrAsse)
!
end subroutine
