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
subroutine ndnpas(fonact, numedd, numins, sddisc, sddyna, &
                  valinc, solalg)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ndmuap.h"
#include "asterfort/ndpred.h"
#include "asterfort/ndynin.h"
#include "asterfort/ndynkk.h"
#include "asterfort/ndynlo.h"
#include "asterfort/ndynre.h"
#include "asterfort/nmdebg.h"
#include "blas/dcopy.h"
!
    integer :: numins
    character(len=24) :: numedd
    character(len=19) :: sddyna, sddisc
    character(len=19) :: solalg(*), valinc(*)
    integer :: fonact(*)
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (DYNAMIQUE)
!
! INITIALISATION DES CHAMPS D'INCONNUES POUR UN NOUVEAU PAS DE TEMPS
!
! ----------------------------------------------------------------------
!
!
! IN  FONACT : FONCTIONNALITES ACTIVEES
! IN  NUMEDD : NUME_DDL
! IN  NUMINS : NUMERO INSTANT COURANT
! IN  SDDISC : SD DISCRETISATION TEMPORELLE
! IN  RESOCO : SD RESOLUTION DU CONTACT
! IN  SDDYNA : SD DYNAMIQUE
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  SOLALG : VARIABLE CHAPEAU POUR INCREMENTS SOLUTIONS
!
! ----------------------------------------------------------------------
!
    real(kind=8) :: zero, un, deux
    parameter(un=1.d0, deux=2.d0)
    parameter(zero=0.d0)
!
    character(len=24) :: cfsc
    real(kind=8), pointer :: coef_sch(:) => null()
    real(kind=8) :: alpha, beta, gamma, phi
    real(kind=8) :: instam, instap, deltat
    aster_logical :: lexge, lctcc, lmuap, lgrot, lexpl, lmpas, lhhtc, limpl
    aster_logical :: ldepl, lvite, lacce
    aster_logical :: lnewma
    real(kind=8) :: coerig, coeamo, coemas
    real(kind=8) :: coeext, coeint, coeequ, coeex2
    integer :: imode
    integer :: neq, nbmodp
    real(kind=8) :: coefd(3), coefv(3), coefa(3)
    real(kind=8) :: coedep, coevit, coeacc
    real(kind=8) :: coerma, coeram, coerri
    real(kind=8) :: coiner
    character(len=19) :: depgem, vitgem, accgem, depgep, vitgep, accgep
    integer :: jdepgm, jvitgm, jaccgm, jdepgp, jvitgp, jaccgp
    integer :: ifm, niv
!
! ----------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> INITIALISATIONS EN DYNAMIQUE'
    end if
!
! --- INITIALISATIONS
!
    call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
    instam = diinst(sddisc, numins-1)
    instap = diinst(sddisc, numins)
    deltat = instap-instam
!
! --- FONCTIONNALITES ACTIVEES
!
    lexge = ndynlo(sddyna, 'EXPL_GENE')
    lctcc = isfonc(fonact, 'CONT_CONTINU')
    lmuap = ndynlo(sddyna, 'MULTI_APPUI')
    lgrot = isfonc(fonact, 'GD_ROTA')
    lexpl = ndynlo(sddyna, 'EXPLICITE')
    lmpas = ndynlo(sddyna, 'MULTI_PAS')
    limpl = ndynlo(sddyna, 'IMPLICITE')
!
! --- ACCES SD DYNA
!
    cfsc = sddyna(1:15)//'.COEF_SCH'
    call jeveuo(cfsc, 'E', vr=coef_sch)
!
! --- TYPE DE FORMULATION SCHEMA DYNAMIQUE GENERAL
!
    ldepl = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 1
    lvite = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 2
    lacce = ndynin(sddyna, 'FORMUL_DYNAMIQUE') .eq. 3
    if (lgrot .and. .not. ldepl) then
        ASSERT(.false.)
    end if
!
! --- TYPE DE SCHEMA: NEWMARK (ET SES DERIVEES)
!
    lnewma = ndynlo(sddyna, 'FAMILLE_NEWMARK')
    if (.not. (lnewma)) then
        ASSERT(.false.)
    end if
!
! --- HHT COMPLET (MULTI-PAS)
!
    lhhtc = ndynlo(sddyna, 'HHT_COMPLET')
!
! --- COEFFICIENTS DU SCHEMA EN TEMPS
!
    beta = ndynre(sddyna, 'BETA')
    gamma = ndynre(sddyna, 'GAMMA')
    phi = ndynre(sddyna, 'PHI')
    alpha = ndynre(sddyna, 'ALPHA')
!
! --- COEFFICIENTS POUR MATRICES
!
    if (lnewma) then
        if (ldepl) then
            coerig = un
            coeamo = gamma/(beta*deltat)
            coemas = un/(beta*deltat*deltat)
        else if (lacce) then
            coerig = beta*deltat*deltat
            coeamo = gamma*deltat
            coemas = un
        else
            ASSERT(.false.)
        end if
        if (lhhtc) then
            coeamo = coeamo/(un+alpha)
            coemas = coemas/(un+alpha)
        end if
    else
        ASSERT(.false.)
    end if
!
    coef_sch(1) = coerig
    coef_sch(2) = coeamo
    coef_sch(3) = coemas
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... COEF. RIGI.: ', coerig
        write (ifm, *) '<MECANONLINE> ... COEF. AMOR.: ', coeamo
        write (ifm, *) '<MECANONLINE> ... COEF. MASS.: ', coemas
    end if
!
! --- COEFFICIENTS POUR MISE A JOUR DEPL/VITE/ACCE
!
    coedep = 1.d0
    if (lnewma) then
        if (ldepl) then
            coedep = un
            coevit = gamma/(beta*deltat)
            coeacc = un/(beta*deltat*deltat)
        else if (lacce) then
            coedep = beta*deltat*deltat
            coevit = gamma*deltat
            coeacc = un
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
    end if
    coef_sch(13) = coedep
    coef_sch(14) = coevit
    coef_sch(15) = coeacc
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... COEF. DEPL.: ', coedep
        write (ifm, *) '<MECANONLINE> ... COEF. VITE.: ', coevit
        write (ifm, *) '<MECANONLINE> ... COEF. ACCE.: ', coeacc
    end if
!
! --- COEFFICIENTS POUR PREDICTEURS
!
    if (lnewma) then
        if (ldepl) then
            coefd(1) = zero
            coefd(2) = zero
            coefd(3) = zero
            coefv(1) = zero
            coefv(2) = (beta-gamma)/beta
            coefv(3) = ((deux*beta-gamma)*deltat)/(deux*beta)
            coefa(1) = zero
            coefa(2) = -un/(beta*deltat)
            coefa(3) = (deux*beta-un)/(deux*beta)
        else if (lacce) then
            if (lexpl) then
                if (ndynlo(sddyna, 'TCHAMWA')) then
                    coefd(1) = un
                    coefd(2) = deltat
                    coefd(3) = deltat*deltat*phi
                    coefv(1) = zero
                    coefv(2) = un
                    coefv(3) = deltat
                    coefa(1) = zero
                    coefa(2) = zero
                    coefa(3) = zero
                else
                    coefd(1) = un
                    coefd(2) = deltat
                    coefd(3) = deltat*deltat/deux
                    coefv(1) = zero
                    coefv(2) = un
                    coefv(3) = deltat*(un-gamma)
                    coefa(1) = zero
                    coefa(2) = zero
                    coefa(3) = zero
                end if
            else
                coefd(1) = un
                coefd(2) = deltat
                coefd(3) = deltat*deltat/deux
                coefv(1) = zero
                coefv(2) = un
                coefv(3) = deltat
                coefa(1) = zero
                coefa(2) = zero
                coefa(3) = un
            end if
        else
            ASSERT(.false.)
        end if
    else
        ASSERT(.false.)
    end if
!
    coef_sch(4) = coefd(1)
    coef_sch(5) = coefd(2)
    coef_sch(6) = coefd(3)
    coef_sch(7) = coefv(1)
    coef_sch(8) = coefv(2)
    coef_sch(9) = coefv(3)
    coef_sch(10) = coefa(1)
    coef_sch(11) = coefa(2)
    coef_sch(12) = coefa(3)
!
! --- CALCUL DES PREDICTEURS
!
    call ndpred(sddyna, valinc, solalg)
!
! --- COEFFICIENTS POUR SCHEMAS A PLUSIEURS PAS
! --- COEEXT: COEF. DE PONDERATION DES FORCES EXTERNES
! --- COEINT: COEF. DE PONDERATION DES FORCES INTERNES
! --- COEEQU: COEF. PERMETTANT DE RESPECTER L'EQUILIBRE SU RLES AUTRES
!             TERMES NON PONDERES
!
    if (lmpas) then
        if (lhhtc) then
            coeext = -alpha/(un+alpha)
            coeint = -alpha/(un+alpha)
            coeequ = un/(un+alpha)
            coeex2 = un
        else
            ASSERT(.false.)
        end if
    else
        coeext = zero
        coeint = zero
        coeequ = un
        coeex2 = un
    end if
    coef_sch(16) = coeext
    coef_sch(17) = coeequ
    coef_sch(18) = coeint
    coef_sch(19) = coeex2
!
    if (lmpas) then
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE> ... MULTI-PAS F. EXT. N-1: ', coeext
            write (ifm, *) '<MECANONLINE> ... MULTI-PAS F. EXT. N  : ', coeex2
            write (ifm, *) '<MECANONLINE> ... MULTI-PAS F. INT. N-1: ', coeint
            write (ifm, *) '<MECANONLINE> ... MULTI-PAS F. EQU.    : ', coeequ
        end if
    end if
!
! --- COEFFICENT POUR CALCUL FORCE D'INERTIE DE REFERENCE (NDINER)
!
    if (lnewma) then
        if (limpl) then
            coiner = un/(beta*deltat)
        else
            coiner = un/deltat
        end if
    else
        coiner = un/deltat
    end if
    coef_sch(23) = coiner
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... COEF. FORC. INERTIE REF: ', &
            coiner
    end if
!
! --- COEFFICIENTS DEVANT MATRICE POUR TERME DE RAPPEL DYNAMIQUE
!
    coerma = un
    coeram = un
    coerri = un
    coef_sch(20) = coerma
    coef_sch(21) = coeram
    coef_sch(22) = coerri
!
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... COEF. FDYNA RIGI: ', coerri
        write (ifm, *) '<MECANONLINE> ... COEF. FDYNA AMOR: ', coeram
        write (ifm, *) '<MECANONLINE> ... COEF. FDYNA MASS: ', coerma
    end if
!
! - Save previous time
!
    coef_sch(24) = instam
!
! --- INITIALISATION DES CHAMPS D'ENTRAINEMENT EN MULTI-APPUI
!
    if (lmuap) then
        call ndmuap(numins, numedd, sddyna, sddisc)
    end if
!
! --- INITIALISATION DES DEPL. GENERALISES SI PROJECTION MODALE
!
    if (lexge) then
        call ndynkk(sddyna, 'PRMO_DEPGEM', depgem)
        call ndynkk(sddyna, 'PRMO_VITGEM', vitgem)
        call ndynkk(sddyna, 'PRMO_ACCGEM', accgem)
        call ndynkk(sddyna, 'PRMO_DEPGEP', depgep)
        call ndynkk(sddyna, 'PRMO_VITGEP', vitgep)
        call ndynkk(sddyna, 'PRMO_ACCGEP', accgep)
        nbmodp = ndynin(sddyna, 'NBRE_MODE_PROJ')
        call jeveuo(accgem, 'E', jaccgm)
        call jeveuo(accgep, 'E', jaccgp)
        call jeveuo(vitgem, 'E', jvitgm)
        call jeveuo(vitgep, 'E', jvitgp)
        call jeveuo(depgem, 'E', jdepgm)
        call jeveuo(depgep, 'E', jdepgp)
        call dcopy(nbmodp, zr(jdepgm), 1, zr(jdepgp), 1)
        call dcopy(nbmodp, zr(jvitgm), 1, zr(jvitgp), 1)
        call dcopy(nbmodp, zr(jaccgm), 1, zr(jaccgp), 1)
!
! --- PREDICTION DEPLACEMENT GENERALISE
!
        do imode = 1, nbmodp
            zr(jdepgp+imode-1) = zr(jdepgm+imode-1)+coefd(2)*zr(jvitgm+imode-1)+coefd(3)*zr(&
                                 &jaccgm+imode-1)
        end do
        if (niv .ge. 2) then
            write (ifm, *) '<MECANONLINE> ...... PRED. DEPL. GENE'
            call nmdebg('VECT', depgep, ifm)
        end if
!
    end if
!
    call jedema()
!
end subroutine
