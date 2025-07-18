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
subroutine digouj(option, rela_comp, nno, nbt, neq, &
                  nc, icodma, dul, sim, varim, &
                  pgl, klv, klc, varip, fono, &
                  sip, nomte)
! ----------------------------------------------------------------------
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mavec.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vecma.h"
    integer(kind=8) :: nbt, neq, icodma, nc
    real(kind=8) :: dul(neq), sim(neq), sip(neq), varim(*)
    real(kind=8) :: pgl(3, 3), klv(nbt), varip(*), fono(neq), klc(neq, neq)
    character(len=16) :: option, rela_comp, nomte
!
!  COMPORTEMENT DIS_GOUJON : APPLICATION : GOUJ2ECH
!           RELATION DE COMPORTEMENT : ELASTIQUE PARTOUT
!           SAUF SUIVANT Y LOCAL : DIS_GOUJON
!       ELEMENTS MECA_DIS_T_L
!
! IN  : NBT    : NOMBRE DE VALEURS POUR LA DEMI-MATRICE
!       NEQ    : NOMBRE DE DDL DE L'ELEMENT
!       NC     : NOMBRE DE DDL PAR NOEUD = 3 OU 6
!       ICODMA : ADRESSE DU MATERIAU CODE
!       DUL    : INCREMENT DE DEPLACEMENT REPERE LOCAL
!       SIM    : EFFORTS GENERALISES A L'INSTANT PRECEDENT
!       TP     : INSTANT ACTUEL
!       VARIM  : VARIABLE INTERNE A L'INSTANT PRECEDENT
!       PGL    : MATRICE DE PASSAGE REPERE GLOBAL -> LOCAL
!
! VAR : KLV    : MATRICE ELASTIQUE REPERE LOCAL EN ENTREE
!              : MATRICE TANGENTE EN SORTIE
! OUT : VARIP  : VARIABLE INTERNE REACTUALISEE
!       FONI   : FORCES NODALES
!       SIP    : EFFORTS INTERNES
!
    integer(kind=8) :: i, nno, jprolp, jvalep, nbvalp, lgpg, jtab(7)
    real(kind=8) :: seuil
    real(kind=8) :: dfl(6), fl(6)
    real(kind=8) :: nu, para_vale, valpap
    character(len=8) :: nompar, para_type
    character(len=24) :: valk(2)
    aster_logical :: plasti
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret
    real(kind=8) :: a, airerp, coef, deps, dp, dsidep, dut
    real(kind=8) :: e, rp, rprim, sieleq, sigdv, sigel, sigeps
    real(kind=8) :: sigy
!-----------------------------------------------------------------------
    valpap = 0.d0
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)
!
    if (nc .ne. 2) then
        valk(1) = nomte
        valk(2) = rela_comp
        call utmess('F', 'ELEMENTS_31', nk=2, valk=valk)
    end if
!
! --- CALCUL ELASTIQUE
!
! --- DEMI-MATRICE KLV TRANSFORMEE EN MATRICE PLEINE KLC
!
    call vecma(klv, nbt, klc, neq)
!
! --- CALCUL DE FL = KLC.DUL (INCREMENT D'EFFORT ELASTIQUE)
!
    call pmavec('ZERO', neq, klc, dul, dfl)
    dut = dul(2+nc)-dul(2)
!
    ASSERT(rela_comp(1:10) .eq. 'DIS_GOUJ2E')
!
    call rctype(icodma, 0, nompar, [valpap], para_vale, &
                para_type)
    call rctrac(icodma, 1, 'SIGM', para_vale, jprolp, &
                jvalep, nbvalp, e)
    if (rela_comp .eq. 'DIS_GOUJ2E_PLAS') then
        call rcfonc('S', 1, jprolp, jvalep, nbvalp, &
                    sigy=sigy)
        call rcfonc('V', 1, jprolp, jvalep, nbvalp, &
                    p=varim(1), rp=rp, rprim=rprim, airerp=airerp)
        plasti = (varim(2) .ge. 0.5d0)
    else if (rela_comp .eq. 'DIS_GOUJ2E_ELAS') then
        sigy = 0.d0
        rp = 0.d0
        plasti = .false.
    end if
!
    deps = dut
!
    sigel = sim(2)+e*deps
    sieleq = abs(sigel)
    seuil = sieleq-rp
!
    dp = 0.d0
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
!
        do i = 1, nc
            sip(i) = -dfl(i)+sim(i)
            sip(i+nc) = dfl(i+nc)+sim(i+nc)
            fl(i) = dfl(i)-sim(i)
            fl(i+nc) = dfl(i+nc)+sim(i+nc)
        end do
!
        if (rela_comp .eq. 'DIS_GOUJ2E_ELAS') then
            sip(2) = sigel
!
        else if (rela_comp .eq. 'DIS_GOUJ2E_PLAS') then
            if (seuil .le. 0.d0) then
                varip(2) = 0.d0
                dp = 0.d0
            else
                varip(2) = 1.d0
                nu = 0.5d0
                call rcfonc('E', 1, jprolp, jvalep, nbvalp, &
                            e=e, nu=nu, p=varim(1), rp=rp, rprim=rprim, &
                            airerp=airerp, sieleq=sieleq, dp=dp)
            end if
            varip(1) = varim(1)+dp
            plasti = (varip(2) .ge. 0.5d0)
!
            sip(2) = sigel*rp/(rp+e*dp)
            varip(1+lgpg) = varip(1)
            varip(2+lgpg) = varip(2)
        end if
        sip(2+nc) = sip(2)
!
!        FL : EFFORTS GENERALISES AUX NOEUDS 1 ET 2 (REPERE LOCAL)
!            ON CHANGE LE SIGNE DES EFFORTS SUR LE PREMIER NOEUD
!        FONO : FORCES NODALES AUX NOEUDS 1 ET 2 (REPERE GLOBAL)
!
        fl(2) = -sip(2)
        fl(2+nc) = sip(2)
        if (nomte .eq. 'MECA_2D_DIS_T_L') then
            call ut2vlg(nno, nc, pgl, fl, fono)
        else
            call utpvlg(nno, nc, pgl, fl, fono)
        end if
    end if
!
    if (option(1:14) .eq. 'RIGI_MECA_TANG' .or. option(1:9) .eq. 'FULL_MECA') then
!
        if (option(1:14) .eq. 'RIGI_MECA_TANG') then
!         - - OPTION='RIGI_MECA_TANG' => SIGMA(T)
            rp = 0.d0
            sigdv = sim(2)
            rp = abs(sigdv)
        else
!         - - OPTION='FULL_MECA' => SIGMA(T+DT)
            sigdv = sip(2)
        end if
!
        a = 1.d0
        dsidep = 0.d0
        if (rela_comp .eq. 'DIS_GOUJ2E_PLAS') then
            sigeps = 0.d0
            sigeps = sigeps+sigdv*deps
            if (plasti .and. (sigeps .ge. 0.d0)) then
                a = 1.d0+e*dp/rp
                coef = -e**2/(e+rprim)/rp**2*(1.d0-dp*rprim/rp)/a
                dsidep = coef*sigdv*sigdv
            end if
        end if
        dsidep = dsidep+e/a
    end if
!
    if (option .eq. 'FULL_MECA' .or. option .eq. 'RIGI_MECA_TANG') then
        if (nc .eq. 2) then
!            KLV(3)  =  DSIDEP
!            KLV(10)  = DSIDEP
            klc(2, 2) = dsidep
            klc(4, 4) = dsidep
            klc(2, 4) = -dsidep
            klc(4, 2) = -dsidep
        else if (nc .eq. 3) then
!            KLV(3)  =  DSIDEP
!            KLV(15)  = DSIDEP
            klc(2, 2) = dsidep
            klc(5, 5) = dsidep
            klc(2, 5) = -dsidep
            klc(5, 2) = -dsidep
        end if
        call mavec(klc, neq, klv, nbt)
    end if
end subroutine
