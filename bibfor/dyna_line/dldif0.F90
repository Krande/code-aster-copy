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
subroutine dldif0(result, force1, neq, istoc, iarchi, &
                  lamort, imat, masse, rigid, amort, &
                  dep0, vit0, acc0, depl1, vite1, &
                  acce1, vite2, fexte, famor, fliai, &
                  nchar, nveca, liad, lifo, modele, &
                  ener, mate, mateco, carele, charge, &
                  infoch, fomult, numedd, dt, temps, &
                  tabwk0, tabwk1, archiv, nbtyar, typear, &
                  numrep, ds_energy)
!
!     CALCUL MECANIQUE TRANSITOIRE PAR INTEGRATION DIRECTE
!     AVEC  METHODE EXPLICITE :  DIFFERENCES CENTREES
!     ------------------------------------------------------------------
!  IN  : NEQ       : NOMBRE D'EQUATIONS
!  IN  : ISTOC     : PILOTAGE DU STOCKAGE DES RESULTATS
!  IN  : IARCHI    : PILOTAGE DE L'ARCHIVAGE DES RESULTATS
!  IN  : LAMORT    : LOGIQUE INDIQUANT SI IL Y A AMORTISSEMENT
!  IN  : IMAT      : TABLEAU D'ADRESSES POUR LES MATRICES
!  IN  : MASSE     : MATRICE DE MASSE
!  IN  : NCHAR     : NOMBRE D'OCCURENCES DU MOT CLE CHARGE
!  IN  : NVECA     : NOMBRE D'OCCURENCES DU MOT CLE VECT_ASSE
!  IN  : LIAD      : LISTE DES ADRESSES DES VECTEURS CHARGEMENT (NVECT)
!  IN  : LIFO      : LISTE DES NOMS DES FONCTIONS EVOLUTION (NVECT)
!  IN  : MODELE    : NOM DU MODELE
!  IN  : MATE      : NOM DU CHAMP DE MATERIAU
!  IN  : CARELE    : CARACTERISTIQUES DES POUTRES ET COQUES
!  IN  : CHARGE    : LISTE DES CHARGES
!  IN  : INFOCH    : INFO SUR LES CHARGES
!  IN  : FOMULT    : LISTE DES FONC_MULT ASSOCIES A DES CHARGES
!  IN  : NUMEDD    : NUME_DDL DE LA MATR_ASSE RIGID
!  VAR : DEP0      : TABLEAU DES DEPLACEMENTS A L'INSTANT N
!  VAR : VIT0      : TABLEAU DES VITESSES A L'INSTANT N
!  VAR : ACC0      : TABLEAU DES ACCELERATIONS A L'INSTANT N
!  IN  : TEMPS     : INSTANT COURANT
! IN  NUMREP : NUMERO DE REUSE POUR LA TABLE PARA_CALC
!
! CORPS DU PROGRAMME
! aslint: disable=W1504
!
    use NonLin_Datastructure_type
!
    implicit none
!
! DECLARATION PARAMETRES D'APPELS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dlarch.h"
#include "asterfort/dlfdyn.h"
#include "asterfort/dlfext.h"
#include "asterfort/enerca.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmarpc.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
    integer(kind=8) :: neq, istoc, iarchi, ivit0r
    integer(kind=8) :: ifnobi, ifcibi
    integer(kind=8) :: archiv, nbtyar
    integer(kind=8) :: imat(3)
    integer(kind=8) :: nchar, nveca, liad(*)
!
    real(kind=8) :: dep0(*), vit0(*), acc0(*)
    real(kind=8) :: depl1(neq), vite1(neq), acce1(neq)
    real(kind=8) :: vite2(neq)
    real(kind=8) :: fexte(*), famor(*), fliai(*)
    real(kind=8) :: tabwk0(neq), tabwk1(neq)
    real(kind=8) :: dt, temps
!
    character(len=8) :: masse, rigid, amort
    character(len=16) :: typear(nbtyar)
    character(len=24) :: modele, mate, mateco, carele, charge, infoch, fomult, numedd
    character(len=24) :: lifo(*)
    type(NL_DS_Energy), intent(inout) :: ds_energy
    character(len=8) :: result
    character(len=19) :: force1
    integer(kind=8) :: numrep
!
    aster_logical :: lamort, ener
!
!
!
!
!
    integer(kind=8) :: iforc1, ieq, alarm
    real(kind=8) :: r8bid
    character(len=19) :: masse1, amort1, rigid1, k19bid
    blas_int :: b_incx, b_incy, b_n
!
! --- CALCUL DES DEPLACEMENTS ET VITESSES
!
    do ieq = 1, neq
        vite1(ieq) = vit0(ieq)+dt*acc0(ieq)
        depl1(ieq) = dep0(ieq)+dt*vite1(ieq)
    end do
!
!====
! 3. CALCUL DU SECOND MEMBRE F*
!====
!
    call jeveuo(force1(1:19)//'.VALE', 'E', iforc1)
!
    call dlfext(nveca, nchar, temps, neq, liad, &
                lifo, charge, infoch, fomult, modele, &
                mate, mateco, carele, numedd, zr(iforc1))
!
    if (ener) then
        do ieq = 1, neq
            fexte(ieq) = fexte(ieq+neq)
            fexte(ieq+neq) = zr(iforc1+ieq-1)
        end do
    end if
!
    call dlfdyn(imat(1), imat(3), lamort, neq, depl1, &
                vite1, zr(iforc1), tabwk0)
!
!====
! 4.  RESOLUTION DU PROBLEME M . A = F ET DE LA VITESSE STOCKEE
!           --- RESOLUTION AVEC FORCE1 COMME SECOND MEMBRE ---
!====
!
    r8bid = dt/2.d0
!
    do ieq = 1, neq
!
        acce1(ieq) = tabwk1(ieq)*zr(iforc1+ieq-1)
!
!        --- VITESSE AUX INSTANTS 'TEMPS + DT' ---
        vite2(ieq) = vite1(ieq)+r8bid*acce1(ieq)
!
    end do
!
!
!====
! 5.  CALCUL DES ENERGIES
!
!====
!
    if (ds_energy%l_comp) then
        masse1 = masse//'           '
        amort1 = amort//'           '
        rigid1 = rigid//'           '
        call wkvect('FNODABID', 'V V R', 2*neq, ifnobi)
        call wkvect('FCINEBID', 'V V R', 2*neq, ifcibi)
! ON CALCULE LA VITESSE A T N-1
        call wkvect('VIT0_TR', 'V V R', neq, ivit0r)
        do ieq = 1, neq
            zr(ivit0r-1+ieq) = vit0(ieq)+r8bid*acc0(ieq)
        end do
        call enerca(k19bid, dep0, zr(ivit0r), depl1, vite2, &
                    masse1, amort1, rigid1, fexte, famor, &
                    fliai, zr(ifnobi), zr(ifcibi), lamort, .true._1, &
                    .false._1, ds_energy, '&&DLDIFF')
        call jedetr('FNODABID')
        call jedetr('FCINEBID')
        call jedetr('VIT0_TR')
    end if
!====
! 5. TRANSFERT DES NOUVELLES VALEURS DANS LES ANCIENNES
!====
!
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, depl1, b_incx, dep0, b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vite1, b_incx, vit0, b_incy)
    b_n = to_blas_int(neq)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, acce1, b_incx, acc0, b_incy)
!
!====
! 7. ARCHIVAGE EVENTUEL DANS L'OBJET SOLUTION
!====
!
    if (archiv .eq. 1) then
!
        istoc = 0
        alarm = 1
!
        call dlarch(result, neq, istoc, iarchi, ' ', &
                    alarm, temps, nbtyar, typear, masse, &
                    depl1, vite2, acce1, fexte(neq+1), famor(neq+1), &
                    fliai(neq+1))
!
    end if
!===
! 8. ARCHIVAGE DES PARAMETRES
!===
    call nmarpc(ds_energy, numrep, temps)
!
end subroutine
