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
subroutine dltali(neq, result, imat, masse, rigid, &
                  liad, lifo, nchar, nveca, lcrea, &
                  lprem, lamort, t0, mate, mateco, &
                  carele, charge, infoch, fomult, modele, &
                  numedd, nume, solveu, criter, dep0, &
                  vit0, acc0, fexte0, famor0, fliai0, &
                  tabwk, force0, force1, ds_energy, kineLoad)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ajlagr.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dlfdyn.h"
#include "asterfort/dlfext.h"
#include "asterfort/dltini.h"
#include "asterfort/getvid.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mrmult.h"
#include "asterfort/preres.h"
#include "asterfort/resoud.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "blas/dcopy.h"
!
    character(len=8), intent(in) :: result
    type(NL_DS_Energy), intent(out) :: ds_energy
!
! --------------------------------------------------------------------------------------------------
!
!       DYNAMIQUE LINEAIRE TRANSITOIRE - ALGORITHME - INITIALISATION
!
! --------------------------------------------------------------------------------------------------
!
!  IN  : NEQ       : NOMBRE D'EQUATIONS
!  IN  : IMAT      : TABLEAU D'ADRESSES POUR LES MATRICES
!  IN  : MASSE     : MATRICE DE MASSE
!  IN  : RIGID     : MATRICE DE RIGIDITE
!  IN  : LIAD      : LISTE DES ADRESSES DES VECTEURS CHARGEMENT (NVECT)
!  IN  : LIFO      : LISTE DES NOMS DES FONCTIONS EVOLUTION (NVECT)
!  IN  : NCHAR     : NOMBRE D'OCCURENCES DU MOT CLE CHARGE
!  IN  : NVECA     : NOMBRE D'OCCURENCES DU MOT CLE VECT_ASSE
!  IN  : LCREA     : LOGIQUE INDIQUANT SI IL Y A REPRISE
!  IN  : LAMORT    : LOGIQUE INDIQUANT SI IL Y A AMORTISSEMENT
!  IN  : MATE      : NOM DU CHAMP DE MATERIAU
!  IN  : CARELE    : CARACTERISTIQUES DES POUTRES ET COQUES
!  IN  : CHARGE    : LISTE DES CHARGES
!  IN  : INFOCH    : INFO SUR LES CHARGES
!  IN  : FOMULT    : LISTE DES FONC_MULT ASSOCIES A DES CHARGES
!  IN  : MODELE    : MODELE
!  IN  : NUMEDD    : NUME_DDL DE LA MATR_ASSE RIGID
!  IN  : NUME      : NUMERO D'ORDRE DE REPRISE
!  IN  : SOLVEU    : NOM DU SOLVEUR
!  VAR : DEP0      : TABLEAU DES DEPLACEMENTS A L'INSTANT N
!  VAR : VIT0      : TABLEAU DES VITESSES A L'INSTANT N
!  VAR : ACC0      : TABLEAU DES ACCELERATIONS A L'INSTANT N
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: neq
    integer(kind=8) :: nveca, nchar
    integer(kind=8) :: liad(*)
    integer(kind=8) :: imat(3), nume
    real(kind=8) :: dep0(*), vit0(*), acc0(*)
    real(kind=8) :: fexte0(*), famor0(*), fliai0(*)
    real(kind=8) :: t0
    real(kind=8) :: tabwk(*)
    character(len=8) :: masse, rigid
    character(len=19) :: solveu
    character(len=24) :: charge, infoch, fomult, mate, mateco, carele
    character(len=24) :: modele, numedd
    character(len=24) :: lifo(*)
    character(len=24) :: criter, kineLoad
    character(len=19) :: force0, force1
    aster_logical :: lcrea, lprem
    aster_logical :: lamort
    complex(kind=8) :: cbid
    integer(kind=8) :: inchac
    integer(kind=8) :: ibid, icode, ieq, ndy, ifextm, ifextc
    character(len=8) :: matrei, maprei, dyna
    character(len=19) :: chsol
    integer(kind=8) :: iforc0, iforc1
    integer(kind=8) :: iret
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
!
! --------------------------------------------------------------------------------------------------
!
    call vtcreb(force0, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
    call jeveuo(force0(1:19)//'.VALE', 'E', iforc0)
    call vtcreb(force1, 'V', 'R', nume_ddlz=numedd, nb_equa_outz=neq)
    call jeveuo(force1(1:19)//'.VALE', 'E', iforc1)
!
! 1.2. ==> NOM DES STRUCTURES DE TRAVAIL
!
    chsol = '&&DLTALI.SOLUTION'
    maprei = ' '
!
!====
! 2. L'INITIALISATION
!====
!
    inchac = 0
    lcrea = .true.
    call dltini(lcrea, nume, result, dep0, vit0, &
                acc0, fexte0, famor0, fliai0, neq, &
                numedd, inchac, ds_energy)
!
!
!====
! 4. --- CHARGEMENT A L'INSTANT INITIAL OU DE REPRISE ---
!====
!
    call dlfext(nveca, nchar, t0, neq, liad, &
                lifo, charge, infoch, fomult, modele, &
                mate, mateco, carele, numedd, zr(iforc0))
!
!====
! 5. --- CALCUL DU CHAMP D'ACCELERATION INITIAL ---
!====
!
    if (inchac .ne. 0) then
!
! 5.1. ==> --- RESOLUTION AVEC FORCE1 COMME SECOND MEMBRE ---
!
        call jeveuo(force1(1:19)//'.VALE', 'E', iforc1)
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(iforc0), b_incx, zr(iforc1), b_incy)
        call dlfdyn(imat(1), imat(3), lamort, neq, dep0, &
                    vit0, zr(iforc1), tabwk)
!
        matrei = '&&MASSI'
        if (lprem) then
            lprem = .false.
            call ajlagr(rigid, masse, matrei)
!
! 5.2. ==> DECOMPOSITION OU CALCUL DE LA MATRICE DE PRECONDITIONEMENT
            call preres(solveu, 'V', icode, maprei, matrei, &
                        ibid, -9999)
        end if
!                                       ..          .
! 5.3. ==> RESOLUTION DU PROBLEME:  M.X  =  F - C.X - K.X
!                                       ..          .
!
        call resoud(matrei, maprei, solveu, kineLoad, 0, &
                    force1, chsol, 'V', [0.d0], [cbid], &
                    criter, .true._1, 0, iret)
!
! 5.4. ==> SAUVEGARDE DU CHAMP SOLUTION CHSOL DANS VDEPL
!
        call copisd('CHAMP_GD', 'V', chsol(1:19), force1(1:19))
        call jeveuo(force1(1:19)//'.VALE', 'L', iforc1)
!
! 5.5. ==> DESTRUCTION DU CHAMP SOLUTION CHSOL
!
        call detrsd('CHAMP_GD', chsol)
!
! 5.6 ==> STOCKAGE DE LA SOLUTION, FORC1, DANS LA STRUCTURE DE RESULTAT
!           EN TANT QUE CHAMP D'ACCELERATION A L'INSTANT COURANT
        b_n = to_blas_int(neq)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(iforc1), b_incx, acc0, b_incy)
!
    end if
!
! CALCUL DE LA FORCE INITIALE SI PAS DE REPRISE A PARTIR D UN RESULTAT
!
!
    call getvid('ETAT_INIT', 'RESULTAT', iocc=1, scal=dyna, nbret=ndy)
    if (ndy .eq. 0) then
        call wkvect('FEXT0M', 'V V R', neq, ifextm)
        call mrmult('ZERO', imat(1), dep0, fexte0, 1, &
                    .true._1)
        call mrmult('ZERO', imat(2), acc0, zr(ifextm), 1, &
                    .true._1)
        call wkvect('FEXT0C', 'V V R', neq, ifextc)
        if (lamort) then
            call mrmult('ZERO', imat(3), vit0, zr(ifextc), 1, &
                        .true._1)
        end if
        do ieq = 1, neq
            fexte0(ieq) = fexte0(ieq)+zr(ifextm-1+ieq)+zr(ifextc-1+ieq)
        end do
    end if
    call jedetr('FEXT0M')
    call jedetr('FEXT0C')
!
    if (ds_energy%l_comp .and. kineLoad .ne. ' ') then
        call utmess('F', 'DYNALINE2_11')
    end if
!
end subroutine
