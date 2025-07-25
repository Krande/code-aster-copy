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
subroutine lc0034(fami, kpg, ksp, imate, &
                  carcri, epsm, &
                  deps, sigm, vim, option, angmas, &
                  sigp, vip, typmod, icomp, &
                  dsidep, codret)
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/nmhuj.h"
#include "asterfort/utlcal.h"
#include "asterfort/tecael.h"
!
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: imate
    real(kind=8) :: carcri(*)
    real(kind=8), intent(in) :: epsm(*)
    real(kind=8), intent(in) :: deps(*)
    real(kind=8), intent(in) :: sigm(6)
    real(kind=8) :: vim(50)
    character(len=16), intent(in) :: option
    real(kind=8), intent(in) :: angmas(3)
    real(kind=8), intent(out) :: sigp(6)
    real(kind=8) :: vip(50)
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: icomp
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
    real(kind=8)     :: npal, crit
    character(len=8) :: nomail
    integer(kind=8)          :: iadzi, iazk24, ndt, ndi
    aster_logical    :: debug, redec
! ----------------------------------------------------------------
    common/meshuj/debug
    common/tdim/ndt, ndi

! ----------------------------------------------------------------
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! HUJEUX
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHinteg         : parameters for integration of behaviour
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: algo_inte
!
! --------------------------------------------------------------------------------------------------
!
    call utlcal('VALE_NOM', algo_inte, carcri(6))
!
    call nmhuj(fami, kpg, ksp, typmod, imate, &
               carcri, &
               angmas, epsm, &
               deps, sigm, vim, option, sigp, &
               vip, dsidep, codret)

!
! M. Kham (18/07/2018) :: Interception de l'erreur
! ===================================================
!
! carcri(3):  RESI_INTE
! carcri(5):  ITER_INTE_RELA
!
! ALGO_INTE = | SPECIFIQUE
!             | SEMI_EXPLICITE
!             | BASCULE_EXPLICIT
!
!    debug=.true.
! -------------------------------------------------------------------------
!
!                          Gestion du cas codret = 0
!
! -------------------------------------------------------------------------
!
! Variable internes disponibles:
!
! XX  V32 : travail du second ordre: valeur instantannee
!         se debrouller pour le stocker en sortie
! XX  V33 : determinant de la matrice tangente
! ->  V34 : mecanismes actifs (inutile)
! ->  V35 : compteur d'iterations locales (inutile)
! -------------------------------------------------------------------
!
    crit = carcri(3)
    redec = (carcri(5) .lt. -1.) .or. (carcri(5) .gt. 1.)
!
    if (algo_inte(1:16) .eq. 'BASCULE_EXPLICIT' .and. (.not. redec)) then
!
! On remet codret a zero si critere non depasse
! et on sauvegarde un indicateur d'erreur dans V34 (realise dans nmhuj)
! -------------------------------------------------------------------
        if (codret .eq. 1) then

            if (vip(34) .gt. crit) then
                sigp(1:ndt) = sigm(1:ndt)
                codret = 2
            else
                codret = 0
            end if
!
            if (debug) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                write (6, *)
                write (6, '(A)') '!!!(o_o)!!! ---------------------------- !!!(o_o)!!!'
                write (6, '(A)') '!!!(o_o)!!! ATTENTION:                   !!!(o_o)!!!'
                write (6, '(A)') '!!!(o_o)!!! CV non atteinte a la maille  !!!(o_o)!!!'
                write (6, '(3A)') '!!!(o_o)!!! ', nomail, '                     !!!(o_o)!!!'
                write (6, '(A,I1,A)') '!!!(o_o)!!! CODRET    =', codret, &
                    '                 !!!(o_o)!!!'
                write (6, '(A,E12.5,A)') '!!!(o_o)!!! ERREUR    = ', vip(34), &
                    '     !!!(o_o)!!!'
                write (6, '(A)') '!!!(o_o)!!! ---------------------------- !!!(o_o)!!!'
            end if
        end if
!
    elseif (algo_inte(1:16) .eq. 'BASCULE_EXPLICIT' .and. icomp .eq. 3) then
!
! initialisation des variables internes utilisees dans le cas suivant
! -------------------------------------------------------------------
        vim(34) = 0.d0
        vim(35) = 0.d0
        vip(34) = 0.d0
        vip(35) = 0.d0
!
    elseif (algo_inte(1:16) .eq. 'BASCULE_EXPLICIT' .and. icomp .gt. 3) then
!
! npal = nombre d'iteration maximal pour icomp=4
        npal = 4.*abs(carcri(5))
!
! dans le cas codret=0, on incremente l'erreur sur le critere (V34)
! -----------------------------------------------------------------
        vip(34) = vim(34)+vip(34)
!
        if (codret .eq. 1) then
!
! On remet codret a zero et on sauvegarde un indicateur d'erreur dans V35
! V35: Nbre d'occurence de l'erreur sur nombre total de redecoupages
! ------------------------------------------------------------------------
            codret = 0
            vip(35) = vim(35)+100./npal
!
            if (debug) then
                call tecael(iadzi, iazk24)
                nomail = zk24(iazk24-1+3) (1:8)
                write (6, *)
                write (6, '(A)') '!!!(o_o)!!! ---------------------------- !!!(o_o)!!!'
                write (6, '(A)') '!!!(o_o)!!! ATTENTION:                   !!!(o_o)!!!'
                write (6, '(A)') '!!!(o_o)!!! CV non atteinte a la maille  !!!(o_o)!!!'
                write (6, '(3A)') '!!!(o_o)!!! ', nomail, '                     !!!(o_o)!!!'
                write (6, '(A,I1,A)') '!!!(o_o)!!! CODRET    =', &
                    codret, '                 !!!(o_o)!!!'
!             write(6,'(A,2(I4,A))') '!!!(o_o)!!! ITERATION =',nint(iter),' SUR ',nint(npal),&
!             '     !!!(o_o)!!!'
                write (6, '(A,E12.5,A)') '!!!(o_o)!!! ERREUR    = ', vip(34), '     !!!(o_o)!!!'
                write (6, '(A)') '!!!(o_o)!!! ---------------------------- !!!(o_o)!!!'
            end if
        end if
!
! evaluation de l'erreur en fin de redecoupage pour icomp=4
! l'erreur est calculee dans nmhuj et stockee dans la variable V34
! stockage du numero d'increment si on n'est pas au dernier pas
!
        if (vip(34) .gt. crit) then
            sigp(1:ndt) = sigm(1:ndt)
            codret = 2
!          codret = 1
        end if
! -------------------------------------------------------------------------
!
!                          Gestion du cas codret = 2
!
! -------------------------------------------------------------------------
!
    elseif (algo_inte(1:14) .eq. 'SEMI_EXPLICITE' .and. codret .eq. 1) then
!
        if (icomp .gt. 3) then
            codret = 2
        end if
!
        if (debug) then
            call tecael(iadzi, iazk24)
            nomail = zk24(iazk24-1+3) (1:8)
            write (6, *)
            write (6, '(A)') '!!!(o_o)!!! ---------------------------- !!!(o_o)!!!'
            write (6, '(A)') '!!!(o_o)!!! ATTENTION:                   !!!(o_o)!!!'
            write (6, '(A)') '!!!(o_o)!!! CV non atteinte a la maille  !!!(o_o)!!!'
            write (6, '(3A)') '!!!(o_o)!!! ', nomail, '                     !!!(o_o)!!!'
            write (6, '(A,I1,A)') '!!!(o_o)!!! CODRET =', codret, '                    !!!(o_o)!!!'
            write (6, '(A)') '!!!(o_o)!!! ---------------------------- !!!(o_o)!!!'
        end if
    end if
!    debug=.false.
!
end subroutine
