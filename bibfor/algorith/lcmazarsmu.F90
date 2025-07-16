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

subroutine lcmazarsmu(fami, kpg, ksp, dimepb, imate, model, epsm, &
                      deps, varim, option, sigp, varip, dsidep)
    implicit none
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/mazarsmu.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/rcexistvarc.h"
#include "asterfort/utmess.h"
#include "asterfort/get_varc.h"
!
    character(len=6)    :: model
    character(len=16)   :: option
    character(len=*)    :: fami
    integer(kind=8)             :: imate, kpg, ksp, dimepb
    real(kind=8)        :: epsm(*), deps(*), varim(*), sigp(*), varip(*), dsidep(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
!                 LOI DE COMPORTEMENT MAZARS UNILATÉRALE : MU_MODEL
!
! --------------------------------------------------------------------------------------------------
!
! IN
!       fami        : famille de point de gauss (rigi, mass, ...)
!       kpg         : numéro du point de gauss
!       ksp         : numéro du sous-point
!       imate       : adresse du matériau
!       epsm        : déformation à temps(-)
!       deps        : incrément de déformation
!       varim       : variables internes à temps(-)
!       option      : option demandée
!                       rigi_meca_tang ->     dsidep      --> rigi
!                       full_meca      -> sig dsidep vip  --> rigi  resi
!                       raph_meca      -> sig        vip  -->       resi
! OUT
!       sigp        : contraintes à temps(+)
!       varip       : variables internes à temps(+)
!       dsidep      : matrice tangente
!
!
! --------------------------------------------------------------------------------------------------
!
! Variables internes
!       1  -> icels  : critere sigma
!       2  -> icelu  : critere epsi
!       3  -> idomm  : endommagement
!       4  -> iepsqt : valeur de epseqp de traction
!       5  -> iepsqc : valeur de epseqp de compression
!       6  -> irsigm : facteur de triaxialite en contrainte
!       7  -> itemp  : temperature maximale atteinte par le materiau
!       8  -> idissd : dissipation d'endommagement
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter   :: nbval = 10
    integer(kind=8)             :: icodre(nbval)
    character(len=16)   :: nomres(nbval)
    character(len=8), parameter :: mazars(nbval) = ['EPSD0   ', 'K       ', 'AC      ', &
                                                    'BC      ', 'AT      ', 'BT      ', &
                                                    'SIGM_LIM', 'EPSI_LIM', 'EPSC0   ', &
                                                    'EPST0   ']
    real(kind=8)        :: valres(nbval+2)
!   Index des coefficients de la loi
!           iepsd0   =  1, ik       = 2, iac    = 3, ibc    =  4, iat = 5, ibt = 6
!           isigmlim =  7, iepsilim = 8, iepsc0 = 9, iepst0 = 10
!           iyoung   = 11, inu      = 12
    integer(kind=8), parameter :: iepsd0 = 1
    integer(kind=8), parameter :: isigmlim = 7, iepsilim = 8, iepsc0 = 9, iepst0 = 10
    integer(kind=8), parameter :: iyoung = 11, inu = 12
!
    integer(kind=8), parameter   :: itemp = 7, nbvarint = 8
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = 1.4142135623731d0
!
    integer(kind=8)             :: iret, nbres, dimloc
    real(kind=8)        :: epsthe, kdess, bendo, epsane, epsmldc(6), depsldc(6), ee, nu
    real(kind=8)        :: valpar, tempm, tempp, tempref, tempmax, valhydr, valsech, valsechref
    character(len=1)    :: poum
    character(len=8)    :: nompar
    character(len=30)   :: valkm(3)
    aster_logical       :: resi, rigi, istemper, ishydrat, issechag, iscplane
!
! --------------------------------------------------------------------------------------------------
!
    rigi = (option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RIGI')
    resi = (option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH')
    !
    iscplane = model .eq. 'C_PLAN'
    if (iscplane) then
        dimloc = 4
    else
        dimloc = 6
    end if
    !
    ! Récupération de la TEMPÉRATURE
    call get_varc(fami, kpg, ksp, 'T', tempm, tempp, tempref, istemper)
    !
    if (.not. istemper) then
        epsthe = 0.0d0
        tempm = 0.0d0
        tempp = 0.0d0
        tempref = 0.0d0
    end if
    !
    ! Y a-t-il de HYDR ou SECH : Par défaut c'est nul
    bendo = 0.0d0
    kdess = 0.0d0
    ! Valeur des champs : par défaut on considère qu'ils sont nuls
    valhydr = 0.0d0
    valsech = 0.0d0
    valsechref = 0.0d0
    !
    ishydrat = rcexistvarc('HYDR')
    issechag = rcexistvarc('SECH')
    if (ishydrat .or. issechag) then
        nomres(1) = 'B_ENDOGE'; nomres(2) = 'K_DESSIC'
        valres(1:2) = 0.0
        call rcvala(imate, ' ', 'ELAS', 0, ' ', [0.d0], 2, nomres, valres, icodre, 0)
        bendo = valres(1)
        kdess = valres(2)
        if ((icodre(1) .eq. 0) .and. (.not. ishydrat)) then
            valkm(1) = 'MAZARS_UNIL'
            valkm(2) = 'ELAS/B_ENDOGE'
            valkm(3) = 'HYDR'
            call utmess('F', 'COMPOR1_74', nk=3, valk=valkm)
        end if
        if ((icodre(2) .eq. 0) .and. (.not. issechag)) then
            valkm(1) = 'MAZARS_UNIL'
            valkm(2) = 'ELAS/K_DESSIC'
            valkm(3) = 'SECH'
            call utmess('F', 'COMPOR1_74', nk=3, valk=valkm)
        end if
    end if
    !
    ! Température maximale
    tempmax = max(tempp, varim(itemp))
    ! Mémorise ou pas la température maximale atteinte
    poum = '-'
    if (resi) then
        varip(itemp) = tempmax
        poum = '+'
    end if
    !
    ! Caractéristiques élastiques
    nomres(1) = 'E'; nomres(2) = 'NU'; nbres = 2
    if (istemper) then
        nbres = nbres+1
        nomres(nbres) = 'ALPHA'
    end if
    nompar = 'TEMP'; valpar = tempmax
    !
    call rcvalb(fami, kpg, ksp, poum, imate, ' ', 'ELAS', 1, nompar, [valpar], &
                nbres, nomres, valres, icodre, 1)
    !
    ee = valres(1)
    nu = valres(2)
    if (istemper) then
        if (icodre(3) .ne. 0) then
            call utmess('F', 'CALCULEL_15')
        end if
        epsthe = valres(3)*(tempp-tempref)
    end if
    !
    if (ishydrat) then
        call rcvarc('F', 'HYDR', poum, fami, kpg, ksp, valhydr, iret)
    end if
    if (issechag) then
        call rcvarc('F', 'SECH', poum, fami, kpg, ksp, valsech, iret)
        call rcvarc('F', 'SECH', 'REF', fami, kpg, ksp, valsechref, iret)
    end if
    ! Déformation sphérique anélastique
    epsane = epsthe+kdess*(valsech-valsechref)-bendo*valhydr
    !
    ! Paramètres matériau
    valres(:) = 0.0
    call rcvalb(fami, kpg, ksp, poum, imate, ' ', 'MAZARS', 1, nompar, [valpar], &
                nbval, mazars, valres, icodre, 0)
    if (icodre(isigmlim)+icodre(iepsilim) .ne. 0) then
        valkm(1) = 'MAZARS_UNIL'
        valkm(2) = mazars(isigmlim)
        valkm(3) = mazars(iepsilim)
        call utmess('F', 'COMPOR1_76', nk=3, valk=valkm)
    end if
    ! C'est soit iepsd0 soit iepst0 et iepsc0
    !     Dans le comportement iepsd0 ne sert plus
    ASSERT(icodre(iepst0) .ne. icodre(iepsd0))
    ASSERT(icodre(iepst0) .eq. icodre(iepsc0))
    if (icodre(iepsd0) .eq. 0) then
        valres(iepst0) = valres(iepsd0)
        valres(iepsc0) = valres(iepst0)/(nu*rac2)
    else
        valres(iepsd0) = valres(iepst0)
    end if
    ! Ajout de Young et nu dans valres
    valres(iyoung) = ee
    valres(inu) = nu
    !
    ! ----------------------------------------------------------------------------------------------
    !
    !   Loi de mazars unilatérale : mu_model
    !
    ! ----------------------------------------------------------------------------------------------
    ! Calcul de la déformation mécanique : Totale - thermique - hydratation - séchage
    depsldc(1:6) = 0.0; epsmldc(1:6) = 0.0
    if (resi) then
        varip(1:nbvarint) = varim(1:nbvarint)
        depsldc(1:dimloc) = deps(1:dimloc)
    end if
    epsmldc(1:3) = epsm(1:3)-epsane
    epsmldc(4) = epsm(4)
    if (.not. iscplane) then
        epsmldc(5) = epsm(5)
        epsmldc(6) = epsm(6)
    end if
    !
    call mazarsmu(option, epsmldc, depsldc, dimloc, valres, varim, varip, sigp, dsidep)
    !
    !  On retourne sigp*rac2.
    sigp(4) = rac2*sigp(4)
    if (.not. iscplane) then
        sigp(5) = rac2*sigp(5)
        sigp(6) = rac2*sigp(6)
    end if
end subroutine
