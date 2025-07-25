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
subroutine nmadat(sddisc, numins, nbiter, valinc)
!
    implicit none
!
#include "asterc/r8maem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/compr8.h"
#include "asterfort/diadap.h"
#include "asterfort/diinst.h"
#include "asterfort/getAdapAction.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/nmcadt.h"
#include "asterfort/nmdcei.h"
#include "asterfort/nmjalo.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "event_def.h"
#include "jeveux.h"
!
    character(len=19) :: valinc(*)
    character(len=19) :: sddisc
    integer(kind=8) :: numins, nbiter
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! GESTION DE L'ADAPTATION DE PAS DE TEMPS
!   CALCUL DE DTPLUS
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! IN  NUMINS : NUMERO D'INSTANT
! IN  NBITER : NOMBRE D'ITERATIONS DE NEWTON
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: dtplus = '&&NMADAP.DTPLUS'
    integer(kind=8) :: nb_adap, i_adap, jdt
    character(len=19) :: metlis
    real(kind=8) :: r8bid, dt, min, pasmin, pasmax, dtm, jalon
    real(kind=8) :: newins, newdt, deltac
    real(kind=8) :: inst, prec, valr(2)
    real(kind=8) :: insfin, insref
    aster_logical :: ladap, uncrok
    character(len=24) :: tpsite
    integer(kind=8) :: jiter
    integer(kind=8) :: nb_inst, nmax, inspas
    character(len=24) :: tpsext
    integer(kind=8) :: jtpsex, actionType
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Get parameters
    call utdidt('L', sddisc, 'LIST', 'PAS_MINI', &
                valr_=pasmin)
    call utdidt('L', sddisc, 'LIST', 'PAS_MAXI', &
                valr_=pasmax)
    call utdidt('L', sddisc, 'LIST', 'METHODE', &
                valk_=metlis)
    call utdidt('L', sddisc, 'LIST', 'NBINST', &
                vali_=nb_inst)
    call utdidt('L', sddisc, 'LIST', 'NB_PAS_MAXI', &
                vali_=nmax)
!
! --- NOM SDS DE LA SDDISC
!
    tpsite = sddisc(1:19)//'.ITER'
!
! --- PRECISION SUR LES INSTANTS
! --- (LIEE A CELLE DE VAL_MIN DE PAS_MINI DANS DEFI_LIST_INST.CAPY)
!
    prec = 1.d-12
!
! --- ON NE FAIT DE L'ADAPTATION DE PAS DE TEMPS QU'EN GESTION AUTO
!
    if (metlis .eq. 'MANUEL') then
        goto 999
    end if
!
! --- INSTANT COURANT
!
    inst = diinst(sddisc, numins)
!
! --- PROCHAIN INSTANT DE PASSAGE OBLIGATOIRE (JALON) ?
!
    call nmjalo(sddisc, inst, prec, jalon)
    if (jalon .eq. r8vide()) goto 999
!
! --- NOMBRE DE SCHEMAS D'ADAPTATION : NADAPT
!
    call utdidt('L', sddisc, 'LIST', 'NADAPT', vali_=nb_adap)
!
! --- LISTE DES NADAPT PAS DE TEMPS POSSIBLES
!
    call wkvect(dtplus, 'V V R', nb_adap, jdt)
!
! --- PAS DE TEMPS PAR DEFAUT (LE DERNIER, SAUF SI JALON) : DTM
!
    call utdidt('L', sddisc, 'LIST', 'DT-', &
                valr_=dtm)
!
! --- STOCKAGE DU NOMBRE D'ITERATIONS DE NEWTON ET EXTENSION
!
    call jeveuo(tpsite, 'L', jiter)
    zi(jiter-1+numins) = nbiter
    call juveca(tpsite, nb_inst+1)
!
! ----------------------------------------------------------------------
!    CALCUL DU PAS DE TEMPS
! ----------------------------------------------------------------------
!
    call utmess('I', 'ADAPTATION_1')
!
    do i_adap = 1, nb_adap
!
        call getAdapAction(sddisc, i_adap, actionType)
!
        zr(jdt-1+i_adap) = r8vide()
!
! ----- DOIT-ON ADAPTER ?
!
        ladap = diadap(sddisc, i_adap)
        if (ladap) then
            call nmcadt(sddisc, i_adap, numins, valinc, zr(jdt-1+i_adap))
        end if
        newdt = zr(jdt-1+i_adap)
!
! ----- AFFICHAGE
!
        if (newdt .ne. r8vide()) then
            call utmess('I', 'ADAPTATION_2', sk=adapActionKeyword(actionType), sr=newdt)
        else
            call utmess('I', 'ADAPTATION_3', sk=adapActionKeyword(actionType))
        end if
    end do
!
! --- ON CHOISIT LE PLUS PETIT DT PARMI LES NADAPT PAS DE TEMPS
! --- POSSIBLES
! --- SI AUCUN CRITERE N'EST VERIFIE, ON PREND LE PAS DE TEMPS "MOINS"
!
    dt = r8maem()
    uncrok = .false.
    do i_adap = 1, nb_adap
        newdt = zr(jdt-1+i_adap)
        if (newdt .ne. r8vide()) then
            dt = min(dt, newdt)
            uncrok = .true.
        end if
    end do
!
    if (uncrok) then
        call utmess('I', 'ADAPTATION_5', sr=dt)
    else
        dt = dtm
        call utmess('I', 'ADAPTATION_4', sr=dt)
    end if
!
! --- PROJECTION SUR LA BORNE SUP (POUR TOUTES LES METHODES)
!
    if (dt .gt. pasmax) then
!       EMISSION DU MESSAGE D'INFO (SAUF POUR IMPLEX)
        if (actionType .ne. ADAP_ACT_IMPLEX) then
            valr(1) = dt
            valr(2) = pasmax
            call utmess('I', 'ADAPTATION_12', nr=2, valr=valr)
        end if
        dt = pasmax
    end if
!
! --- PROJECTION SUR LA BORNE INF POUR IMPLEX
! --- (ATTENTION : A FAIRE AVANT L'AJUSTEMENT / JALON)
!
    if (actionType .eq. ADAP_ACT_IMPLEX) then
        if (dt .lt. pasmin) dt = pasmin
    end if
!
! --- LA DECOUPE DU PAS DE TEMPS PEUT DONNER UN DELTAT MAXI A RESPECTER
!
    tpsext = sddisc(1:19)//'.AEXT'
    call jeveuo(tpsext, 'E', jtpsex)
    insref = zr(jtpsex-1+1)
    deltac = zr(jtpsex-1+2)
    insfin = zr(jtpsex-1+3)
    if (insref .ne. r8vide()) then
        if (inst .le. insfin) then
            dt = deltac
            call utmess('I', 'ADAPTATION_10', sr=dt)
        else
            zr(jtpsex-1+1) = r8vide()
        end if
    end if
!
! --- AJUSTEMENT DE DT EN FONCTION DU PROCHAIN JALON
!
    if (compr8(inst+dt, 'GT', jalon, prec, 1)) then
!       LE NOUVEAU PAS DEPASSE LE PROCHAIN IPO :
!       ON FORCE A Y PASSER ET ON N'ENREGISTRE PAS DT
        dt = jalon-inst
    else if (compr8(inst+dt, 'GT', jalon-pasmin, prec, 1)) then
!       NOUVEAU DE PAS INFERIEUR A JALON, MAIS TROP PROCHE DE JALON :
!       ON FORCE A Y PASSER ET ON ENREGISTRE DT
        dt = jalon-inst
        call utdidt('E', sddisc, 'LIST', 'DT-', &
                    valr_=dt)
    else
!       NOUVEAU PAS DE TEMPS OK
!       ON ENREGISTRE DT
        call utdidt('E', sddisc, 'LIST', 'DT-', &
                    valr_=dt)
    end if
!
    call utmess('I', 'ADAPTATION_6', sr=dt)
!
! --- ON VERIFIE LES GARDE FOUS
!
    if (actionType .ne. ADAP_ACT_IMPLEX) then
        if (dt .lt. pasmin) then
            call utmess('F', 'ADAPTATION_11', sr=dt)
        end if
    end if

!
! --- INSERTION DU NOUVEL INSTANT
!
    inspas = 1
    newins = inst+dt
    call nmdcei(sddisc, numins, [newins], nb_inst, inspas, &
                'ADAP', r8bid)
!
999 continue
!
    call jedetr(dtplus)
    call jedema()
end subroutine
