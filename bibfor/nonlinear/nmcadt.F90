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
subroutine nmcadt(sddisc, i_adap, nume_inst, hval_incr, dtp)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8vide.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/extdch.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
#include "asterfort/getAdapAction.h"
!
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: i_adap
    integer(kind=8), intent(in) :: nume_inst
    character(len=19), intent(in) :: hval_incr(*)
    real(kind=8), intent(out) :: dtp
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! CALCUL DU NOUVEAU PAS DE TEMPS EN CAS D'ADAPTATION
!
! ----------------------------------------------------------------------
!
!
! In  sddisc           : datastructure for time discretization
! IN  IADAPT : NUMERO DE LA METHODE D ADAPTATION TRAITEE
! IN  NUMINS : NUMERO D'INSTANT
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! OUT DTP    : NOUVEAU PAS DE TEMPS (DT+)
!
!
!
!
    integer(kind=8) :: deb, fin, etat_loca
    integer(kind=8), pointer:: loca(:) => null()

    integer(kind=8) :: nit, nbiter, action_type
    real(kind=8) :: dtm, pcent, valref, dval
    character(len=8) :: typext
    character(len=16) :: nocham, nocmp
    real(kind=8) :: eta, etad
    character(len=24) :: tpsite
    integer(kind=8) :: jiter
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! --- ACCES SD
!
    tpsite = sddisc(1:19)//'.ITER'
    call jeveuo(tpsite, 'L', jiter)
!
! --- METHODE DE CALCUL DE DT+
!
    call getAdapAction(sddisc, i_adap, action_type)
!
! --- PAS DE TEMPS PAR DEFAUT (LE DERNIER, SAUF SI JALON)
!
    call utdidt('L', sddisc, 'LIST', 'DT-', valr_=dtm)
!
!     ------------------------------------------------------------------
    if (action_type .eq. ADAP_ACT_FIXE) then
!     ------------------------------------------------------------------
!
        call utdidt('L', sddisc, 'ADAP', 'PCENT_AUGM', index_=i_adap, &
                    valr_=pcent)
        dtp = dtm*(1.d0+pcent/100.d0)
!
!     ------------------------------------------------------------------
    else if (action_type .eq. ADAP_ACT_INCR_QUANT) then
!     ------------------------------------------------------------------
!
        call utdidt('L', sddisc, 'ADAP', 'NOM_CHAM', index_=i_adap, &
                    valk_=nocham)
        call utdidt('L', sddisc, 'ADAP', 'NOM_CMP', index_=i_adap, &
                    valk_=nocmp)
        call utdidt('L', sddisc, 'ADAP', 'VALE_REF', index_=i_adap, &
                    valr_=valref)
        typext = 'MAX_ABS'
!
! ----- CALCUL DE C = MIN (VREF / |DELTA(CHAMP+CMP)| )
! -----             = VREF / MAX ( |DELTA(CHAMP+CMP)| )
! ----- DVAL :MAX EN VALEUR ABSOLUE DU DELTA(CHAMP+CMP)
!
        ! Extraction du filtre sur la liste des mailles
        call jeveuo(sddisc//'.ALOC', 'L', vi=loca)
        etat_loca = loca(SIZE_LALOCA*(i_adap-1)+1)

        if (etat_loca .eq. LOCA_VIDE) then
            dval = 0
        else if (etat_loca .eq. LOCA_PARTIEL) then
            deb = loca(SIZE_LALOCA*(i_adap-1)+2)
            fin = loca(SIZE_LALOCA*(i_adap-1)+3)
            call extdch(typext, hval_incr, nocham, nocmp, dval, lst_loca=loca(deb:fin))
        else if (etat_loca .eq. LOCA_TOUT) then
            call extdch(typext, hval_incr, nocham, nocmp, dval)
        else
            ASSERT(.false.)
        end if
!
! ----- LE CHAMP DE VARIATION EST IDENTIQUEMENT NUL : ON SORT
!
        if (dval .eq. 0.d0) then
            dtp = r8vide()
        else
            dtp = dtm*valref/dval
        end if
!
!     ------------------------------------------------------------------
    else if (action_type .eq. ADAP_ACT_ITER) then
!     ------------------------------------------------------------------
!
        call utdidt('L', sddisc, 'ADAP', 'NB_ITER_NEWTON_REF', index_=i_adap, &
                    vali_=nit)
        nbiter = zi(jiter-1+nume_inst)
        dtp = dtm*sqrt(dble(nit)/dble(nbiter+1))
!
!     ------------------------------------------------------------------
    else if (action_type .eq. ADAP_ACT_IMPLEX) then
!     ------------------------------------------------------------------
!
! ----- FACTEUR D'ACCELERATION ETA
!
        eta = 1.2d0
!
! ----- FACTEUR DE DECELERATION ETAD
!
        etad = 0.5d0
        typext = 'MIN_VAR'
!
! ----- CALCUL DE C = MIN (VREF / |DELTA(CHAMP+CMP)| )
!                   = VREF / MAX ( |DELTA(CHAMP+CMP)| )
! ----- DVAL :MAX EN VALEUR ABSOLUE DU DELTA(CHAMP+CMP)
!
        nocham = 'VARI_ELGA'
        nocmp = 'V1'
        call extdch(typext, hval_incr, nocham, nocmp, dval)
!
! ----- LE CHAMP DE VARIATION EST IDENTIQUEMENT NUL : ON SORT
!
        if (dval .ge. r8maem()) dval = eta
        dtp = dtm*dval
!
! ----- ON IMPOSE QUE LE DT SOIT COMPRIS ENTRE ETAD*DTM ET ETA*DTM
!
        if (dtp/dtm .ge. eta) then
            dtp = eta*dtm
        else if (dtp/dtm .le. etad) then
            dtp = etad*dtm
        end if
!
    else
        ASSERT(.false.)
    end if
!
    call asmpi_comm_vect('MPI_MIN', 'R', scr=dtp)
!
    call jedema()
end subroutine
