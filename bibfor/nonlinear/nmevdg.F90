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
subroutine nmevdg(sddisc, vale, i_echec, i_echec_acti)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/extdch.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbacce.h"
#include "asterfort/tbliva.h"
#include "asterfort/utdidt.h"
!
    integer(kind=8) :: i_echec, i_echec_acti
    character(len=19) :: sddisc, vale(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - EVENEMENTS)
!
! GESTION DE L'EVENEMENT DELTA_GRANDEUR
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization TEMPORELLE
! IN  VALE   : INCREMENTS DES VARIABLES
!               OP0070: VARIABLE CHAPEAU
!               OP0033: TABLE
! IN  i_echec : OCCURRENCE DE L'ECHEC
! OUT i_echec_acti : VAUT i_echec SI EVENEMENT DECLENCHE
!                   0 SINON
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv, ier
    integer(kind=8) :: deb, fin, etat_loca
    integer(kind=8), pointer:: loca(:) => null()
    real(kind=8) :: valref, dval, r8bid
    integer(kind=8) :: ibid
    character(len=8) :: k8bid, crit, typext
    complex(kind=8) :: c16bid
    character(len=16) :: nocham, nocmp
    parameter(typext='MAX_ABS')
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... DELTA_GRANDEUR'
    end if
!
! --- INITIALISATIONS
!
    i_echec_acti = 0
    r8bid = 0.d0
!
! --- PARAMETRES
!
    call utdidt('L', sddisc, 'ECHE', 'NOM_CHAM', index_=i_echec, &
                valk_=nocham)
    call utdidt('L', sddisc, 'ECHE', 'NOM_CMP', index_=i_echec, &
                valk_=nocmp)
    call utdidt('L', sddisc, 'ECHE', 'VALE_REF', index_=i_echec, &
                valr_=valref)
    call utdidt('L', sddisc, 'ECHE', 'CRIT_COMP', index_=i_echec, &
                valk_=crit)
!
! --- DVAL :MAX EN VALEUR ABSOLUE DU DELTA(CHAMP+CMP)
!
    if (vale(1) (1:8) .eq. '&&OP0033') then
!
!       RESULTAT DE CALC_POINT_MAT OP0033, FORMAT_TABLE='CMP_COLONNE',
        call tbacce(vale(1) (1:16), 1, nocmp, 'L', ibid, &
                    dval, c16bid, k8bid)
    else if (vale(1) (1:8) .eq. '&&OPB033') then
!
!       RESULTAT DE CALC_POINT_MAT OP0033, FORMAT_TABLE='CMP_LIGNE',
        call tbliva(vale(1) (1:16), 1, 'CMP', [ibid], [r8bid], &
                    [c16bid], nocmp, 'EGAL', [0.d0], 'VALEUR', &
                    k8bid, ibid, dval, c16bid, k8bid, &
                    ier)
        if (ier .ne. 0) then
            dval = 0.d0
        end if
    else
!
!       RESULTAT DE STAT_NON_LINE

        ! Extraction du filtre sur la liste des mailles
        call jeveuo(sddisc//'.ELOC', 'L', vi=loca)
        etat_loca = loca(SIZE_LELOCA*(i_echec-1)+1)

        if (etat_loca .eq. LOCA_VIDE) then
            dval = 0
        else if (etat_loca .eq. LOCA_PARTIEL) then
            deb = loca(SIZE_LELOCA*(i_echec-1)+2)
            fin = loca(SIZE_LELOCA*(i_echec-1)+3)
            call extdch(typext, vale, nocham, nocmp, dval, lst_loca=loca(deb:fin))
        else if (etat_loca .eq. LOCA_TOUT) then
            call extdch(typext, vale, nocham, nocmp, dval)
        else
            ASSERT(.false.)
        end if
    end if
!
    ASSERT(crit .eq. 'GT')
!
    if (dval .gt. valref) then
        i_echec_acti = i_echec
    end if
!
    call jedema()
end subroutine
