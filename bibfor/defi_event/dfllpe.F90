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
subroutine dfllpe(keywf, i_fail, event_typek, &
                  vale_ref, nom_cham, nom_cmp, crit_cmp, lst_loca, &
                  etat_loca, pene_maxi, resi_glob_maxi)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/dismoi.h"
#include "asterfort/reliem.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: keywf
    integer(kind=8), intent(in) :: i_fail
    character(len=16), intent(in) :: event_typek
    real(kind=8), intent(out) :: vale_ref
    character(len=16), intent(out) :: nom_cham
    character(len=16), intent(out) :: nom_cmp
    character(len=16), intent(out) :: crit_cmp
    character(len=24), intent(out) :: lst_loca
    integer(kind=8), intent(out):: etat_loca
    real(kind=8), intent(out) :: pene_maxi
    real(kind=8), intent(out) :: resi_glob_maxi
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_LIST_INST - Read parameters
!
! Get parameters of EVENEMENT for current failure keyword
!
! --------------------------------------------------------------------------------------------------
!
! In  keywf            : factor keyword to read failures
! In  i_fail           : index of current factor keyword to read failure
! In  event_typek      : type of event
! Out vale_ref         : value of VALE_REF for EVENEMENT=DELTA_GRANDEUR
! Out nom_cham         : value of NOM_CHAM for EVENEMENT=DELTA_GRANDEUR
! Out nom_cmp          : value of NOM_CMP for EVENEMENT=DELTA_GRANDEUR
! Out crit_cmp         : value of CRIT_CMP for EVENEMENT=DELTA_GRANDEUR
! Out lst_loca         : vecteur jeveux contenant la liste des mailles si DELTA_GRANDEUR
! out etat_loca        : en lien avec lst_loca 0=vide, 1=partiel (cf. lst_loca), 2=tout
! Out pene_maxi        : value of PENE_MAXI for EVENEMENT=INTERPENETRATION
! Out resi_glob_maxi   : value of RESI_GLOB_MAXI for EVENEMENT=RESI_MAXI
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nocc
    integer(kind=8)           :: nb_loca
    character(len=8)  :: mesh
    character(len=24) :: model
!
! --------------------------------------------------------------------------------------------------
!
    pene_maxi = 0.d0
    vale_ref = 0.d0
    resi_glob_maxi = 0.d0
    nom_cham = ' '
    nom_cmp = ' '
    crit_cmp = ' '
    etat_loca = 0
    lst_loca = ' '
!
! - Read parameters
!
    if (event_typek .eq. failEventKeyword(FAIL_EVT_INCR_QUANT)) then
        call getvr8(keywf, 'VALE_REF', iocc=i_fail, scal=vale_ref, nbret=nocc)
        ASSERT(nocc .gt. 0)
        call getvtx(keywf, 'NOM_CHAM', iocc=i_fail, scal=nom_cham, nbret=nocc)
        ASSERT(nocc .gt. 0)
        call getvtx(keywf, 'NOM_CMP', iocc=i_fail, scal=nom_cmp, nbret=nocc)
        ASSERT(nocc .gt. 0)
        crit_cmp = 'GT'

        if (nom_cham .eq. 'DEPL') then
            call getvtx(keywf, 'GROUP_NO', iocc=i_fail, nbret=nocc)
        else if (nom_cham .eq. 'SIEF_ELGA' .or. nom_cham .eq. 'VARI_ELGA') then
            call getvtx(keywf, 'GROUP_MA', iocc=i_fail, nbret=nocc)
        else
            ASSERT(.false.)
        end if
        etat_loca = merge(LOCA_TOUT, LOCA_PARTIEL, nocc .eq. 0)

        if (etat_loca .eq. LOCA_PARTIEL) then
            call getvid(' ', 'MODELE', scal=model, nbret=nocc)
            if (nocc .ne. 1) call utmess('F', 'LISTINST_4')
            call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
            write (lst_loca, '(A19,I5.5)') '&&OP0028.ECHE.LOCA.', i_fail
            if (nom_cham .eq. 'DEPL') then
                call reliem(model, mesh, 'NU_NOEUD', keywf, i_fail, 1, ['GROUP_NO'], &
                            ['GROUP_NO'], lst_loca, nb_loca)
            else if (nom_cham .eq. 'SIEF_ELGA' .or. nom_cham .eq. 'VARI_ELGA') then
                call reliem(model, mesh, 'NU_MAILLE', keywf, i_fail, 1, ['GROUP_MA'], &
                            ['GROUP_MA'], lst_loca, nb_loca)
            else
                ASSERT(.false.)
            end if
            if (nb_loca .eq. 0) etat_loca = LOCA_VIDE
        end if

    else if (event_typek .eq. failEventKeyword(FAIL_EVT_INTERPENE)) then
        call getvr8(keywf, 'PENE_MAXI', iocc=i_fail, scal=pene_maxi, nbret=nocc)
        ASSERT(nocc .gt. 0)
    else if (event_typek .eq. failEventKeyword(FAIL_EVT_RESI_MAXI)) then
        call getvr8(keywf, 'RESI_GLOB_MAXI', iocc=i_fail, scal=resi_glob_maxi, nbret=nocc)
        ASSERT(nocc .gt. 0)
    end if
!
end subroutine
