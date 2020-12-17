! --------------------------------------------------------------------
! Copyright (C) 1991 - 2020 - EDF R&D - www.code-aster.org
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
! person_in_charge: tanguy.mathieu at edf.fr
!
subroutine op0027()
!
use calcG_type
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/ccbcop.h"
#include "asterfort/cgComputeGtheta.h"
#include "asterfort/cgComputeTheta.h"
#include "asterfort/cgTableG.h"
#include "asterfort/cgExportTableG.h"
#include "asterfort/cgVerification.h"
#include "asterfort/deprecated_algom.h"
#include "asterfort/detrsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "jeveux.h"
#include "asterfort/rsmena.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsrusd.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/utmess.h"
#include "asterfort/xcourb.h"
#include "asterfort/wkvect.h"
! --------------------------------------------------------------------------------------------------
!
!      OPERATEUR :     CALC_H
!
!      BUT:CALCUL DU TAUX DE RESTITUTION D'ENERGIE PAR LA METHODE THETA
!          CALCUL DES FACTEURS D'INTENSITE DE CONTRAINTES
!
!---------------------------------------------------------------------------------------------------
!
    type(CalcG_field) :: cgField
    type(CalcG_theta) :: cgTheta
    type(CalcG_study) :: cgStudy
    type(CalcG_InfoTe) :: cgInfoTe
!
    integer           :: iopt, istore, nume_ordre, iret, ipuls, jinst
    integer           :: jopt, nbropt
    real(kind=8)      :: time, puls, gth(4)
    character(len=8)  :: option, k8bid, resuc2
    character(len=19) :: lisopt
    character(len=24) :: depla, chvite, chacce, basloc, courb
    aster_logical     :: lmoda, exitim
!
    integer, pointer :: ordr(:) => null()
!---------------------------------------------------------------------------------------------------
    call jemarq()
    call infmaj()
    call deprecated_algom('CALC_H')
!
! Fiches concernées par le chantier (A supprimer à la fin)
! A Faire: #29573, #27931, #29703, #30288
!
!-- Initialisation des champs et des paramètres
    call cgField%initialize()
    call cgTheta%initialize()
    call cgInfoTe%initialize()
    call cgStudy%initialize(cgField%result_in, cgField%list_nume(1))
!
    courb = '&&OP0027.COURB' 
    puls  =0.0d0   
    lmoda = .false.
    exitim = .false.
!
!-- Calcul de la courbure 
    if (cgField%ndim .eq. 3) then
        basloc=cgTheta%crack//'.BASLOC'
        call xcourb(basloc, cgStudy%mesh, cgStudy%model, courb)
    endif
!
!-- Verification (A nettoyer)
    call cgVerification(cgField, cgTheta)
!
!-- Compute Theta factors
    call cgComputeTheta(cgField, cgTheta)
!
!-- Création table G
    call tbcrsd(cgField%table_g, 'G')
!
!-- ELAS INCR
    if (cgField%l_incr) then

        lisopt = '&&OP0027.LISOPT' 
        nbropt = 2
!
        call wkvect(lisopt, 'V V K16', nbropt, jopt)
        zk16(jopt) = 'VARI_ELNO'
        zk16(jopt+1) = 'EPSP_ELNO'
!
        call ccbcop(cgField%result_in, resuc2, cgField%list_nume_name,&
                    cgField%nb_nume, lisopt, nbropt)
    endif
!
!-- Loop on option
    do iopt = 1, cgField%nb_option
!
        option = cgField%list_option(iopt)
        call cgStudy%setOption(option)
!
        if (cgField%isModeMeca()) then
            if (option .eq. 'K') then
                lmoda = .true.
            else
                call utmess('F', 'RUPTURE0_27')
            endif
        endif
!
        do istore = 1, cgField%nb_nume
!
            nume_ordre = cgField%list_nume(istore)
!
            call cgStudy%initialize(cgField%result_in, nume_ordre)
!
! --------  Maillage similaire sd_fond_fissure et sd_resu
            ASSERT(cgTheta%mesh == cgStudy%mesh)
!            
! --------  Récupération des champs utiles pour l'appel à calcul
            call rsexch('F', cgField%result_in, 'DEPL', nume_ordre, depla, iret)
            call rsexch(' ', cgField%result_in, 'VITE', nume_ordre, chvite, iret)
            if (iret .ne. 0) then
                chvite = ' '
            else
                call rsexch(' ', cgField%result_in, 'ACCE', nume_ordre, chacce, iret)
            endif
!
            if (lmoda) then
                call rsadpa(cgField%result_in, 'L', 1, 'OMEGA2', nume_ordre,&
                            0, sjv=ipuls, styp=k8bid)
                puls = zr(ipuls)
                puls = sqrt(puls)
                time = 0.d0
            else
                call rsadpa(cgField%result_in, 'L', 1, 'INST', nume_ordre,&
                            0, sjv=jinst, styp=k8bid)
                time = zr(jinst)
                exitim = .true.
            endif      
!
!---------- Calcul de G(theta) pour les éléments 2D/3D option G et K
            call cgComputeGtheta(cgField, cgTheta, cgStudy, nume_ordre, depla, &
                                 chvite, chacce, time, courb, option, puls, &
                                 lmoda, gth)         
!
!---------- Création de la table de G et des SIFS
            call cgTableG(cgField, cgTheta, nume_ordre, option, time, lmoda, gth)
!
        end do
!
!------ Print fields 
        call cgTheta%print()
        call cgField%print()
!
    end do
!
    if (cgField%l_incr) then
!
        call jeexin(resuc2//'           .ORDR', iret)
        if (iret .ne. 0) then
            call jeveuo(resuc2//'           .ORDR', 'L', vi=ordr)
            call rsrusd(resuc2, ordr(1))
            call detrsd('RESULTAT', resuc2)
        endif
!
        call jedetr(cgField%list_nume_name)
        call jedetr(resuc2)
        call rsmena(cgField%result_in)
    endif
!
!-- Création de la table container 
    call cgExportTableG(cgField, cgTheta)

    call jedema()
!
end subroutine
