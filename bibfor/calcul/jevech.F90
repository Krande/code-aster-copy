! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine jevech(nmparz, louez, itab)
!
    use calcul_module, only : ca_caindz_, ca_capoiz_, ca_iaoppa_, ca_iawlo2_, ca_iawloc_, ca_iel_, ca_igr_, ca_nbgr_, ca_nomte_, ca_nparin_, ca_npario_, ca_option_
!
    implicit none
!
! person_in_charge: jacques.pellet at edf.fr
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/chloet.h"
#include "asterfort/contex_param.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: nmparz, louez
    character(len=8) :: nompar, nommai
    character(len=1) :: loue
    integer :: itab
!-----------------------------------------------------------------
!  entrees:
!     nompar  : nom du parametre de l'option
!     louez   : 'L' ou 'E'  ( lecture/ecriture )
!
!  sorties:
!     itab     : adresse du champ local correspondant a nompar
!-----------------------------------------------------------------
    integer :: iachlo
    integer :: ilchlo, k, kk, debugr
    integer :: iparg, lgcata
    integer :: jceld, adiel
    integer :: debgr2, lonchl, decael, iadzi, iazk24
    integer :: opt, iaopd2, iaoplo, iapara, ipara, npari2
    aster_logical :: etendu
    character(len=24) :: valk(5)
! ------------------------------------------------------------------
    nompar = nmparz
    loue = louez
!
    ASSERT(loue.eq.'L' .or. loue.eq.'E')
!
!   -- recherche de la chaine nompar avec memoire sur tout 'calcul'
    ca_capoiz_ = ca_capoiz_ + 1
    if (ca_capoiz_ .gt. 512) then
        iparg = indik8(zk8(ca_iaoppa_),nompar,1,ca_npario_)
    else
        if (zk8(ca_iaoppa_-1+ca_caindz_(ca_capoiz_)) .eq. nompar) then
            iparg = ca_caindz_(ca_capoiz_)
        else
            iparg = indik8(zk8(ca_iaoppa_),nompar,1,ca_npario_)
            ca_caindz_(ca_capoiz_) = iparg
        endif
    endif
!
!
    if (iparg .eq. 0) then
        valk(1) = nompar
        valk(2) = ca_option_
        call utmess('E', 'CALCUL_15', nk=2, valk=valk)
        call contex_param(ca_option_, ' ')
    endif
!
!   -- on verifie que les parametre in sont en lecture
!      et que les parametres out sont en ecriture
    if (iparg .gt. ca_nparin_ .and. loue .eq. 'L') then
        write(6,*)'PARAMETRE OUT EN LECTURE : ',nompar
        ASSERT(.false.)
    else if (iparg.le.ca_nparin_ .and. loue.eq.'E') then
        write(6,*)'PARAMETRE IN EN ECRITURE : ',nompar
        ASSERT(.false.)
    endif
!
    iachlo=zi(ca_iawloc_-1+3*(iparg-1)+1)
    ilchlo=zi(ca_iawloc_-1+3*(iparg-1)+2)
    lgcata=zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2)
    debugr=zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)
!
    if (lgcata .eq. -1) then
        valk(1) = nompar
        valk(2) = ca_option_
        valk(3) = ca_nomte_
        call utmess('E', 'CALCUL_16', nk=3, valk=valk)
        call contex_param(ca_option_, nompar)
    endif
!
!
    if (iachlo .eq. -1) then
!        on ajoute cela pour emettre un message plus clair dans
!        le cas ou il manque un champ lie a un parametre
        call jenonu(jexnom('&CATA.OP.NOMOPT', ca_option_), opt)
        call jeveuo(jexnum('&CATA.OP.DESCOPT', opt), 'L', iaopd2)
        call jeveuo(jexnum('&CATA.OP.LOCALIS', opt), 'L', iaoplo)
        call jeveuo(jexnum('&CATA.OP.OPTPARA', opt), 'L', iapara)
        npari2 = zi(iaopd2-1+2)
        do ipara = 1, npari2
            if (zk8(iapara+ipara-1) .eq. nompar) goto 30
        end do
        goto 40
 30     continue
        valk(1) = ca_option_
!        on peut trouver d'ou vient le probleme dans 3 cas
        if (zk24(iaoplo+3*ipara-3) .eq. 'CARA') then
            call utmess('E', 'CALCUL_10', sk=valk(1))
        else if (zk24(iaoplo+3*ipara-3).eq.'CHMA') then
            call utmess('E', 'CALCUL_11', sk=valk(1))
        else if (zk24(iaoplo+3*ipara-3).eq.'MODL') then
            call utmess('E', 'CALCUL_12', sk=valk(1))
        endif
 40     continue
        valk(1) = nompar
        valk(2) = ca_option_
        valk(3) = ca_nomte_
        call utmess('E', 'CALCUL_17', nk=3, valk=valk)
        call contex_param(ca_option_, nompar)
!
    endif
    ASSERT(iachlo.ne.-2)
!
!
!   -- calcul de itab, lonchl, decael :
!   -----------------------------------
    call chloet(iparg, etendu, jceld)
    if (etendu) then
        adiel = zi(jceld-1+zi(jceld-1+4+ca_igr_)+4+4* (ca_iel_-1)+4)
        debgr2 = zi(jceld-1+zi(jceld-1+4+ca_igr_)+8)
        ASSERT(lgcata.eq.zi(jceld-1+zi(jceld-1+4+ca_igr_)+3))
        decael = (adiel-debgr2)
        lonchl = zi(jceld-1+zi(jceld-1+4+ca_igr_)+4+4* (ca_iel_-1)+3)
    else
        decael = (ca_iel_-1)*lgcata
        lonchl = lgcata
    endif
    itab = iachlo+debugr-1+decael
!
!   -- pour les champs "in" on verifie que l'extraction est
!      complete sur l'element:
!   ----------------------------------------------------------
    if (ilchlo .ne. -1) then
        do k = 1, lonchl
            if (.not.zl(ilchlo+debugr-1+decael-1+k)) then
                call tecael(iadzi, iazk24)
                nommai=zk24(iazk24-1+3)(1:8)
                valk(1) = nompar
                valk(2) = ca_option_
                valk(3) = ca_nomte_
                valk(4) = nommai
!
!               -- pour certains parametres "courants" on emet
!                  un message plus clair :
                if (nompar .eq. 'PMATERC') then
                    call utmess('F', 'CALCUL_20', nk=4, valk=valk)
                else if (nompar.eq.'PCACOQU') then
                    call utmess('F', 'CALCUL_21', nk=4, valk=valk)
                else if (nompar.eq.'PCAGNPO') then
                    call utmess('F', 'CALCUL_22', nk=4, valk=valk)
                else if (nompar.eq.'PCAORIE') then
                    call utmess('F', 'CALCUL_23', nk=4, valk=valk)
!
                else
!
                    write (6,*) 'ERREUR JEVECH ZL :',nompar, (zl(&
                    ilchlo+debugr-1+decael-1+kk),kk=1,lonchl)
                    write (6,*) 'MAILLE: ',zk24(iazk24-1+3)
                    call utmess('E', 'CALCUL_19', nk=4, valk=valk)
                    call contex_param(ca_option_, nompar)
                endif
            endif
        end do
    endif
!
end subroutine
