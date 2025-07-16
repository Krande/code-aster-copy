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
subroutine jevech(nmparz, louez, itab, vl, vi, &
                  vi4, vr, vc, vk8, vk16, &
                  vk24, vk32, vk80)
!
    use calcul_module, only: ca_caindz_, ca_capoiz_, ca_iaoppa_, ca_iawlo2_, ca_iawloc_, ca_iel_, &
                             ca_igr_, ca_nbgr_, ca_nomte_, ca_nparin_, ca_npario_, ca_option_
    use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
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
#include "asterfort/jgetptc.h"
!
    character(len=*), intent(in) :: nmparz, louez
    integer(kind=8), optional :: itab
    aster_logical, pointer, optional :: vl(:)
    integer(kind=8), pointer, optional :: vi(:)
    integer(kind=4), pointer, optional :: vi4(:)
    real(kind=8), pointer, optional :: vr(:)
    complex(kind=8), pointer, optional :: vc(:)
    character(len=8), pointer, optional :: vk8(:)
    character(len=16), pointer, optional :: vk16(:)
    character(len=24), pointer, optional :: vk24(:)
    character(len=32), pointer, optional :: vk32(:)
    character(len=80), pointer, optional :: vk80(:)
!-----------------------------------------------------------------
!  entrees:
!     nompar  : nom du parametre de l'option
!     louez   : 'L' ou 'E'  ( lecture/ecriture )
!
!  sorties:
!     itab     : adresse du champ local correspondant a nompar
!-----------------------------------------------------------------
    integer(kind=8) :: iachlo, jtab
    integer(kind=8) :: ilchlo, k, kk, debugr
    integer(kind=8) :: iparg, lgcata
    integer(kind=8) :: jceld, adiel, numCell
    integer(kind=8) :: debgr2, lonchl, decael, iadzi, iazk24
    integer(kind=8) :: opt, iaopd2, iaoplo, iapara, ipara, npari2
    aster_logical :: etendu
    character(len=8) :: nompar
    character(len=1) :: loue
    character(len=24) :: valk(5)
    type(c_ptr) :: pc
! ------------------------------------------------------------------
    nompar = nmparz
    loue = louez
!
    ASSERT(loue .eq. 'L' .or. loue .eq. 'E')
!
!   -- recherche de la chaine nompar avec memoire sur tout 'calcul'
    ca_capoiz_ = ca_capoiz_+1
    if (ca_capoiz_ .gt. 512) then
        iparg = indik8(zk8(ca_iaoppa_), nompar, 1, ca_npario_)
    else
        if (zk8(ca_iaoppa_-1+ca_caindz_(ca_capoiz_)) .eq. nompar) then
            iparg = ca_caindz_(ca_capoiz_)
        else
            iparg = indik8(zk8(ca_iaoppa_), nompar, 1, ca_npario_)
            ca_caindz_(ca_capoiz_) = iparg
        end if
    end if
!
!
    if (iparg .eq. 0) then
        valk(1) = nompar
        valk(2) = ca_option_
        call utmess('E', 'CALCUL_15', nk=2, valk=valk)
        call contex_param(ca_option_, ' ')
    end if
!
!   -- on verifie que les parametre in sont en lecture
!      et que les parametres out sont en ecriture
    if (iparg .gt. ca_nparin_ .and. loue .eq. 'L') then
        write (6, *) 'PARAMETRE OUT EN LECTURE : ', nompar
        ASSERT(.false.)
    else if (iparg .le. ca_nparin_ .and. loue .eq. 'E') then
        write (6, *) 'PARAMETRE IN EN ECRITURE : ', nompar
        ASSERT(.false.)
    end if
!
    iachlo = zi(ca_iawloc_-1+3*(iparg-1)+1)
    ilchlo = zi(ca_iawloc_-1+3*(iparg-1)+2)
    lgcata = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+2)
    debugr = zi(ca_iawlo2_-1+5*(ca_nbgr_*(iparg-1)+ca_igr_-1)+5)
!
    if (lgcata .eq. -1) then
        valk(1) = nompar
        valk(2) = ca_option_
        valk(3) = ca_nomte_
        call utmess('E', 'CALCUL_16', nk=3, valk=valk)
        call contex_param(ca_option_, nompar)
    end if
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
30      continue
        valk(1) = ca_option_
!        on peut trouver d'ou vient le probleme dans 3 cas
        if (zk24(iaoplo+3*ipara-3) .eq. 'CARA') then
            call utmess('E', 'CALCUL_10', sk=valk(1))
        else if (zk24(iaoplo+3*ipara-3) .eq. 'CHMA') then
            call utmess('E', 'CALCUL_11', sk=valk(1))
        else if (zk24(iaoplo+3*ipara-3) .eq. 'MODL') then
            call utmess('E', 'CALCUL_12', sk=valk(1))
        end if
40      continue
        valk(1) = nompar
        valk(2) = ca_option_
        valk(3) = ca_nomte_
        call utmess('E', 'CALCUL_17', nk=3, valk=valk)
        call contex_param(ca_option_, nompar)
!
    end if
    ASSERT(iachlo .ne. -2)
!
!
!   -- calcul de jtab, lonchl, decael :
!   -----------------------------------
    call chloet(iparg, etendu, jceld)
    if (etendu) then
        adiel = zi(jceld-1+zi(jceld-1+4+ca_igr_)+4+4*(ca_iel_-1)+4)
        debgr2 = zi(jceld-1+zi(jceld-1+4+ca_igr_)+8)
        ASSERT(lgcata .eq. zi(jceld-1+zi(jceld-1+4+ca_igr_)+3))
        decael = (adiel-debgr2)
        lonchl = zi(jceld-1+zi(jceld-1+4+ca_igr_)+4+4*(ca_iel_-1)+3)
    else
        decael = (ca_iel_-1)*lgcata
        lonchl = lgcata
    end if
    jtab = iachlo+debugr-1+decael
!
!   -- pour les champs "in" on verifie que l'extraction est
!      complete sur l'element:
!   ----------------------------------------------------------
    if (ilchlo .ne. -1) then
        do k = 1, lonchl
            if (.not. zl(ilchlo+debugr-1+decael-1+k)) then
                call tecael(iadzi, iazk24)
                valk(1) = nompar
                valk(2) = ca_option_
                valk(3) = ca_nomte_
                numCell = zi(iadzi)
!
!               -- pour certains parametres "courants" on emet
!                  un message plus clair :
                if (nompar .eq. 'PMATERC') then
                    call utmess('F', 'CALCUL_20', nk=3, valk=valk, si=numCell)
                else if (nompar .eq. 'PCACOQU') then
                    call utmess('F', 'CALCUL_21', nk=3, valk=valk, si=numCell)
                else if (nompar .eq. 'PCAGNPO') then
                    call utmess('F', 'CALCUL_22', nk=3, valk=valk, si=numCell)
                else if (nompar .eq. 'PCAORIE') then
                    call utmess('F', 'CALCUL_23', nk=3, valk=valk, si=numCell)
!
                else
!
                    write (6, *) 'ERREUR JEVECH ZL :', nompar, (zl(ilchlo+debugr-1+decael-1+kk), &
                                                                kk=1, lonchl)
                    write (6, *) 'MAILLE: ', zk24(iazk24-1+3)
                    call utmess('E', 'CALCUL_19', nk=3, valk=valk, si=numCell)
                    call contex_param(ca_option_, nompar)
                end if
            end if
        end do
    end if
!
    if (present(itab)) then
        itab = jtab
    end if
    if (present(vl)) then
        call jgetptc(jtab, pc, vl=zl(1))
        call c_f_pointer(pc, vl, [lonchl])
!
    else if (present(vi)) then
        call jgetptc(jtab, pc, vi=zi(1))
        call c_f_pointer(pc, vi, [lonchl])
!
    else if (present(vi4)) then
        call jgetptc(jtab, pc, vi4=zi4(1))
        call c_f_pointer(pc, vi4, [lonchl])
!
    else if (present(vr)) then
        call jgetptc(jtab, pc, vr=zr(1))
        call c_f_pointer(pc, vr, [lonchl])
!
    else if (present(vc)) then
        call jgetptc(jtab, pc, vc=zc(1))
        call c_f_pointer(pc, vc, [lonchl])
!
    else if (present(vk8)) then
        call jgetptc(jtab, pc, vk8=zk8(1))
        call c_f_pointer(pc, vk8, [lonchl])
!
    else if (present(vk16)) then
        call jgetptc(jtab, pc, vk16=zk16(1))
        call c_f_pointer(pc, vk16, [lonchl])
!
    else if (present(vk24)) then
        call jgetptc(jtab, pc, vk24=zk24(1))
        call c_f_pointer(pc, vk24, [lonchl])
!
    else if (present(vk32)) then
        call jgetptc(jtab, pc, vk32=zk32(1))
        call c_f_pointer(pc, vk32, [lonchl])
!
    else if (present(vk80)) then
        call jgetptc(jtab, pc, vk80=zk80(1))
        call c_f_pointer(pc, vk80, [lonchl])
!
    else
        ASSERT(present(itab))
    end if
!
end subroutine
