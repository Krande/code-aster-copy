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
subroutine ascavc(lchar , infcha   , fomult, numedd, inst, vci, dlci_, &
                  l_hho_, hhoField_, basez)
!
use HHO_type
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/ascova.h"
#include "asterfort/assert.h"
#include "asterfort/calvci.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
!
character(len=24) :: lchar, infcha, fomult
character(len=*) :: vci, numedd
real(kind=8) :: inst
character(len=*), optional :: dlci_
aster_logical, intent(in), optional :: l_hho_
type(HHO_Field), intent(in), optional :: hhoField_
character(len=1), intent(in), optional :: basez
! ----------------------------------------------------------------------
! BUT  :  CALCUL DU CHAM_NO CONTENANT LE VECTEUR LE CINEMATIQUE
! ---     ASSOCIE A LA LISTE DE CHAR_CINE_* LCHAR A UN INSTANT INST
!         AVEC LES FONCTIONS MULTIPLICATIVES FOMULT.
! ----------------------------------------------------------------------
! IN  K*24 LCHAR : NOM DE L'OJB S V K24 CONTENANT LA LISTE DES CHARGES
! IN  K*19 INFCHA : NOM DE L'OJB S V I CONTENANT LA LISTE DES INFO.
! IN  K*24 FOMULT : NOM DE L'OJB S V K24 CONTENANT LA LISTE DES FONC.
! IN  K*14 NUMEDD  : NOM DE LA NUMEROTATION SUPPORTANT LE CHAM_NO
! IN  R*8  INST   : VALE DU PARAMETRE INST.
! VAR/JXOUT  K*19 VCI    :  CHAM_NO RESULTAT
!   -------------------------------------------------------------------
!-----------------------------------------------------------------------
!----------------------------------------------------------------------
!     VARIABLES LOCALES
!----------------------------------------------------------------------
    integer :: idchar, jinfc, idfomu, nchtot, nchci, ichar, icine, ilchno
    integer :: ichci, ifm, niv, neq, ieq, jdlci2,  ieqmul
    character(len=1) :: base
    character(len=8) :: newnom
    character(len=19) :: charci, chamno, vci2
    character(len=24) :: vachci, dlci
    character(len=8) :: charge
    aster_logical :: l_hho
    integer, pointer :: v_dlci(:) => null()
    data chamno/'&&ASCAVC.???????'/
    data vachci/'&&ASCAVC.LISTE_CI'/
    data charci/'&&ASCAVC.LISTE_CHI'/
!----------------------------------------------------------------------
!
!
    call jemarq()
    if (vci .eq. ' ') vci='&&ASCAVC.VCI'
    vci2=vci
!
    call infniv(ifm, niv)
!
!
    if(present(basez)) then
        base = basez
    else
        base = 'V'
    end if
!
    if(present(dlci_)) then
        dlci = dlci_
    else
        dlci = '&&ASCAVC.DLCI'
    end if
!
    newnom='.0000000'
!
! - For HHO
!
    l_hho = ASTER_FALSE
    if (present(l_hho_)) then
        l_hho = l_hho_
    endif
!
    call jedetr(vachci)
    call jedetr(vci2//'.DLCI')
    call jeveuo(lchar, 'L', idchar)
    call jeveuo(infcha, 'L', jinfc)
    call jeveuo(fomult, 'L', idfomu)
!
    nchtot = zi(jinfc)
    nchci = 0
    ieqmul=0
!
    do ichar = 1, nchtot
        icine = zi(jinfc+ichar)
        if (icine .lt. 0) nchci=nchci+1
!       -- UNE CHARGE NON "CINEMATIQUE" PEUT EN CONTENIR UNE :
        charge=zk24(idchar-1+ichar)(1:8)
    end do
!
!
    call wkvect(vachci, 'V V K24', max(nchci, 1), ilchno)
!
!     -- S'IL N'Y A PAS DE CHARGES CINEMATIQUES, ON CREE UN CHAMP NUL:
    if (nchci .eq. 0) then
        call gcnco2(newnom)
        chamno(10:16) = newnom(2:8)
        call corich('E', chamno, ichin_ = -2)
        call vtcreb(chamno, 'V', 'R',&
                    nume_ddlz = numedd,&
                    nb_equa_outz = neq)
        zk24(ilchno-1+1) = chamno
!
!
!     -- S'IL Y A DES CHARGES CINEMATIQUES :
    else
!
        ichci = 0
        call dismoi('NB_EQUA', numedd, 'NUME_DDL', repi=neq)
        call wkvect(dlci, base//' V I', neq, jdlci2)
        do ichar = 1, nchtot
            charge=zk24(idchar-1+ichar)(1:8)
            icine = zi(jinfc+ichar)
            if (icine .lt. 0) then
                ichci = ichci + 1
                call gcnco2(newnom)
                chamno(10:16) = newnom(2:8)
                call corich('E', chamno, ichin_ = ichar)
                zk24(ilchno-1+ichci) = chamno
                if (l_hho) then
                    call calvci(chamno, numedd, 1, charge, inst, &
                                'V',l_hho, hhoField_)
                else
                    call calvci(chamno, numedd, 1, charge, inst, &
                                'V',l_hho)
                endif
                call jeveuo(chamno//'.DLCI', 'L', vi=v_dlci)
!           --- COMBINAISON DES DLCI (OBJET CONTENANT DES 0 OU DES 1),
!           --- LES 1 ETANT POUR LES DDL CONTRAINT
!           --- LE RESTE DE L OBJECT VCI2 EST CREE PAR ASCOVA
                do ieq = 1, neq
!             -- ON REGARDE SI UN DDL N'EST PAS ELIMINE PLUSIEURS FOIS:
                    if (zi(jdlci2-1+ieq) .gt. 0 .and. v_dlci(ieq) .gt. 0) ieqmul=ieq
!
                    zi(jdlci2-1+ieq)=max(zi(jdlci2-1+ieq),v_dlci(&
                    ieq))
                end do
                call jedetr(chamno//'.DLCI')
            endif
        end do
    endif

!
!     -- ON COMBINE LES CHAMPS CALCULES :
    call ascova('D', vachci, fomult, 'INST', inst,&
                'R', vci2, base)
!
!     --SI ON A PAS DE CHARGE CINEMATIQUE, IL FAUT QUAND MEME
!        FAIRE LE MENAGE
    if (nchci .eq. 0) call detrsd('CHAMP_GD', chamno(1:19))
    call jedetr(charci)
!
    if(.not.present(dlci_)) then
        call jedetr(dlci)
    end if
!
    call jedema()
end subroutine
