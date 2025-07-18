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
subroutine op0147()
    implicit none
!   CALCUL DES INTERSPECTRES DE REPONSE MODALE (DYNA_SPEC_MODAL)
!      LE CONCEPT PRODUIT EST UN INTERSPECTRE
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/calcsp.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/ordis.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ifreq, im, imasg, imod1, inumo
    integer(kind=8) :: inuor, ivite, jnuor, nbm, nbmr
    integer(kind=8) :: nnn, npv, i1, i3, ivitef
!-----------------------------------------------------------------------
    aster_logical :: casint
    character(len=8) :: table, nomu, option
    character(len=16) :: concep, cmd
    character(len=19) :: base
    character(len=24) :: freq, masg, vite, numo, nomobj, chnumi
    integer(kind=8) :: lnumi, lrefes
    real(kind=8) :: epsi, val, vitef
    character(len=16), pointer :: refe(:) => null()
!
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
!
    call getres(nomu, concep, cmd)
!
! --- 1.RECUPERATION DES OBJETS DE LA BASE MODALE PERTURBEE ---
!
    call getvid(' ', 'BASE_ELAS_FLUI', scal=base)
!
    masg = base//'.MASG'
    vite = base//'.VITE'
    freq = base//'.FREQ'
    numo = base//'.NUMO'
!
    call jeveuo(masg, 'L', imasg)
!
    call jeveuo(vite, 'L', ivite)
    call jelira(vite, 'LONUTI', npv)
    call getvr8(' ', 'VITE_FLUI', scal=vitef)
    call getvr8(' ', 'PRECISION', scal=epsi)
!
    ivitef = 0
    do i3 = 1, npv
        val = zr(ivite-1+i3)-vitef
        if (abs(val) .lt. epsi) then
            ivitef = i3
        end if
    end do
    if (ivitef .eq. 0) then
        call utmess('F', 'ALGELINE3_25', sr=vitef)
    end if
!
    call jeveuo(freq, 'L', ifreq)
    call jelira(freq, 'LONUTI', nbm)
    nbm = nbm/(2*npv)
!
    call jeveuo(numo, 'L', inumo)
!
! --- 2.RECUPERATION DU NOM DE LA TABLE ---
!
    call getvid('EXCIT ', 'INTE_SPEC_GENE', iocc=1, scal=table)
!
!     VERIFICATION DES PARAMETRES DE LA TABLE
!
    chnumi = table//'.NUMI'
    call jeveuo(chnumi, 'L', lnumi)
    call jelira(chnumi, 'LONMAX', nbmr)
!
    nomobj = '&&OP0147.TEMP.NUOR'
    call wkvect(nomobj, 'V V I', nbmr, jnuor)
    do i1 = 1, nbmr
        zi(jnuor-1+i1) = zi(lnumi-1+i1)
    end do
    call ordis(zi(jnuor), nbmr)
    call wkvect('&&OP0147.MODE', 'V V I', nbmr, inuor)
    nnn = 1
    zi(inuor) = zi(jnuor)
    do i = 2, nbmr
        if (zi(jnuor+i-1) .eq. zi(inuor+nnn-1)) goto 20
        nnn = nnn+1
        zi(inuor+nnn-1) = zi(jnuor+i-1)
20      continue
    end do
    nbmr = nnn
    do im = 1, nbm
        if (zi(inumo+im-1) .eq. zi(inuor)) then
            imod1 = im
            goto 31
        end if
    end do
    call utmess('F', 'MODELISA5_78')
31  continue
!
! --- 3.RECUPERATION DE L'OPTION DE CALCUL ---
!
    casint = .true.
    call getvtx(' ', 'OPTION', scal=option)
    if (option(1:4) .eq. 'DIAG') casint = .false.
!
    call jeveuo(table//'.REFE', 'L', vk16=refe)
    if (refe(2) (1:4) .eq. 'DIAG' .and. casint) then
        call utmess('F', 'MODELISA5_79')
    end if
!
! --- 4.CREATION DE LA STRUCTURE RESULTAT ET CALCUL DE LA REPONSE ---
! ---   PAR CALCSP                                                ---
!
    call wkvect(nomu//'.REFE', 'G V K16', 3, lrefes)
    zk16(lrefes) = 'DEPL_GENE'
    zk16(lrefes+1) = option
    zk16(lrefes+2) = 'FREQ'
!
    call calcsp(casint, nomu, table, zr(ifreq), zr(imasg), &
                nbm, nbmr, imod1, zi(inuor), ivitef)
!
    call titre()
!
!
    call jedema()
end subroutine
