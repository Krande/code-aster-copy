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

subroutine op0130()
    implicit none
!
!     OPERATEUR "POST_DYNA_MODA_T"
!
! ----------------------------------------------------------------------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedetr.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pochoc.h"
#include "asterfort/pochpv.h"
#include "asterfort/porefd.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/titre.h"
    integer(kind=8) :: nbbloc, nbclas, i, nbind, nbno, jnomno
    real(kind=8) :: tdebut, tfin, offset, trepos
    character(len=8) :: base, modele, maillage
    character(len=8) :: trange, noeu, cmp, nomres
    character(len=16) :: nomcmd, concep, koptio, motcle(2), typmcl(2)
    character(len=24) :: nomno
    aster_logical :: loptio
    integer(kind=8), pointer :: desc(:) => null()
!   ------------------------------------------------------------------
    data motcle/'NOEUD', 'GROUP_NO'/
    data typmcl/'NOEUD', 'GROUP_NO'/
!   ------------------------------------------------------------------

!
    call jemarq()
    nomres = ' '
    call getres(nomres, concep, nomcmd)
    call infmaj()
!
    call getvid(' ', 'RESU_GENE', scal=trange)
    call jeveuo(trange//'           .DESC', 'L', vi=desc)
!
    call getfac('CHOC', nbind)
    if (nbind .ne. 0) then
        do i = 1, nbind
            call getvis('CHOC', 'NB_BLOC', iocc=i, scal=nbbloc)
            call getvr8('CHOC', 'INST_INIT', iocc=i, scal=tdebut)
            call getvr8('CHOC', 'INST_FIN', iocc=i, scal=tfin)
            call getvr8('CHOC', 'SEUIL_FORCE', iocc=i, scal=offset)
            call getvr8('CHOC', 'DUREE_REPOS', iocc=i, scal=trepos)
            call getvtx('CHOC', 'OPTION', iocc=i, scal=koptio)
            call getvis('CHOC', 'NB_CLASSE', iocc=i, scal=nbclas)
            if (koptio(1:6) .eq. 'USURE') then
                loptio = .true.
            else
                loptio = .false.
            end if
!           --- Constant integration step
            if (desc(1) .eq. 2) then
                call pochoc(trange, nbbloc, tdebut, tfin, offset, &
                            trepos, nbclas, nomres, loptio)
!           --- Non constant integration
            else if (desc(1) .eq. 3) then
                call pochpv(trange, nbbloc, tdebut, tfin, offset, &
                            trepos, nbclas, nomres, loptio)
            end if
        end do
    end if
!
    call getfac('RELA_EFFO_DEPL', nbind)
    nomno = '&&OP0130.MES_NOEUDS'
    if (nbind .ne. 0 .and. desc(4) .ne. 0) then
        do i = 1, nbind
            call dismoi('BASE_MODALE', trange, 'RESU_DYNA', repk=base, arret='F')
            call dismoi('NOM_MODELE', base, 'RESULTAT', repk=modele)
            call dismoi('NOM_MAILLA', base, 'RESULTAT', repk=maillage)
!
            call reliem(modele, maillage, 'NO_NOEUD', 'RELA_EFFO_DEPL', i, &
                        2, motcle, typmcl, nomno, nbno)
            call jeveuo(nomno, 'L', jnomno)
            if (nbno .ne. 1) then
                call utmess('F', 'MODELISA5_67', sk='RELA_EFFO_DEPL', si=nbno)
            end if
            noeu = zk8(jnomno)
            call getvtx('RELA_EFFO_DEPL', 'NOM_CMP', iocc=i, scal=cmp)
!
            call porefd(trange, noeu, cmp, nomres)
        end do
    end if
    call jedetr(nomno)
!
    call titre()
!
    call jedema()
end subroutine
