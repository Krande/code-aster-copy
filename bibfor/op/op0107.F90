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

subroutine op0107()
    implicit none
!     OPERATEUR   POST_ELEM
!     ------------------------------------------------------------------
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/chpve2.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medomp.h"
#include "asterfort/peaire.h"
#include "asterfort/pecage.h"
#include "asterfort/pecapo.h"
#include "asterfort/pechli.h"
#include "asterfort/peecin.h"
#include "asterfort/peeint.h"
#include "asterfort/peepot.h"
#include "asterfort/peingl.h"
#include "asterfort/pemain.h"
#include "asterfort/pemima.h"
#include "asterfort/penorm.h"
#include "asterfort/peritr.h"
#include "asterfort/pevolu.h"
#include "asterfort/peweib.h"
#include "asterfort/pewext.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    integer(kind=8) :: nh, iret, jordr, n1, n2, nbocc, nbordr, nc, np, nr, ier
    real(kind=8) :: prec
    character(len=8) :: k8b, modele, carele, deform, resuco, crit, mesh
    character(len=16) :: concep, nomcmd
    character(len=19) :: resu, knum, tabtyp(3)
    character(len=24) :: mate, mateco, chdef
!
!     ------------------------------------------------------------------
!
    call jemarq()
!
    call getres(resu, concep, nomcmd)
    call getvid(' ', 'RESULTAT', scal=resuco, nbret=nr)
!
    if (nr .eq. 0) resuco = ' '
!
    call infmaj()
!
    call getfac('TRAV_EXT', nbocc)
    if (nbocc .ne. 0) then
        call pewext(resu)
    end if
!
    call getfac('CHAR_LIMITE', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mateco=mateco)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pechli(resu, modele, mateco)
    end if
!
    call getfac('AIRE_INTERNE', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peaire(resu, mesh, nbocc)
    end if
!
    call getfac('MASS_INER', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        chdef = ' '
        call getvtx(' ', 'GEOMETRIE', scal=deform, nbret=n1)
        if (deform .eq. 'DEFORMEE') then
            call getvid(' ', 'CHAM_GD', scal=chdef, nbret=n2)
            if (n2 .eq. 0) then
                tabtyp(1) = 'NOEU#DEPL_R'
                tabtyp(2) = 'NOEU#TEMP_R'
                tabtyp(3) = 'ELEM#ENER_R'
                knum = '&&OP0107.NUME_ORDRE'
                call getvid(' ', 'RESULTAT', scal=resuco, nbret=nr)
                call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
                call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
                call rsutnu(resuco, ' ', 0, knum, nbordr, &
                            prec, crit, iret)
                if (nbordr .ne. 1) then
                    call utmess('F', 'POSTELEM_10')
                end if
                if (iret .ne. 0) goto 999
                call jeveuo(knum, 'L', jordr)
                call rsexch('F', resuco, 'DEPL', zi(jordr), chdef, &
                            iret)
                call chpve2(chdef, 3, tabtyp, ier)
            end if
        end if
        call pemain(resu, modele, mate, mateco, carele, nh, &
                    nbocc, chdef)
!
    end if
!
    call getfac('ENER_POT', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peepot(resu, modele, mate, mateco, carele, nh, &
                    nbocc)
!
    end if
!
    call getfac('ENER_CIN', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peecin(resu, modele, mate, mateco, carele, nh, &
                    nbocc)
!
    end if
!
    call getfac('INTEGRALE', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele)
        call peeint(resu, modele, nbocc)
    end if
!
    call getfac('NORME', nbocc)
    if (nbocc .ne. 0) then
!         --- ON RECUPERE LE MODELE
        call getvid('NORME', 'CHAM_GD', iocc=1, scal=chdef, nbret=n1)
        if (n1 .ne. 0) then
            call getvid('NORME', 'MODELE', iocc=1, scal=modele, nbret=n2)
        else
            call getvid('NORME', 'RESULTAT', iocc=1, scal=resuco, nbret=nr)
            call medomp(resuco, modele)
        end if
        call penorm(resu, modele)
    end if
!
    call getfac('VOLUMOGRAMME', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, carele=carele)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pevolu(resu, modele, carele, nbocc)
    end if
!
    call getfac('MINMAX', nbocc)
    if (nbocc .ne. 0) then
        call getvid('MINMAX', 'CHAM_GD', iocc=1, scal=chdef, nbret=n1)
        if (n1 .ne. 0) then
            call getvid('MINMAX', 'MODELE', iocc=1, scal=modele, nbret=n2)
        else
            call getvid('MINMAX', 'RESULTAT', iocc=1, scal=resuco, nbret=nr)
            call medomp(resuco, modele)
        end if
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ! ASSERT(.not. isParallelMesh(mesh))
        call pemima(n1, chdef, resu, modele, nbocc)
    end if
!
    call getfac('WEIBULL', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peweib(resu, modele, mate, mateco, carele, k8b, &
                    nh, nbocc, 0, nomcmd)
    end if
!
    call getfac('RICE_TRACEY', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, carele=carele, nh=nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peritr(resu, modele, carele, nh, nbocc)
    end if
!
    call getfac('CARA_GEOM', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pecage(resu, modele, nbocc)
    end if
!
    call getfac('CARA_POUTRE', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, carele=carele, nh=nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pecapo(resu, modele, carele, nh)
    end if
!
    call getfac('INDIC_ENER', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(resu, modele, mate, mateco, carele, nh, &
                    nbocc, 'INDIC_ENER')
    end if
!
    call getfac('INDIC_SEUIL', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(resu, modele, mate, mateco, carele, nh, &
                    nbocc, 'INDIC_SEUIL')
    end if
!
    call getfac('ENER_ELAS', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(resu, modele, mate, mateco, carele, nh, &
                    nbocc, 'ENER_ELAS')
    end if
!
    call getfac('ENER_ELTR', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(resu, modele, mate, mateco, carele, nh, &
                    nbocc, 'ENER_ELTR')
    end if

!
    call getfac('ENER_TOTALE', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(resu, modele, mate, mateco, carele, nh, &
                    nbocc, 'ENER_TOTALE')
    end if
!
    call getfac('ENER_DISS', nbocc)
    if (nbocc .ne. 0) then
        call medomp(resuco, modele, mate, mateco, carele, nh)
        call dismoi('NOM_MAILLA', modele, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(resu, modele, mate, mateco, carele, nh, &
                    nbocc, 'ENER_DISS')
    end if
!
999 continue
    call titre()
!
    call jedema()
!
end subroutine
