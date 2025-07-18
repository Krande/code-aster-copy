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
subroutine giecma(nfic, trouve, nbele, nomobj, tymail, &
                  nbno, ecrma, icoma)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/codent.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8) :: nfic, nbele, nbno, icoma, ibid
    character(len=8) :: tymail, nomobj
    aster_logical :: trouve, ecrma(*)
! ----------------------------------------------------------------------
!     BUT: ECRIRE SUR LE FICHIER DE MAILLAGE ASTER
!          LES MAILLES CORRESPONDANT A L'OBJET GIBI
!     ( SI CES MAILLES SONT EN DOUBLE, ON NE LES REECRIT PAS)
!
!     IN : NFIC : UNITE D'ECRITURE
!          TROUVE : LE GROUP_MA COURANT EST-IL A TRAITE
!            SI OUI : ON ECRIT LA CONNECTIVITE ET ON INCREMENTE ICOMA
!            SI NON : ON INCREMENTE ICOMA MAIS ON N'ECRIT RIEN.
!          NBELE: NOMBRE DE MAILLES DANS L'OBJET.
!          NOMBOJ:NOM DE LA SD CONTENANT LA CONNECTIVITE.
!          TYMAIL:TYPE_MAILLE (GIBI)
!          NBNO : NOMBRE DE NOEUD DE TYMAIL.
!     VAR: ICOMA : COMPTEUR DE MAILLE.
!          (ICOMA EST INCREMENTE DE NBELE A CHAQUE APPEL)
!
! ----------------------------------------------------------------------
!
!     VARIABLES LOCALES:
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iacorr, ibvec, icoj
    integer(kind=8) :: icok, ii, itymai, ivect, j, k, l
    integer(kind=8) :: maili, maille, nbelem, nbfois, nbrest, nmtot, numno
!
!-----------------------------------------------------------------------
    parameter(nbelem=18)
    character(len=7) :: k7nom(8)
    character(len=8) :: k8nom(8), tymagi(nbelem), tymaas(nbelem)
!
!
!     COGIAS EST UN TAMPON POUR ECRIRE LA CONNECTIVITE DES MAILLES
!        DANS L'ORDRE ASTER . (27 EST LE MAX DE NOEUDS POSSIBLE).
    integer(kind=8) :: cogias(27)
    integer(kind=8), pointer :: indirect(:) => null()
    integer(kind=8), pointer :: connex(:) => null()
    integer(kind=8), pointer :: numanew(:) => null()
    data tymaas/'POI1    ', 'SEG2    ', 'SEG3    ', 'TRIA3   ',&
     &     'TRIA6   ', 'QUAD4   ', 'QUAD8   ', 'QUAD9   ', 'TETRA4  ',&
     &     'TETRA10 ', 'PENTA6  ', 'PENTA15 ', 'HEXA8   ', 'HEXA20  ',&
     &     'HEXA27  ', 'PYRAM5  ', 'PYRAM13 ', '????    '/
    data tymagi/'POI1    ', 'SEG2    ', 'SEG3    ', 'TRI3    ',&
     &     'TRI6    ', 'QUA4    ', 'QUA8    ', 'QUA9    ', 'TET4    ',&
     &     'TE10    ', 'PRI6    ', 'PR15    ', 'CUB8    ', 'CU20    ',&
     &     'CU27    ', 'PYR5    ', 'PY13    ', '????    '/
!
!
    call jemarq()
    if (nbno .gt. 27) then
        call utmess('F', 'PREPOST_54')
    end if
!
    call jeveuo('&&GILIRE'//nomobj//'.CONNEX', 'L', vi=connex)
    call jeveuo('&&GILIRE.NUMANEW', 'E', vi=numanew)
    call jelira('&&GILIRE.NUMANEW', 'LONUTI', nmtot)
    call jeexin('&&GILIRE.INDIRECT', ibid)
    if (ibid .eq. 0) then
        call utmess('F', 'PREPOST_55')
    end if
    call jeveuo('&&GILIRE.INDIRECT', 'L', vi=indirect)
!
    call jeexin('&&GILIRE.VECT', ibvec)
    if (ibvec .eq. 0) then
        call wkvect('&&GILIRE.VECT', 'V V I', nmtot, ivect)
        do i = 1, nmtot
            zi(ivect+i-1) = 0
        end do
    else
        call jeveuo('&&GILIRE.VECT', 'L', ivect)
    end if
!
!
!  -- ON VERIFIE QUE LE GROUPE COURANT EST NOMME OU SOUS GROUPE
!
    if (.not. trouve) then
        icoma = icoma+nbele
        goto 999
    end if
!
    call jeveuo(jexnom('&&GILIRE.CORR_GIBI_ASTER', tymail), 'L', iacorr)
    itymai = indik8(tymagi(1), tymail, 1, nbelem)
    if (itymai .eq. 0) then
        call utmess('F', 'PREPOST_56', sk=tymail)
    end if
!
    write (nfic, *) tymaas(itymai)
!
!     -- BOUCLE SUR LES MAILLES DE L'OBJET SIMPLE:
!     --------------------------------------------
    do i = 1, nbele
!
!
        icoma = icoma+1
        maille = numanew(icoma)
!
!        -- SI LA MAILLE N'A PAS SON NUMERO INITIAL
!           ET SI ELLE EST DEJA ECRITE ON SORT
        if ((maille .ne. icoma) .and. (ecrma(maille))) goto 1
!
! SI LA MAILLE N A PAS LE NUMERO COURANT ET QU'ELLE
! N' A PAS ETE ECRITE ON ECRIT LE NOEUD COURANT
!
        if ((maille .ne. icoma) .and. (.not. (ecrma(maille)))) then
            if (zi(ivect+maille-1) .eq. 0) zi(ivect+maille-1) = icoma
            do ii = 1, nmtot
                maili = numanew(ii)
                if (maili .eq. maille) then
                    numanew(ii) = zi(ivect+maille-1)
                    ecrma(numanew(ii)) = .true.
                end if
            end do
        end if
!
!
        ecrma(maille) = .true.
!
        call codent(icoma, 'G', k7nom(1))
        k8nom(1) = 'M'//k7nom(1)
!
!        -- REMPLISSAGE DE COGIAS:
        do j = 1, nbno
            numno = connex(nbno*(i-1)+j)
            cogias(j) = indirect(numno)
        end do
!
        nbfois = nbno/7
        nbrest = nbno-7*nbfois
        icoj = 0
        icok = 0
!
        do j = 1, nbfois
            do k = 1, 7
                icok = icok+1
                numno = cogias(zi(iacorr-1+icok))
                call codent(numno, 'G', k7nom(1+k))
                k8nom(1+k) = 'N'//k7nom(1+k)
            end do
            write (nfic, 1001) (k8nom(l), l=1, 8)
            k8nom(1) = ' '
            icoj = icoj+7
        end do
!
        do k = 1, nbrest
            icok = icok+1
            numno = cogias(zi(iacorr-1+icok))
            call codent(numno, 'G', k7nom(1+k))
            k8nom(1+k) = 'N'//k7nom(1+k)
        end do
        write (nfic, 1001) (k8nom(l), l=1, nbrest+1)
1       continue
    end do
    write (nfic, *) 'FINSF'
    write (nfic, *) '%'
!
999 continue
    call jelibe('&&GILIRE.VECT')
!
1001 format(2x, a8, 7(1x, a8), 1x)
!
    call jedema()
end subroutine
