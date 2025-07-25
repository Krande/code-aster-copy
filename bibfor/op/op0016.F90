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
subroutine op0016()
! person_in_charge: j-pierre.lefebvre at edf.fr
    implicit none
!
!     DIRECTIVE IMPR_JEVEUX
!
#include "jeveux.h"
#include "asterfort/getvis.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/iunifi.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeimpa.h"
#include "asterfort/jeimpd.h"
#include "asterfort/jeimpm.h"
#include "asterfort/jeimpo.h"
#include "asterfort/jeimpr.h"
#include "asterfort/jelira.h"
#include "asterfort/jeprat.h"
#include "asterfort/jepreg.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/jjvern.h"
#include "asterfort/uldefi.h"
!
    character(len=8) :: fich
    character(len=24) :: nomobj, nom
    character(len=32) :: noml32
    character(len=80) :: txt
    character(len=16) :: noment
    character(len=1) :: nomcla
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, imes, info, ires, iret, iuni
    integer(kind=8) :: n, n1, n2, n3, nfic, nif, noc
    integer(kind=8) :: nrg, num, numerg, nuni
!-----------------------------------------------------------------------
    call infmaj()
    ires = iunifi('RESULTAT')
    imes = iunifi('MESSAGE')
    fich = 'MESSAGE'
    iuni = imes
!
    call getvis('IMPRESSION', 'UNITE', iocc=1, scal=iuni, nbret=nuni)
    call uldefi(iuni, ' ', ' ', 'A', 'N', &
                'O')
    call getvtx('IMPRESSION', 'NOM', iocc=1, scal=fich, nbret=nfic)
    if (nfic .ne. 0) then
!
        if (fich(1:8) .eq. 'RESULTAT') then
            iuni = ires
        else if (fich(1:7) .eq. 'MESSAGE') then
            iuni = imes
        else
            write (imes, *) fich, ': FICHIER INCONNU'
            write (imes, *) 'LES ECRITURES SONT SUR LE FICHIER MESSAGE'
        end if
        if (nuni .ne. 0) then
            write (imes, *) 'CHOISIR ENTRE LES MOTS-CLE NOM ET UNITE'
            write (imes, *) 'LES ECRITURES SONT SUR LE FICHIER MESSAGE'
        end if
    end if
!
!
    call getvtx(' ', 'ENTITE', scal=noment, nbret=n)
    if (n .eq. 0) goto 999
!
    call getvtx(' ', 'COMMENTAIRE', scal=txt, nbret=n)
    if (n .eq. 0) txt = ' '
!
    if (noment(1:6) .eq. 'DISQUE') then
!
        nomcla = ' '
        call getvtx(' ', 'CLASSE', scal=nomcla, nbret=n)
        write (iuni, *) ' '
        call jeimpd(iuni, nomcla, txt)
        write (iuni, *) ' '
!
    else if (noment(1:7) .eq. 'MEMOIRE') then
!
        write (iuni, *) ' '
        call jeimpm(iuni)
        write (iuni, *) ' '
!
    else if (noment .eq. 'REPERTOIRE') then
!
        nomcla = ' '
        call getvtx(' ', 'CLASSE', scal=nomcla, nbret=n)
        write (iuni, *) ' '
        call jeimpr(iuni, nomcla, txt)
        write (iuni, *) ' '
!
    else if (noment(1:5) .eq. 'OBJET') then
!
        call getvtx(' ', 'NOMOBJ', scal=nomobj, nbret=n)
        noml32 = nomobj
        call jjvern(noml32, 0, iret)
        write (iuni, *) ' '
        write (iuni, *) ' '
        if (iret .eq. 0) then
            write (iuni, *) ' DIRECTIVE IMPR_JEVEUX '
            write (iuni, *) ' L''OBJET "', nomobj, '" N''EXISTE PAS'
            goto 999
        else
            write (iuni, *) ' '
            write (iuni, *) ' ECRITURE DE L''OBJET : "', nomobj, '"'
        end if
        call getvtx(' ', 'COMMENTAIRE', scal=txt, nbret=n)
        if (n .eq. 0) txt = ' '
        write (iuni, *) ' '
        call jeimpa(iuni, nomobj, txt)
        write (iuni, *) ' '
        write (iuni, *) ' '
        if (iret .eq. 2) then
            call getvis(' ', 'NUMOC', scal=num, nbret=n1)
            call getvtx(' ', 'NOMOC', scal=nom, nbret=n2)
            call getvtx(' ', 'NOMATR', scal=nom, nbret=n3)
!           CALL LXCAPS(NOM)
            if (n1 .ne. 0) then
                call jeexin(jexnum(nomobj, num), iret)
                if (iret .eq. 0) then
                    write (iuni, *) ' L''OBJET : "', num,&
     &              '" DE LA COLLECTION : "', nomobj, '" N''EXISTE PAS'
                else
                    write (iuni, *) ' CONTENU DE L''OBJET : "', num, &
                        '" DE LA COLLECTION : "', nomobj, '"'
                    write (iuni, *) ' '
                    call jeimpa(iuni, jexnum(nomobj, num), txt)
                    call jeimpo(iuni, jexnum(nomobj, num), txt)
                    write (iuni, *) ' '
                end if
            else if (n2 .ne. 0) then
                call jeexin(jexnom(nomobj, nom), iret)
                if (iret .eq. 0) then
                    write (iuni, *) ' L''OBJET : "', nom,&
     &              '" DE LA COLLECTION : "', nomobj, '" N''EXISTE PAS'
                else
                    write (iuni, *) ' CONTENU DE L''OBJET : "', nom, &
                        '" DE LA COLLECTION : "', nomobj, '"'
                    write (iuni, *) ' '
                    call jeimpa(iuni, jexnom(nomobj, nom), txt)
                    call jeimpo(iuni, jexnom(nomobj, nom), txt)
                    write (iuni, *) ' '
                end if
            else if (n3 .ne. 0) then
                noml32 = nomobj
                call jeprat(iuni, noml32, nom(1:8), txt)
                write (iuni, *) ' '
            else
                call jelira(nomobj, 'NMAXOC', noc)
                do i = 1, noc
                    call jeexin(jexnum(nomobj, i), iret)
                    if (iret .ne. 0) then
                        write (iuni, *) ' CONTENU DE L''OBJET : "', i, &
                            '" DE LA COLLECTION : "', nomobj, '"'
                        write (iuni, *) ' '
                        call jeimpa(iuni, jexnum(nomobj, i), txt)
                        call jeimpo(iuni, jexnum(nomobj, i), txt)
                        write (iuni, *) ' '
                    end if
                end do
            end if
        else
            write (iuni, *) ' '
            write (iuni, *) ' CONTENU DE L''OBJET : "', nomobj, '"'
            call jeimpo(iuni, nomobj, txt)
            write (iuni, *) ' '
        end if
        write (iuni, *) ' '
        write (iuni, *) ' FIN DE L''OBJET : "', nomobj, '"'
        write (iuni, *) ' '
!
    else if (noment(1:7) .eq. 'SYSTEME') then
!
        nomcla = ' '
        call getvtx(' ', 'CLASSE', scal=nomcla, nbret=n)
        call getvtx(' ', 'NOMATR', scal=nom, nbret=n3)
        if (n3 .ne. 0) then
            call jeprat(iuni, '$'//nomcla(1:1), nom(1:8), txt)
            write (iuni, *) ' '
        end if
!
    else if (noment(1:8) .eq. 'ATTRIBUT') then
!
        call getvtx(' ', 'NOMOBJ', scal=nomobj, nbret=n)
        call jeexin(nomobj, iret)
        write (iuni, *) ' '
        if (iret .eq. 0) then
            write (iuni, *) ' DIRECTIVE IMPR_JEVEUX '
            write (iuni, *) ' L''OBJET "', nomobj, '" N''EXISTE PAS'
            goto 999
        else
            write (iuni, *) ' '
            write (iuni, *) ' ECRITURE DES ATTRIBUTS DE "', nomobj, '"'
        end if
        write (iuni, *) ' '
        call jeimpa(iuni, nomobj, txt)
        write (iuni, *) ' '
!
    else if (noment(1:14) .eq. 'ENREGISTREMENT') then
!
        nomcla = ' '
        call getvtx(' ', 'CLASSE', scal=nomcla, nbret=n)
        call getvis(' ', 'NUMERO', scal=numerg, nbret=nrg)
        call getvis(' ', 'INFO', scal=info, nbret=nif)
        write (iuni, *) ' '
        call jepreg(fich, nomcla, numerg, txt, info)
        write (iuni, *) ' '
!
    end if
999 continue
end subroutine
