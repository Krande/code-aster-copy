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

subroutine chprec(chou)
!     TRAITEMENT DE COMMANDE:   CREA_CHAMP / OPTION: 'EXTR'
!
!     ------------------------------------------------------------------
!
    implicit none
!
!
! 0.1. ==> ARGUMENTS
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterfort/gettco.h"
#include "asterfort/assert.h"
#include "asterfort/chmima.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsinch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "asterfort/chpchd.h"
    character(len=*) :: chou
!
! 0.2. ==> COMMUNS
!
!
! 0.3. ==> VARIABLES LOCALES
!
    character(len=6) :: nompro
    parameter(nompro='CHPREC')
!
    integer(kind=8) :: ibid, icoret, iret, jordr, n1, n2, n3, n4, n5, nbordr, nc, np
    integer(kind=8) :: ifm, niv, iexi
    real(kind=8) :: inst, epsi
    character(len=1) :: base
    character(len=8) :: resuco, interp, crit, proldr, prolga, typmax, cara
    character(len=8) :: nomgd, charme
    character(len=16) :: k16bid, nomcmd, nomch, acces, tysd, tychlu, tych
    character(len=19) :: chextr, noch19, knum, chnos
    character(len=24) :: valk(3)
    character(len=8) :: k8bid, ma, fis
    aster_logical :: grille
!     ------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
!
!
    base = 'G'
    call getres(k8bid, k16bid, nomcmd)
    noch19 = chou
    chnos = "&&CHPREC.CHNOS"
!
    call getvtx(' ', 'NOEUD_CMP', nbval=0, nbret=n1)
    if (n1 .ne. 0 .and. n1 .ne. -2) then
        call utmess('F', 'MODELISA4_16')
    end if
    nomch = ' '
    call getvtx(' ', 'NOM_CHAM', scal=nomch, nbret=n2)
    tychlu = ' '
    call getvtx(' ', 'TYPE_CHAM', scal=tychlu, nbret=n2)
!
!
!     1a. CAS DE LA RECUPERATION DU CHAMP DE GEOMETRIE D'UN MAILLAGE
!     ==============================================================
    if (nomch .eq. 'GEOMETRIE') then
        call getvid(' ', 'MAILLAGE', scal=ma, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'MODELISA4_17')
        end if
!
!       ON VERIFIE QUE LE MOT-CLE TYPE_CHAMP EST COHERENT AVEC LE
!       TYPE DU CHAMP EXTRAIT.
        call dismoi('TYPE_CHAMP', ma//'.COORDO', 'CHAMP', repk=tych)
        call dismoi('NOM_GD', ma//'.COORDO', 'CHAMP', repk=nomgd)
!
        if (tychlu(6:) .ne. nomgd) then
            valk(1) = tychlu
            valk(2) = tych(1:4)
            valk(3) = nomgd
            call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
        end if

        if (tychlu(1:4) .eq. "GEOM") then
            call copisd('CHAMP_GD', 'G', ma//'.COORDO', noch19)
        elseif (tychlu(1:4) .eq. "NOEU") then
            call chpchd(ma//'.COORDO', "NOEU", ma, "NON", "G", noch19)
        else
            valk(1) = tychlu
            valk(2) = tych(1:4)
            valk(3) = nomgd
            call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
        end if
        goto 20
    end if
!
!
!     1b. CAS DE LA RECUPERATION DU CHAMP D'ABSC_CURV D'UN MAILLAGE
!     ==============================================================
    if (nomch .eq. 'ABSC_CURV') then
        call getvid(' ', 'MAILLAGE', scal=ma, nbret=n1)
        if (n1 .eq. 0) then
            call utmess('F', 'MODELISA4_17')
        end if
!
!       ON VERIFIE QUE LE MOT-CLE TYPE_CHAMP EST COHERENT AVEC LE
!       TYPE DU CHAMP EXTRAIT.
        call dismoi('TYPE_CHAMP', ma//'.ABSC_CURV', 'CHAMP', repk=tych)
        call dismoi('NOM_GD', ma//'.ABSC_CURV', 'CHAMP', repk=nomgd)
!
        if ((tychlu(1:4) .ne. tych) .or. (tychlu(6:) .ne. nomgd)) then
            valk(1) = tychlu
            valk(2) = tych(1:4)
            valk(3) = nomgd
            call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
        end if
        call copisd('CHAMP_GD', 'G', ma//'.ABSC_CURV', noch19)
        goto 20
    end if
!
!
!     2. CAS DE LA RECUPERATION D'UN CHAMP DANS UNE SD FISS_XFEM
!     ==============================================================
    call getvid(' ', 'FISSURE', scal=fis, nbret=n1)
    if (n1 .eq. 1) then
!
!       VERIFIE SI UNE GRILLE AUXILIAIRE EST DEFINIE POUR LA FISS
        call jeexin(fis//'.GRI.MAILLA', ibid)
        if (ibid .eq. 0) then
            grille = .false.
        else
            grille = .true.
        end if
!
        if (nomch .eq. 'LTNO') then
            chextr = fis//'.LTNO'
        else if (nomch .eq. 'LNNO') then
            chextr = fis//'.LNNO'
        else if (nomch .eq. 'GRLNNO') then
            chextr = fis//'.GRLNNO'
        else if (nomch .eq. 'GRLTNO') then
            chextr = fis//'.GRLTNO'
        else if (nomch .eq. 'STNO') then
            chextr = fis//'.STNO'
        else if (nomch .eq. 'STNOR') then
            chextr = fis//'.STNOR'
        else if (nomch .eq. 'BASLOC') then
            chextr = fis//'.BASLOC'
        else
            if (grille) then
                if (nomch .eq. 'GRI.LTNO') then
                    chextr = fis//'.GRI.LTNO'
                else if (nomch .eq. 'GRI.LNNO') then
                    chextr = fis//'.GRI.LNNO'
                else if (nomch .eq. 'GRI.GRLNNO') then
                    chextr = fis//'.GRI.GRLNNO'
                else if (nomch .eq. 'GRI.GRLTNO') then
                    chextr = fis//'.GRI.GRLTNO'
                end if
            else
                call utmess('F', 'XFEM2_98')
            end if
        end if
!
!       ON VERIFIE QUE LE MOT-CLE TYPE_CHAMP EST COHERENT AVEC LE
!       TYPE DU CHAMP EXTRAIT.
        call dismoi('TYPE_CHAMP', chextr, 'CHAMP', repk=tych)
        call dismoi('NOM_GD', chextr, 'CHAMP', repk=nomgd)
!
        if ((tychlu(1:4) .ne. tych) .or. (tychlu(6:) .ne. nomgd)) then
            valk(1) = tychlu
            valk(2) = tych(1:4)
            valk(3) = nomgd
            call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
        end if
        call copisd('CHAMP_GD', 'G', chextr, noch19)
        goto 20
    end if
!
!
!     3. CAS DE LA RECUPERATION D'UN CHAMP DANS UNE SD CARA_ELEM
!     ==============================================================
    call getvid(' ', 'CARA_ELEM', scal=cara, nbret=n1)
    if (n1 .eq. 1) then
!
        ASSERT(nomch(1:1) .eq. '.')
        chextr = cara//nomch(1:11)
        call exisd('CHAMP', chextr, iexi)
        if (iexi .eq. 0) then
            call utmess('F', 'CALCULEL3_17', sk=chextr)
        end if
!
!       ON VERIFIE QUE LE MOT-CLE TYPE_CHAMP EST COHERENT AVEC LE
!       TYPE DU CHAMP EXTRAIT.
        call dismoi('TYPE_CHAMP', chextr, 'CHAMP', repk=tych)
        call dismoi('NOM_GD', chextr, 'CHAMP', repk=nomgd)
!
        if ((tychlu(1:4) .ne. tych) .or. (tychlu(6:) .ne. nomgd)) then
            valk(1) = tychlu
            valk(2) = tych(1:4)
            valk(3) = nomgd
            call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
        end if
        call copisd('CHAMP_GD', 'G', chextr, noch19)
        goto 20
    end if
!
!
!     4. CAS DE LA RECUPERATION D'UN CHAMP DANS UNE SD CHAR_MECA
!     ==============================================================
    call getvid(' ', 'CHARGE', scal=charme, nbret=n1)
    if (n1 .eq. 1) then
!
        ASSERT(nomch(1:1) .eq. '.')
        chextr = charme//nomch(1:11)
        call exisd('CHAMP', chextr, iexi)
        if (iexi .eq. 0) then
            call utmess('F', 'CALCULEL3_17', sk=chextr)
        end if
!
!       ON VERIFIE QUE LE MOT-CLE TYPE_CHAMP EST COHERENT AVEC LE
!       TYPE DU CHAMP EXTRAIT.
        call dismoi('TYPE_CHAMP', chextr, 'CHAMP', repk=tych)
        call dismoi('NOM_GD', chextr, 'CHAMP', repk=nomgd)
!
        if ((tychlu(1:4) .ne. tych) .or. (tychlu(6:) .ne. nomgd)) then
            valk(1) = tychlu
            valk(2) = tych(1:4)
            valk(3) = nomgd
            call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
        end if
        call copisd('CHAMP_GD', 'G', chextr, noch19)
        goto 20
    end if
!
!
!     5. CAS DE LA RECUPERATION D'UN CHAMP D'UNE SD RESULTAT
!     ==============================================================
    call getvid(' ', 'RESULTAT', scal=resuco, nbret=n1)
    interp = ' '
    call getvtx(' ', 'INTERPOL', scal=interp, nbret=n3)
    typmax = ' '
    call getvtx(' ', 'TYPE_MAXI', scal=typmax, nbret=n5)
    call gettco(resuco, tysd)
!
!     --- ON PEUT FAIRE UNE INTERPOLATION ---
!         ===============================
    if (tysd .eq. 'EVOL_THER' .or. tysd .eq. 'EVOL_ELAS' .or. tysd .eq. 'EVOL_NOLI' .or. &
        tysd .eq. 'DYNA_TRANS' .or. tysd .eq. 'EVOL_VARC' .or. tysd .eq. 'EVOL_SECH') then
!
        if (interp(1:3) .eq. 'LIN') then
            call getvr8(' ', 'INST', scal=inst, nbret=n4)
            ASSERT(n4 .eq. 1)
            proldr = 'EXCLUS'
            prolga = 'EXCLUS'
            acces = 'INST'
            call getvr8(' ', 'PRECISION', scal=epsi, nbret=np)
            ASSERT(np .eq. 1)
            call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
            ASSERT(nc .eq. 1)
            call rsinch(resuco, nomch, acces, inst, noch19, &
                        proldr, prolga, 2, base, epsi, crit, icoret)
        else
            if (n5 .ne. 0) then
                call chmima(resuco, nomch, tychlu, typmax, noch19)
            else
                knum = '&&'//nompro//'.NUME_ORDRE'
                call getvr8(' ', 'PRECISION', scal=epsi, nbret=np)
                call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
                call rsutnu(resuco, ' ', 0, knum, nbordr, &
                            epsi, crit, iret)
                if ((iret .ne. 0) .or. (nbordr .gt. 1)) goto 10
                if (nbordr .eq. 0) then
                    call utmess('F', 'UTILITAI_23')
                end if
                call jeveuo(knum, 'L', jordr)
                call rsexch('F', resuco, nomch, zi(jordr), chextr, &
                            iret)
!
!
!           ON VERIFIE QUE LE MOT-CLE TYPE_CHAMP EST COHERENT AVEC LE
!           TYPE DU CHAMP EXTRAIT.
!
                call dismoi('TYPE_CHAMP', chextr, 'CHAMP', repk=tych)
                call dismoi('NOM_GD', chextr, 'CHAMP', repk=nomgd)
!
                if ((tychlu(1:4) .ne. tych) .or. (tychlu(6:) .ne. nomgd)) then
                    valk(1) = tychlu
                    valk(2) = tych(1:4)
                    valk(3) = nomgd
                    call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
                end if
                call copisd('CHAMP_GD', 'G', chextr, noch19)
                call jedetr(knum)
            end if
        end if
!
!     --- ON NE FAIT QU'UNE EXTRACTION ---
!         ===========================
    else
        if (interp(1:3) .eq. 'LIN') then
            valk(1) = tysd
            call utmess('F', 'MODELISA8_55', sk=valk(1))
        elseif (n5 .ne. 0) then
            call chmima(resuco, nomch, tychlu, typmax, noch19)
        else
            knum = '&&'//nompro//'.NUME_ORDRE'
            call getvr8(' ', 'PRECISION', scal=epsi, nbret=np)
            call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
            call rsutnu(resuco, ' ', 0, knum, nbordr, &
                        epsi, crit, iret)
            if ((iret .ne. 0) .or. (nbordr .gt. 1)) goto 10
            if (nbordr .eq. 0) then
                call utmess('F', 'UTILITAI_23')
            end if
            call jeveuo(knum, 'L', jordr)
            call rsexch('F', resuco, nomch, zi(jordr), chextr, &
                        iret)
            call dismoi('TYPE_CHAMP', chextr, 'CHAMP', repk=tych)
            call dismoi('NOM_GD', chextr, 'CHAMP', repk=nomgd)
!
            if ((tychlu(1:4) .ne. tych) .or. (tychlu(6:) .ne. nomgd)) then
                valk(1) = tychlu
                valk(2) = tych(1:4)
                valk(3) = nomgd
                call utmess('F', 'MODELISA4_18', nk=3, valk=valk)
            end if
!
            call copisd('CHAMP_GD', 'G', chextr, noch19)
            call jedetr(knum)
        end if
    end if
!
!
    goto 20
10  continue
    call utmess('F', 'MODELISA4_19')
!
20  continue

    call titre()
!
!
! --- MENAGE
!     ======
    call jedetr('&&'//nompro//'.NUME_ORDRE')
!
    call jedema()
end subroutine
