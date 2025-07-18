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

subroutine rsacpa(nomsdz, numva, icode, nomva, ctype, &
                  ival, rval, kval, ier)
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexpa.h"
!
    integer(kind=8) :: numva, icode, ctype, ival(*), ier
    real(kind=8) :: rval(*)
    character(len=80) :: kval(*)
    character(len=16) :: nomva
    character(len=*) :: nomsdz
! ---------------------------------------------------------------------
!  DETERMINE LE NOM D'UNE VARIABLE D'ACCES ET SES VALEURS
!  CONNAISSANT SON NUMERO D'ACCES DANS LA COLLECTION
! ---------------------------------------------------------------------
! IN  NOMSDZ K*  NOM DE LA SD
! IN  NUMVA  I   NUMERO DE LA VARIABLE D'ACCES A RECHERCHER
!                OU 0 POUR LES NUMEROS D'ORDRE  (NUME_ORDRE)
! IN  ICODE  I   CODE POUR FILTRER LES PARAMETRES (VOIR RSEXPA)
! OUT NOMVA  K16 NOM DE LA VARIABLE D'ACCES
! OUT CTYPE  I   LE TYPE : CTYPE (CF. GETCON)
! OUT IVAL   I   LISTE DES VALEURS DU PARAMETRE (CAS ENTIER)
! OUT RVAL   R   LISTE DES VALEURS DU PARAMETRE (CAS REEL)
! OUT IER    I   1 EN CAS D'ERREUR
!                0 SI OK
! ---------------------------------------------------------------------
!
!
    integer(kind=8) :: iret, iord, nbord, i, iad, numord
    character(len=8) :: ktype
    character(len=19) :: nomsd
! ---------------------------------------------------------------------
    call jemarq()
    nomsd = nomsdz
    ier = 0
    nomva = ' '
    ctype = -1
!     TRAITEMENT DE LA DEMANDE PROPRE A NUME_ORDRE
    if (numva .eq. 0) then
        nomva = 'NUME_ORDRE'
        ctype = 2
        call jeveuo(nomsd//'.ORDR', 'L', iord)
        call jelira(nomsd//'.ORDR', 'LONUTI', nbord)
        do i = 1, nbord
            ival(i) = zi(iord-1+i)
        end do
        goto 999
    end if
!     ACCES AU NOM DU CHAMP
    call jenuno(jexnum(nomsd//'.NOVA', numva), nomva)
!     S'AGIT-IL D'UNE VARIABLE D'ACCES
    call rsexpa(nomsd, icode, nomva, iret)
    if (iret .eq. 0) then
        ier = 1
        goto 999
    end if
!     ACCES AUX VALEURS DE LA VARIABLE
    call jeveuo(nomsd//'.ORDR', 'L', iord)
    call jelira(nomsd//'.ORDR', 'LONUTI', nbord)
    do i = 1, nbord
        numord = zi(iord-1+i)
        call rsadpa(nomsd, 'L', 1, nomva, numord, &
                    1, sjv=iad, styp=ktype, istop=0)
!                  123456789.123456789.1234
        kval(i) = '                        '
        if (ktype .eq. 'R') then
!        LES VALEURS SONT REELLES
            ctype = 1
            rval(i) = zr(iad)
        else if (ktype .eq. 'I') then
!        LES VALEURS SONT ENTIERES
            ctype = 2
            ival(i) = zi(iad)
        else if (ktype .eq. 'K8') then
!        LES VALEURS SONT DES CHAINES DE K8
            ctype = 4
            kval(i) = zk8(iad)
        else if (ktype .eq. 'K16') then
!        LES VALEURS SONT DES CHAINES DE K16
            ctype = 5
            kval(i) = zk16(iad)
        else if (ktype .eq. 'K24') then
!        LES VALEURS SONT DES CHAINES DE K24
            ctype = 6
            kval(i) = zk24(iad)
        else if (ktype .eq. 'K32') then
!        LES VALEURS SONT DES CHAINES DE K32
            ctype = 7
            kval(i) = zk32(iad)
        else if (ktype .eq. 'K80') then
!        LES VALEURS SONT DES CHAINES DE K80
            ctype = 8
            kval(i) = zk80(iad)
        else
            ier = 1
            goto 999
        end if
    end do
999 continue
    call jedema()
end subroutine
