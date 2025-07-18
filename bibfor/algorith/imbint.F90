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

subroutine imbint(nomres, ifm)
!    P. RICHARD     DATE 21/02/1991
!-----------------------------------------------------------------------
!  BUT:  IMPRIMER LES RESULTATS RELATIFS A LA BASE MODALE
    implicit none
!
!-----------------------------------------------------------------------
!
! NOMRES   /I/: NOM DU CONCEPT RESULTAT
! MAILLA   /I/: NOM DU MAILLA
! IFM      /I/: UMITE DU FICHIER MESSAGE
!
!
!
!
!
#include "jeveux.h"
#include "asterfort/bmnodi.h"
#include "asterfort/dismoi.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/int_to_char8.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ibid(1), idau, idcb, idda, idha, idmn
    integer(kind=8) :: ifau, ifcb, ifha, ifmn, ino, ipoin
    integer(kind=8) :: ityp, j, k, llact, lldes, llncmp
    integer(kind=8) :: llnoe, lltyp, nbcmp, nbcpmx, nbdef, nbec, nbint
    integer(kind=8) :: nbno, nbnot, ncomp, numgd
!-----------------------------------------------------------------------
    parameter(nbcpmx=300)
    character(len=8) :: nomcou, typcou, nomnoe, typ, typbas(3), nomtyp
    character(len=8) :: nomres, mailla, flec, craigb, mneal, aucun, cbharm
    character(len=16) :: tydef
    character(len=11) :: dactif
    character(len=24) :: nomint, typint, noeint, desdef, ddact
    character(len=80) :: chaine
    integer(kind=8) :: idec(nbcpmx), ifm, i1
!
    integer(kind=8) :: ibid1
    integer(kind=8), pointer :: idc_desc(:) => null()
    data ibid1/0/
!
!-----------------------------------------------------------------------
    data mneal, craigb, aucun, flec/'MNEAL', 'CRAIGB', 'AUCUN', '--->'/
    data cbharm/'CB_HARMO'/
    data dactif/'DDL_ACTIF: '/
    data typbas/'CONNUE', 'CYCLIQUE', 'RITZ'/
!-----------------------------------------------------------------------
!
    call jemarq()
!
    write (ifm, *) ' '
    write (ifm, *) '----------------------------------------------------'
    write (ifm, *) ' '
    write (ifm, *) '                DEFI_INTERF_DYNA '
    write (ifm, *) ' '
    write (ifm, *) '  IMPRESSIONS NIVEAU: 2 '
    write (ifm, *) ' '
!
!--------------RECUPERATION DU NOM DU MAILLA--------------------------
!
    call dismoi('NOM_MAILLA', nomres, 'INTERF_DYNA', repk=mailla)
!
!--------------RECUPERATION TYPE LIST_INTERFACE-------------------------
!
    call jeveuo(nomres//'.IDC_DESC', 'L', vi=idc_desc)
    ityp = idc_desc(1)
    call jelibe(nomres//'.IDC_DESC')
!
    nomtyp = typbas(ityp)
!
!----RECUPERATION DES DONNEES RELATIVES A LA GRANDEUR SOUS-JACENTE------
!
    call dismoi('NB_CMP_MAX', nomres, 'INTERF_DYNA', repi=nbcmp)
    call dismoi('NB_EC', nomres, 'INTERF_DYNA', repi=nbec)
    call dismoi('NUM_GD', nomres, 'INTERF_DYNA', repi=numgd)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', numgd), 'L', llncmp)
!
!--------------------INITIALISATION DES MOTS USUELS---------------------
!
    desdef = nomres//'.IDC_DEFO'
    nomint = nomres//'.IDC_NOMS'
    typint = nomres//'.IDC_TYPE'
    noeint = nomres//'.IDC_LINO'
    ddact = nomres//'.IDC_DDAC'
!
    call jeveuo(typint, 'L', lltyp)
    call jeveuo(desdef, 'L', lldes)
    call jelira(nomint, 'NOMMAX', nbint)
    call jelira(desdef, 'LONMAX', nbnot)
    nbnot = nbnot/(2+nbec)
!
    idau = nbnot+1
    idmn = nbnot+1
    idcb = nbnot+1
    idha = nbnot+1
    ifcb = 0
    ifmn = 0
    ifau = 0
    ifha = 0
!
    write (ifm, *) ' '
    write (ifm, *) ' NOM DE L'' INTERF_DYNA: ', nomres
    write (ifm, *) '-------------------------- '
!
    write (ifm, *) ' '
    write (ifm, *) ' TYPE : ', nomtyp
    write (ifm, *) '------ '
!
    write (ifm, *) ' '
    write (ifm, *) ' DEFINITION DES INTERFACES'
    write (ifm, *) '--------------------------- '
!
!  BOUCLE SUR LES INTERFACES
!
    do i = 1, nbint
        write (ifm, *) ' '
        write (ifm, *) ' '
        typcou = zk8(lltyp+i-1)
        call jeveuo(jexnum(ddact, i), 'L', llact)
        call jenuno(jexnum(nomint, i), nomcou)
        call jeveuo(jexnum(noeint, i), 'L', llnoe)
        call jelira(jexnum(noeint, i), 'LONMAX', nbno)
        write (ifm, *) ' INTERFACE: ', nomcou
        write (ifm, *) '---------- '
        write (ifm, *) '              TYPE: ', typcou
        write (ifm, *) ' '
        write (ifm, *) ' LISTE DES NOEUDS:  NOMBRE: ', nbno
        write (ifm, *) ' '
!
!  BOUCLE SUR LES NOEUDS DE L'INTERFACE COURANTE
!
        do j = 1, nbno
            ipoin = zi(llnoe+j-1)
            call isdeco(zi(llact+(j-1)*nbec+1-1), idec, nbcmp)
            idda = 1
            do k = 1, nbcmp
                if (idec(k) .gt. 0) then
                    chaine(idda:idda+7) = zk8(llncmp+k-1)
                    idda = idda+8
                end if
            end do
            idda = idda-1
            ino = zi(lldes+ipoin-1)
            nomnoe = int_to_char8(ino)
            if (idda .lt. 1) then
                write (ifm, *) 'NOEUD: ', j, flec, nomnoe, ' ', dactif, &
                    'PAS DE DDL ACTIF'
                goto 100
            end if
            write (ifm, *) 'NOEUD: ', j, flec, nomnoe, ' ', dactif, chaine(1: &
                                                                           idda)
100         continue
!
!  STOCKAGE DU NUMERO DU PREMIER ET DERNIER NOEUD DE CHAQUE TYPE
!            D'INTERFACE
!
            if (typcou .eq. mneal) then
                idmn = min(idmn, ipoin)
                ifmn = max(ifmn, ipoin)
            end if
!
            if (typcou .eq. craigb) then
                idcb = min(idcb, ipoin)
                ifcb = max(ifcb, ipoin)
            end if
!
            if (typcou .eq. cbharm) then
                idha = min(idha, ipoin)
                ifha = max(ifha, ipoin)
            end if
!
            if (typcou .eq. aucun) then
                idau = min(idau, ipoin)
                ifau = max(ifau, ipoin)
            end if
!
        end do
        write (ifm, *) ' '
        i1 = i
        call bmnodi(' ', nomres, ' ', i1, 0, &
                    ibid(1), nbdef)
        write (ifm, *) '  '
        write (ifm, *) ' NOMBRE DE DEFORMEES STATIQUES ASSOCIES: ', nbdef
        write (ifm, *) '  '
        call jelibe(jexnum(ddact, i))
        call jelibe(jexnum(noeint, i))
    end do
!
    write (ifm, *) ' '
    write (ifm, *) ' '
    write (ifm, *) ' '
    write (ifm, *) ' DEFINITION DES DEFORMEES A CALCULER'
    write (ifm, *) '------------------------------------'
!
!  CAS OU IL Y A DES DEFORMEES STATIQUES
!
    if (idau .eq. 1) then
        write (ifm, *) ' PAS DE DEFORMEES STATIQUES A CALCULER'
        goto 999
    end if
    write (ifm, *) ' '
    ncomp = 0
    do i = 1, nbnot
        write (ifm, *) ' '
        if (i .ge. idmn .and. i .le. ifmn) tydef = 'MODE D''ATTACHE'
        if (i .ge. idcb .and. i .le. ifcb) tydef = 'MODE CONTRAINT'
        if (i .ge. idha .and. i .le. ifha) tydef = 'MODE CONT-HARM'
        ino = zi(lldes+i-1)
        nomnoe = int_to_char8(ino)
        call isdeco(zi(lldes+nbnot*2+(i-1)*nbec+1-1), idec, nbcmp)
        do j = 1, nbcmp
            if (idec(j) .gt. 0) then
                typ = zk8(llncmp+j-1)
                ncomp = ncomp+1
                write (ifm, *) 'DEFORMEE: ', ncomp, flec, nomnoe, ' ', typ, ' ', tydef
            end if
        end do
    end do
!
    write (ifm, *) ' '
    write (ifm, *) '----------------------------------------------------'
    write (ifm, *) ' '
!
999 continue
    call jedema()
end subroutine
