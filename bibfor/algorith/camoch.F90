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

subroutine camoch(nomres, numref, intf, raid, mass, &
                  raildl, inord)
    implicit none
!
!***********************************************************************
!    P. RICHARD     DATE 19/02/91
!-----------------------------------------------------------------------
!  BUT:  CALCUL DES MODES CONTRAINTS (DEPLACEMENT UNITAIRE IMPOSE)
!      ET STOCKAGE DANS LE CONCEPT MODE MECA A PARTIR D'UN
!                  NUMERO D'ORDRE
!-----------------------------------------------------------------------
!
! NOMRES    /I/: NOM DU CONCEPT RESULTAT
! NUMREF    /I/: NOM UT DU NUM_DDL DE REFERENCE
! INTF      /I/: NOM UT DE L'INTERF_DYNA EN AMONT
! RAID      /I/: NOM DE LA MATRICE RAIDEUR
! MASS      /I/: NOM DE LA MATRICE MASS
! RAILDL    /M/: NOM DE LA MATRICE RAIDEUR FACTORISEE LDLT OU BLANC
! INORD     /M/: DERNIER NUMERO D'ORDRE UTILISE
!
!
!
!
!
#include "jeveux.h"
#include "asterfort/cheddl.h"
#include "asterfort/defsta.h"
#include "asterfort/dismoi.h"
#include "asterfort/facmtr.h"
#include "asterfort/isdeco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mtcmbl.h"
#include "asterfort/mtdefs.h"
#include "asterfort/mtdscr.h"
#include "asterfort/wkvect.h"
#include "asterfort/int_to_char8.h"
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iad, ier, ik, ino, inord
    integer(kind=8) :: j, jj, lldeeq, lldes, llncmp
    integer(kind=8) :: llnoin, lmat, ltddl, ltpar, nbcb, nbcmp
    integer(kind=8) :: nbcont, nbcpmx, nbdeb, nbec, nbfin, nbint, nbnoe
    integer(kind=8) :: nbnot, neq, ntail1, ntail2, numgd
    real(kind=8) :: freq, omega, pi
!-----------------------------------------------------------------------
    parameter(nbcpmx=300)
    character(len=1) :: typcst(2)
    character(len=6) :: pgc
    character(len=8) :: nomddl
    character(len=8) :: nomres, intf, typcou, nomnoe, nomcmp, mailla
    character(len=19) :: numref, numddl
    character(len=19) :: raildl, raid, mass
    character(len=16) :: typdef
    character(len=24) :: desdef, deeq, temddl, tempar
    real(kind=8) :: const(2)
    integer(kind=8) :: idec(nbcpmx)
    character(len=24) :: nmat(2), nlmat
    character(len=8), pointer :: idc_type(:) => null()
    real(kind=8), pointer :: idc_dy_freq(:) => null()
!
!-----------------------------------------------------------------------
    data nomddl/'        '/
    data pgc/'CAMOCH'/
!-----------------------------------------------------------------------
!
!
    call jemarq()
    pi = 4.d0*atan2(1.d0, 1.d0)
    typdef = 'CONT-HARM'
!
!---------------------RECHERCHE DU NUMDDL ASSOCIE A LA MATRICE----------
!
    call dismoi('NOM_NUME_DDL', raid, 'MATR_ASSE', repk=numddl)
!
!---------------------REQUETTE DU DEEQ DU NUMDDL------------------------
!
    numddl(15:19) = '.NUME'
    deeq = numddl//'.DEEQ'
    call jeveuo(deeq, 'L', lldeeq)
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!
!--------------------RECUPERATION DU MAILLAGE---------------------------
!
    call dismoi('NOM_MAILLA', numddl, 'NUME_DDL', repk=mailla)
!
!----RECUPERATION DES DONNEES RELATIVES A LA GRANDEUR SOUS-JACENTE------
!
    call dismoi('NB_CMP_MAX', intf, 'INTERF_DYNA', repi=nbcmp)
    call dismoi('NB_EC', intf, 'INTERF_DYNA', repi=nbec)
    call dismoi('NUM_GD', intf, 'INTERF_DYNA', repi=numgd)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', numgd), 'L', llncmp)
!
!
!-----------REQUETTE ADRESSE DE LA TABLE DESCRIPTION DES DEFORMEES------
!
    desdef = intf//'.IDC_DEFO'
    call jeveuo(desdef, 'L', lldes)
    call jelira(desdef, 'LONMAX', nbnot)
!**************************************************************
    nbnot = nbnot/(2+nbec)
!      NBNOT=NBNOT/3
!**************************************************************
!
!
!------------REQUETTE ADRESSE DEFINITION INTERFACE ET TYPE--------------
!
    call jelira(intf//'.IDC_LINO', 'NMAXOC', nbint)
    call jeveuo(intf//'.IDC_TYPE', 'L', vk8=idc_type)
!
!-----------COMPTAGE DU NOMBRE DE NOEUDS CB-HARM------------------------
!
    nbdeb = nbnot
    nbfin = 0
!
    do j = 1, nbint
        call jelira(jexnum(intf//'.IDC_LINO', j), 'LONMAX', nbnoe)
        typcou = idc_type(j)
        if (typcou .eq. 'CB_HARMO') then
            call jeveuo(jexnum(intf//'.IDC_LINO', j), 'L', llnoin)
            do i = 1, nbnoe
                ik = zi(llnoin+i-1)
                nbfin = max(nbfin, ik)
                nbdeb = min(nbdeb, ik)
            end do
            call jelibe(jexnum(intf//'.IDC_LINO', j))
        end if
    end do
!
    call jelibe(intf//'.IDC_TYPE')
!
    if (nbfin .gt. 0) then
        nbcb = nbfin-nbdeb+1
    else
        nbcb = 0
    end if
!
!
!----------ALLOCATION DU VECTEUR DES DDL A IMPOSER A 1------------------
!
    ntail1 = (nbcb*nbcmp)*2
    ntail2 = (nbcb*nbcmp)
!
!  TAILLE DOUBLE CAR PRESENCE EVENTUELLE DE DOUBLE LAGRANGE POUR LE
!   BLOCAGE
!
    if (ntail1 .eq. 0) goto 999
!
!
    call jeveuo(intf//'.IDC_DY_FREQ', 'L', vr=idc_dy_freq)
    freq = idc_dy_freq(1)
    call jelibe(intf//'.IDC_DY_FREQ')
    omega = 2.d0*pi*freq
!
    temddl = '&&'//pgc//'.LISTE.DDL'
    tempar = '&&'//pgc//'.PARA.NOCMP'
    call wkvect(temddl, 'V V I', ntail1, ltddl)
    call wkvect(tempar, 'V V K16', ntail2, ltpar)
    if (raildl .eq. '                  ') then
        call mtdscr(raid)
        nmat(1) = raid(1:19)//'.&INT'
        call mtdscr(mass)
        nmat(2) = mass(1:19)//'.&INT'
        raildl = '&&'//pgc//'.RAID.HARM'
        call mtdefs(raildl, raid, 'V', 'R')
        call mtdscr(raildl)
        nlmat = raildl(1:19)//'.&INT'
        call jeveuo(nlmat, 'E', lmat)
        typcst(1) = 'R'
        typcst(2) = 'R'
        const(1) = 1.d0
        const(2) = -(omega**2)
        call mtcmbl(2, typcst, const, nmat, nlmat, &
                    nomddl, ' ', 'ELIM=')
        call facmtr(raildl, raildl, ier)
    end if
!
!-------------COMPTAGE ET REPERAGE DES DEFORMEES A CALCULER-------------
!
    nbcont = 0
!
    if (nbcb .gt. 0) then
        do i = nbdeb, nbfin
!**************************************************************
!          ICOD=ZI(LLDES+2*NBNOT+I-1)
            call isdeco(zi(lldes+2*nbnot+(i-1)*nbec+1-1), idec, nbcmp)
!**************************************************************
            ino = zi(lldes+i-1)
            nomnoe = int_to_char8(ino)
            do j = 1, nbcmp
                if (idec(j) .eq. 1) then
                    jj = -j
                    nbcont = nbcont+1
                    nomcmp = zk8(llncmp+j-1)
                    zk16(ltpar+nbcont-1) = nomnoe//nomcmp
                    iad = ltddl+(nbcont-1)*2
                    call cheddl(zi(lldeeq), neq, ino, jj, zi(iad), &
                                2)
                end if
            end do
        end do
    end if
!
!
!
!------------------CALCUL DES MODES CONTRAINTS--------------------------
!
    call defsta(nomres, numref, raildl, zi(ltddl), zk16(ltpar), &
                2, nbcont, typdef, inord)
!
!----------------------MENAGE-------------------------------------------
!
    call jedetr(temddl)
    call jedetr(tempar)
    call jelibe(deeq)
!
!
999 continue
    call jedema()
end subroutine
