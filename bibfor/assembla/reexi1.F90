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
subroutine reexi1(nu, mo, ma, nlili, nm, &
                  nl, nbntt)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "jeveux.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: mo, ma
    character(len=14) :: nu
    integer(kind=8) :: nlili, nm, nl, nbntt
! ----------------------------------------------------------------------
!     BUT:  CETTE ROUTINE SERT A :
!
!           1) RENDRE : MO, MA, NLILI, NM, NL, NBNTT
!
!           2) RENDRE DANS LE NUME_DDL L' OBJET .EXI1
!              (POUR DIRE QUELS SONT LES NOEUDS (TARDIFS OU NON)
!               INTERVENANT REELLEMENT DANS LA NUMEROTATION
!               PLUS L'EXISTENCE D'UN NOEUD FACTICE 0 )
!             ATTENTION : CET OBJET A UN CONTENU PLUS RICHE QUE
!                         NECESSAIRE EN PREVISION DE 'RCMK'
!
!     IN:
!     ---
!       NU : NOM DU NUME_DDL  AUQUEL ON VA AJOUTER  L' OBJET .EXI1
!           (ON SE SERT EN ENTREE DU SEUL OBJET NU//'.NUME.LILI')
!
!     OUT:
!     ---- MO : NOM DU MODELE SOUS-JACENT AU NUME_DDL
!          MA : NOM DU MAILLAGE SOUS-JACENT AU NUME_DDL
!          NLILI: NOMBRE DE LIGREL DE L'OBJET .LILI
!
!         SOIT NM LE NOMBRE DE NOEUDS PHYSIQUES DU MAILLAGE (LILI(1))
!              NL LE NOMBRE DE NOEUDS TARDIFS DU MAILLAGE
!              N2 LE NOMBRE DE NOEUDS TARDIFS DU MODELE (LILI(2))
!              N3 LE NOMBRE DE NOEUDS TARDIFS DE LA 1ERE CHARGE(LILI(3))
!              .....
!              NP LE NOMBRE DE NOEUDS TARDIFS DE LA DERE CHARGE(LILI(P))
!
!          NBNOM = NM+NL      (NOMBRE MAX DE NOEUDS DU MAILLAGE)
!          NBNOT = N2+...+NP  (NOMBRE MAX DE NOEUDS TARDIFS DU MODELE
!                              ET DE LA LISTE DE CHARGES)
!          NBNTT = NBNOM+NBNOT
!
!        NU EST COMPLETE PAR .EXI1
!        -------------------------
!        .EXI1(*) EST DIMENSIONNE A NBNTT+1
!
!        SOIT LA NUMEROTATION IMPLICITE TOTALE :
!          1- LE NOEUD FACTICE 0
!          1- LES NOEUDS PHYSIQUES DU MAILLAGE (NI)
!          3- LES NOEUDS TARDIFS DU MAILLAGE (&I)   /
!          4- LES NOEUDS TARDIFS DU MODELE   (&LMI)
!          5- LES NOEUDS TARDIFS DE LA CHARGE 1 (&LCH1I)
!           - ...
!           - LES NOEUDS TARDIFS DE LA CHARGE P (&LCH1P)
!
!        .EXI1(0)=1
!        POUR I=1,NBNTT
!        .EXI1(I) >0 SI LE NOEUD I EXISTE DANS LE NUME_DDL
!                 =0 SINON
!                  (EN FAIT EXI1(I) CONTIENT UN MAJORANT DU NOMBRE
!                  DE NOEUDS CONNECTES AU NOEUD I. CET OBJET SERVIRA
!                  DANS 'RCMK' A ALLOUER LA TABLE DE CONNECTIVITE
!                  INVERSE)
! ----------------------------------------------------------------------
!     VARIABLES LOCALES:
!     ------------------
    character(len=8) :: exiele
    character(len=24) :: nomli2, mo2
    character(len=19) :: nomlig
!
!
!     -- RECUPERATION DU NOM DU MODELE SOUS-JACENT A LA NUMEROTATION :
!     ----------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iaconx, iaexi1, iagrel, ialiel, iamail, ianbno
    integer(kind=8) :: ianema, iel, igrel
    integer(kind=8) :: iino, ilconx, ili, illiel, ilnema, ima, ino
    integer(kind=8) :: iret, j, jjno, jno, nbel, nbgrel, nbnm
    integer(kind=8) :: nbnom, nbnot, nbsma, nbssa, nma
    integer(kind=8), pointer :: sssa(:) => null()
!-----------------------------------------------------------------------
    call jemarq()
    mo = ' '
    call jelira(nu//'.NUME.LILI', 'NOMMAX', nlili)
    do i = 1, nlili
        call jenuno(jexnum(nu//'.NUME.LILI', i), mo2)
        if (mo2(9:15) .eq. '.MODELE') then
            mo = mo2(1:8)
            goto 42
        end if
    end do
!    call utmess('F', 'ASSEMBLA_35', sk=nu)
!
42  continue
    if (mo .ne. ' ') then
        call dismoi('NOM_MAILLA', mo, 'MODELE', repk=ma)
    else
        call jenuno(jexnum(nu//'.NUME.LILI', nlili), nomli2)
        nomlig = nomli2(1:19)
        call dismoi('NOM_MAILLA', nomlig, 'LIGREL', repk=ma)
    end if
!
!     -- CALCUL DE NBNOM:
!     -------------------
    call dismoi('NB_NO_MAILLA', ma, 'MAILLAGE', repi=nm)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nma)
    nl = 0
    if (mo .ne. ' ') then
        call dismoi('NB_NL_MAILLA', mo, 'MODELE', repi=nl)
    end if
    nbnom = nm+nl
!
!     -- CALCUL DE NBNOT ET NBNTT :
!     -----------------------------
    nbnot = 0
    do 1, ili = 2, nlili
        call jenuno(jexnum(nu//'.NUME.LILI', ili), nomli2)
        nomlig = nomli2(1:19)
        call jeveuo(nomlig//'.NBNO', 'L', ianbno)
        nbnot = nbnot+zi(ianbno)
1   end do
    nbntt = nbnom+nbnot
!
!
!     -- ALLOCATION DE .EXI1
!     ----------------------
    call wkvect(nu//'.EXI1', 'V V I', nbntt+1, iaexi1)
!
    zi(iaexi1) = 1
!
    if (nma .gt. 0) then
        call jeveuo(ma//'.CONNEX', 'L', iaconx)
        call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilconx)
    end if
!
!
!     -- 1ERE ETAPE : (SUPER)MAILLES DU MAILLAGE:
!     -------------------------------------------
    nbssa = 0
    nbsma = 0
    if (mo .ne. ' ') then
        call dismoi('NB_SS_ACTI', mo, 'MODELE', repi=nbssa)
        call dismoi('NB_SM_MAILLA', mo, 'MODELE', repi=nbsma)
    end if
    if (nbssa .gt. 0) then
        call jeveuo(mo//'.MODELE    .SSSA', 'L', vi=sssa)
    else
        goto 22
    end if
!
    do ima = 1, nbsma
        if (sssa(ima) .eq. 1) then
            call jeveuo(jexnum(ma//'.SUPMAIL', ima), 'L', iamail)
            call jelira(jexnum(ma//'.SUPMAIL', ima), 'LONMAX', nbnm)
            do i = 1, nbnm
                ino = zi(iamail-1+i)
                if (ino .eq. 0) then
                    call utmess('F', 'ASSEMBLA_36')
                end if
                do j = i+1, nbnm
                    jno = zi(iamail-1+j)
                    zi(iaexi1+ino) = zi(iaexi1+ino)+1
                    zi(iaexi1+jno) = zi(iaexi1+jno)+1
                end do
            end do
        end if
    end do
!
22  continue
!
!     -- 2EME ETAPE : MAILLES TARDIVES (OU NON) DES LIGRELS:
!                     (MODELE + LISTE DE CHARGES)
!     ------------------------------------------------------
    nbnot = 0
    do 31, ili = 2, nlili
        call jenuno(jexnum(nu//'.NUME.LILI', ili), nomli2)
        nomlig = nomli2(1:19)
        call dismoi('EXI_ELEM', nomlig, 'LIGREL', repk=exiele)
        if (exiele(1:3) .eq. 'NON') goto 31
!
        call jeveuo(nomlig//'.LIEL', 'L', ialiel)
        call jeveuo(jexatr(nomlig//'.LIEL', 'LONCUM'), 'L', illiel)
        call jelira(nomlig//'.LIEL', 'NMAXOC', nbgrel)
!
!
        call jeexin(nomlig//'.NEMA', iret)
        if (iret .gt. 0) then
            call jeveuo(nomlig//'.NEMA', 'L', ianema)
            call jeveuo(jexatr(nomlig//'.NEMA', 'LONCUM'), 'L', ilnema)
        end if
!
        do igrel = 1, nbgrel
            nbel = zi(illiel-1+igrel+1)-zi(illiel-1+igrel)-1
            iagrel = ialiel+zi(illiel-1+igrel)-1
!
            do iel = 1, nbel
                ima = zi(iagrel-1+iel)
                if (ima .gt. 0) then
                    nbnm = zi(ilconx-1+ima+1)-zi(ilconx-1+ima)
                    iamail = iaconx+zi(ilconx-1+ima)-1
                else
                    nbnm = zi(ilnema-1-ima+1)-zi(ilnema-1-ima)-1
                    iamail = ianema+zi(ilnema-1-ima)-1
                end if
!
                do i = 1, nbnm
                    ino = zi(iamail-1+i)
                    if (ino .eq. 0) then
                        call utmess('F', 'ASSEMBLA_36')
                    end if
                    iino = ino
                    if (ino .lt. 0) iino = nbnom+nbnot-ino
                    if (nbnm .eq. 1) zi(iaexi1+iino) = zi(iaexi1+iino)+1
                    do j = i+1, nbnm
                        jno = zi(iamail-1+j)
                        jjno = jno
                        if (jno .lt. 0) jjno = nbnom+nbnot-jno
                        zi(iaexi1+iino) = zi(iaexi1+iino)+1
                        zi(iaexi1+jjno) = zi(iaexi1+jjno)+1
                    end do
                end do
            end do
        end do
        call jeveuo(nomlig//'.NBNO', 'L', ianbno)
        nbnot = nbnot+zi(ianbno)
31  end do
!
!
!
    call jedema()
end subroutine
