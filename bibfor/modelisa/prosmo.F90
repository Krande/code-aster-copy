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

subroutine prosmo(matrez, limat, nbmat, basez, numedd, &
                  lsym, rouc)
!.======================================================================
    implicit none
!
!     PROSMO  --  LE BUT DE CETTE ROUTINE EST DE CONSTRUIRE LA MATR_ASSE
!                 DE NOM MATRES QUI VA RESULTER DE LA COMBINAISON
!                 LINEAIRE DES NBMAT MATR_ASSE DE LA LISTE LISMAT
!                 DE NOMS DE MATR_ASSE. LES MATRICES ONT UN STOCKAGE
!                 MORSE
!
!
!   ARGUMENT        E/S  TYPE         ROLE
!    MATREZ         OUT    K*     NOM DE LA MATR_ASSE RESULTANT DE LA
!                                 COMBINAISON LINEAIRE DES MATR_ASSE
!                                 DE LA LISTE LISMAT.
!    LIMAT          IN    K24     LISTE DES MATR_ASSE A COMBINER
!                                 DES MATR_ASSE A COMBINER.
!    NBMAT          IN    I       ON FAIT LA COMBINAISON LINEAIRE
!                                 DES NBMAT PREMIERS MATR_ASSE DE LA
!                                 LISTE LIMAT.
!    BASEZ          IN    K*      NOM DE LA BASE SUR LAQUELLE ON
!                                 CONSTRUIT LA MATR_ASSE.
!    NUMEDD         IN    K14    NOM DU NUME_DDL SUR LEQUEL S'APPUIERA
!                                 LA MATR_ASSE MATREZ
!        SI NUMEDD  =' ', LE NOM DU NUME_DDL SERA OBTENU PAR GCNCON
!        SI NUMEDD /=' ', ON PRENDRA NUMEDD COMME NOM DE NUME_DDL
!
!    LSYM           IN    L      /.TRUE.  : MATRICE SYMETRIQUE
!                                /.FALSE. : MATRICE NON-SYMETRIQUE
!    ROUC           IN    K1     /'R ' : MATRICE REELLE
!                                /'C'  : MATRICE COMPLEXE
!
!.========================= DEBUT DES DECLARATIONS ====================
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcncon.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jevtbl.h"
#include "asterfort/jexnum.h"
#include "asterfort/uttrii.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    real(kind=8) :: tmax
! -----  ARGUMENTS
    integer(kind=8) :: nbmat
    aster_logical :: lsym
    character(len=*) :: matrez, basez, numedd
    character(len=*) :: limat(nbmat)
    character(len=1) :: rouc
! -----  VARIABLES LOCALES
    character(len=1) :: base
    character(len=14) :: numddl, numdd1, numddi
    character(len=19) :: matres, mat1, mati, nume_equa1, nume_equal
    character(len=24) :: ksmhc, ksmdi, krefa, kconl, kvalm
    character(len=24) :: krefi, kliste
    integer(kind=8) :: lgbl, jhtc, i, iadi, jeq, nbter, ibl1, lcumu, kbl
    integer(kind=8) :: jbl1
    integer(kind=8) :: iblav, idhcoi, icum, ismdi, lsmhc, nterm, idsmhc, l, jsmde
    integer(kind=8) :: itbloc, nbloc, kbloc, jrefa, idrefi, idconl, ieq
    integer(kind=8) :: ier, neq, k, htc
    integer(kind=8), pointer :: ibl(:) => null()
    integer(kind=8), pointer :: pbl(:) => null()
    integer(kind=8), pointer :: smde(:) => null()
!
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    call jemarq()
!
! --- INITIALISATIONS :
!     ---------------
    base = basez
    matres = matrez
!
!
! --- NOM DU NUME_DDL A CONSTRUIRE :
!     ----------------------------
    if (numedd .eq. ' ') then
        call gcncon('_', numddl(1:8))
        numddl(9:14) = '.NUDDL'
!
    else
        numddl = numedd
    end if
!
! --- NOM DE LA PREMIERE MATR_ASSE :
!     ----------------------------
    mat1 = limat(1)
!
! --- RECUPERATION DU NUME_DDL ATTACHE A LA PREMIERE MATR_ASSE  :
!     --------------------------------------------------------
    call dismoi('NOM_NUME_DDL', mat1, 'MATR_ASSE', repk=numdd1)
!
!
! --- RECOPIE DU NUME_EQUA DE LA PREMIERE MATRICE SUR LA MATRICE
! --- RESULTANTE :
!     ---------
    nume_equa1 = numdd1//'.NUME'
    nume_equal = numddl//'.NUME'
    call copisd('NUME_EQUA', base, nume_equa1, nume_equal)
!
! --- RECUPERATION DU NOMBRE D'EQUATIONS DE LA PREMIERE MATRICE
! --- A COMBINER (C'EST LE MEME POUR TOUTES LES MATRICES) :
!     ---------------------------------------------------
    call jeveuo(numdd1//'.SMOS.SMDE', 'L', vi=smde)
    neq = smde(1)
!
!
!     7) CONSTRUCTION DE L'OBJET KLISTE QUI CONTIENDRA LES DIFFERENTS
!        SMHC(MAT_I) MIS BOUT A BOUT (EQUATION PAR EQUATION) :
!        KLISTE(JEQ)=SMHC(IMAT_1)(JEQ)//SMHC(IMAT_2)(JEQ)//...
!        KLISTE EST ALLOUé PAR "BLOC" POUR EVITER D'UTILISER
!        TROP DE MEMOIRE
!     =========================================================
    kliste = '&&PROSMO.KLISTE'
!     RECUPERATION DE LA TAILLE DES BLOCS DONNEE DANS LA COMMANDE DEBUT:
    tmax = jevtbl('TAILLE_BLOC')
    lgbl = int(tmax*1024)
!
!     7-1) HTC : HAUTEUR CUMULEE DE KLISTE(JEQ)
!     --------------------------------------------------------
    call wkvect('&&PROSMO.HTC', 'V V I', neq, jhtc)
    do i = 1, nbmat
        mati = limat(i)
        call dismoi('NOM_NUME_DDL', mati, 'MATR_ASSE', repk=numddi)
        call jeveuo(numddi//'.SMOS.SMDI', 'L', iadi)
        do jeq = 1, neq
            if (jeq .eq. 1) then
                nbter = 1
!
            else
                nbter = zi(iadi-1+jeq)-zi(iadi-1+jeq-1)
            end if
            zi(jhtc-1+jeq) = zi(jhtc-1+jeq)+nbter
        end do
        call jelibe(numddi//'.SMOS.SMDI')
    end do
!
!     7-2) IBL : NUMERO DU BLOC DE KLISTE(JEQ) :
!          PBL : POSITION DE L'EQUATION JEQ DANS LE BLOC IBL  :
!     ------------------------------------------------------------
    AS_ALLOCATE(vi=ibl, size=neq)
    AS_ALLOCATE(vi=pbl, size=neq)
    ibl1 = 1
    lcumu = 0
    do jeq = 1, neq
        htc = zi(jhtc-1+jeq)
        ASSERT(htc .le. lgbl)
!       -- SI ON CHANGE DE BLOC :
        if (lcumu+htc .gt. lgbl) then
            ibl1 = ibl1+1
            lcumu = 0
        end if
        ibl(jeq) = ibl1
        pbl(jeq) = lcumu
        lcumu = lcumu+htc
    end do
!
!     7-3) ALLOCATION DE KLISTE :
!     ------------------------------------------------------------
    if (ibl1 .eq. 1) lgbl = lcumu
    call jecrec(kliste, 'V V I', 'NU', 'DISPERSE', 'CONSTANT', &
                ibl1)
    call jeecra(kliste, 'LONMAX', lgbl)
    do kbl = 1, ibl1
        call jeveuo(jexnum(kliste, kbl), 'E', jbl1)
        call jelibe(jexnum(kliste, kbl))
    end do
!
!     7-4) REMPLISSAGE DE KLISTE :
!     ------------------------------------------------------------
    call jedetr('&&PROSMO.HTC')
    call wkvect('&&PROSMO.HTC', 'V V I', neq, jhtc)
    iblav = 1
    call jeveuo(jexnum(kliste, iblav), 'E', jbl1)
    do i = 1, nbmat
        mati = limat(i)
        call dismoi('NOM_NUME_DDL', mati, 'MATR_ASSE', repk=numddi)
        call jeveuo(numddi//'.SMOS.SMDI', 'L', iadi)
        call jeveuo(numddi//'.SMOS.SMHC', 'L', idhcoi)
        icum = 0
        do jeq = 1, neq
            if (jeq .eq. 1) then
                nbter = 1
!
            else
                nbter = zi(iadi-1+jeq)-zi(iadi-1+jeq-1)
            end if
!
!         LE BLOC CONTENANT J DOIT-IL ETRE RAMENE EN MEMOIRE ?
            ibl1 = ibl(jeq)
            if (iblav .ne. ibl1) then
                call jelibe(jexnum(kliste, iblav))
                call jeveuo(jexnum(kliste, ibl1), 'E', jbl1)
                iblav = ibl1
            end if
            do k = 1, nbter
                zi(jbl1+pbl(jeq)+zi(jhtc-1+jeq)+k-1) = zi4(idhcoi+ &
                                                           icum+(k-1))
            end do
            icum = icum+nbter
            zi(jhtc-1+jeq) = zi(jhtc-1+jeq)+nbter
        end do
        call jelibe(numddi//'.SMOS.SMDI')
        call jelibe(numddi//'.SMOS.SMHC')
    end do
    call jelibe(jexnum(kliste, iblav))
!
!
!     7-5) COMPACTAGE DE L'OBJET KLISTE
!     8)   ET CREATION  DU TABLEAU .SMDI
!     ===================================
    ksmdi = numddl//'.SMOS.SMDI'
    call wkvect(ksmdi, base//' V I', neq, ismdi)
!
    lsmhc = 0
    iblav = 1
    call jeveuo(jexnum(kliste, iblav), 'E', jbl1)
    do jeq = 1, neq
!       LE BLOC CONTENANT JEQ DOIT-IL ETRE RAMENE EN MEMOIRE ?
        ibl1 = ibl(jeq)
        if (iblav .ne. ibl1) then
            call jelibe(jexnum(kliste, iblav))
            call jeveuo(jexnum(kliste, ibl1), 'E', jbl1)
            iblav = ibl1
        end if
!
!       ON TRIE ET ORDONNE LA COLONNE (EN PLACE)
        nterm = zi(jhtc-1+jeq)
        call uttrii(zi(jbl1+pbl(jeq)), nterm)
        zi(jhtc-1+jeq) = nterm
        if (jeq .eq. 1) then
            ASSERT(nterm .eq. 1)
            zi(ismdi+1-1) = nterm
!
        else
            zi(ismdi+jeq-1) = zi(ismdi+(jeq-1)-1)+nterm
        end if
        lsmhc = lsmhc+nterm
    end do
    call jelibe(jexnum(kliste, iblav))
!
!
!     9) CREATION ET AFFECTATION DU TABLEAU .SMHC
!     ====================================================
    ksmhc = numddl//'.SMOS.SMHC'
    call wkvect(ksmhc, base//' V S', lsmhc, idsmhc)
    iblav = 1
    call jeveuo(jexnum(kliste, iblav), 'E', jbl1)
    l = 0
    do jeq = 1, neq
!       LE BLOC CONTENANT JEQ DOIT-IL ETRE RAMENE EN MEMOIRE ?
        ibl1 = ibl(jeq)
        if (iblav .ne. ibl1) then
            call jelibe(jexnum(kliste, iblav))
            call jeveuo(jexnum(kliste, ibl1), 'E', jbl1)
            iblav = ibl1
        end if
!
        nterm = zi(jhtc-1+jeq)
        do k = 1, nterm
            l = l+1
            zi4(idsmhc-1+l) = zi(jbl1+pbl(jeq)-1+k)
        end do
    end do
    call jelibe(jexnum(kliste, iblav))
    call jedetr(kliste)
    call jedetr('&&PROSMO.HTC')
    AS_DEALLOCATE(vi=ibl)
    AS_DEALLOCATE(vi=pbl)
!
!
!     10) CREATION ET AFFECTATION DU TABLEAU .IABL
!     ========================================================
!
!
!     11) CREATION ET AFFECTATION DU TABLEAU .SMDE
!     =============================================
!
    call wkvect(numddl//'.SMOS.SMDE', base//' V I', 6, jsmde)
!
! --- RECUPERATION DE LA TAILLE DU BLOC DE LA MATRICE RESULTANTE
! --- (NOMBRE DE TERMES NON NULS DE LA MATRICE)
    itbloc = zi(ismdi+neq-1)
!
    zi(jsmde-1+1) = neq
    zi(jsmde-1+2) = itbloc
    zi(jsmde-1+3) = 1
!
!
!     12) CREATION ET AFFECTATION DE LA COLLECTION .VALM
!     ===================================================
    kvalm = matres//'.VALM'
    call jedetr(kvalm)
    if (lsym) then
        nbloc = 1
!
    else
        nbloc = 2
    end if
    call jecrec(kvalm, base//' V '//rouc, 'NU', 'DISPERSE', 'CONSTANT', &
                nbloc)
    call jeecra(kvalm, 'LONMAX', itbloc)
    do kbloc = 1, nbloc
        call jecroc(jexnum(kvalm, kbloc))
    end do
!
!
!     13) CREATION ET AFFECTATION DU TABLEAU .REFA
!     ============================================================
    krefa = matres//'.REFA'
    kconl = matres//'.CONL'
    call jeexin(krefa, ier)
    if (ier .eq. 0) then
        call wkvect(krefa, base//' V K24', 20, jrefa)
!
    else
        call jeveuo(krefa, 'E', jrefa)
    end if
    zk24(jrefa-1+2) = numddl
    if (lsym) then
        zk24(jrefa-1+9) = 'MS'
!
    else
        zk24(jrefa-1+9) = 'MR'
    end if
    zk24(jrefa-1+10) = 'NOEU'
!
    do i = 1, nbmat
        mati = limat(i)
        krefi = mati//'.REFA'
        call jeveuo(krefi, 'L', idrefi)
        if (zk24(idrefi+1-1) .ne. ' ') then
            zk24(jrefa-1+1) = zk24(idrefi+1-1)
            goto 130
!
        end if
    end do
130 continue
!
!
!     15) CREATION ET AFFECTATION DU VECTEUR .CONL
!     =============================================
    call wkvect(kconl, base//' V R', neq, idconl)
    do ieq = 1, neq
        zr(idconl+ieq-1) = 1.d0
    end do
!
    call jedema()
!
end subroutine
