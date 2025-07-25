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

subroutine defsta(nmresz, numrfz, raildz, lddl, nocmp, &
                  nbfor, nbdef, tydef, inord)
    implicit none
!  P. RICHARD     DATE 09/07/91
!-----------------------------------------------------------------------
!  BUT : CALCULER LES DEFORMEES STATIQUES CORRESPONDANT A UNE FORCE
!        UNITAIRE IMPOSEE SUR LES DDL D'UNE LISTE.
!
!        REMARQUE: UNE FORCE UNITAIRE SUR UN LAGRANGE REVIENT A UN
!                  DEPLACEMENT UNITAIRE.
!
!  LA FORCE UNITAIRE PEUT ETRE IMPOSEE SOIT A UN SEUL DDL DANS LE CAS
!  D'UN DDL PHYSIQUE, SOIT A DEUX DDL DANS LE CAS DE DOUBLES DDL DE
!  LAGRANGE DE BLOCAGE, SOIT A N DDL POUR ETRE BETON.
!
!-----------------------------------------------------------------------
!
! NMRESZ /I/ : NOM UTILISATEUR DU CONCEPT RESULTAT
! NUMRFZ /I/ : NOM UTILISATEUR DU NUME_DDL DE REFERENCE
! RAILDZ /I/ : NOM UTILISATEUR DE LA MATRICE RAIDEUR FACTORISEE
! LDDL   /I/ : LISTE DES DDL A FORCER SUCCESSIVEMENT
! NOCMP  /I/ : VALEURS DES PARAMETRES "NOEUD_CMP" CORRESPONDANTES
! NBFOR  /I/ : NOMBRE DE FORCES UNITAIRES A IMPOSER PAR DEFORMEES
! NBDEF  /I/ : NOMBRE DE DEFORMEES A CALCULER
! TYDEF  /I/ : VALEURS DES PARAMETRES "TYPE_DEFO" CORRESPONDANTS
! INORD  /M/ : NUMERO D'ORDRE DE DEPART POUR STOCKAGE DANS NMRESZ
!
!
!
!
#include "jeveux.h"
#include "asterfort/convnu.h"
#include "asterfort/copcvn.h"
#include "asterfort/dismoi.h"
#include "asterfort/inivec.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/resoud.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtcrea.h"
#include "asterfort/wkvect.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: i, iadval, ier, inord, lmat
    integer(kind=8) :: ltcham, ltcvn, nbdef, nbfor, nbpabm, neq, neqr
!
!-----------------------------------------------------------------------
    parameter(nbpabm=9)
!
    integer(kind=8) :: lddl(nbfor, nbdef), ldpar(nbpabm)
    character(len=6) :: pgc
    character(len=8) :: nomres, mailla, kbid
    character(len=16) :: bmpara(nbpabm), depl, tydef, nocmp(nbdef)
    character(len=19) :: raildl, numref, numddl, matpre, solveu
    character(len=24) :: chamno, nomcvn, crefe(2), blanc
    character(len=*) :: nmresz, numrfz, raildz
    complex(kind=8) :: cbid
    integer(kind=8) :: iret
    cbid = dcmplx(0.d0, 0.d0)
!
!-----------------------------------------------------------------------
!
    data bmpara/&
     &  'NUME_MODE  ', 'FREQ', 'NORME',&
     &  'NOEUD_CMP', 'TYPE_DEFO', 'OMEGA2',&
     &  'MASS_GENE', 'RIGI_GENE', 'TYPE_MODE'/
    data pgc/'DEFSTA'/
    data blanc/'                        '/
    data depl/'DEPL'/
!
    call jemarq()
    nomres = nmresz
    numref = numrfz
    raildl = raildz
!
!-----------------------------------------------------------------------
!
    if (nbdef .eq. 0) goto 999
!
! --- RECUPERATION DU MODELE DE REFERENCE
!
    call dismoi('NOM_MAILLA', numref, 'NUME_DDL', repk=mailla)
    crefe(1) = mailla
    crefe(2) = numref
!
! --- CONVERSION DU NUMDDL ASSOCIE A LA MATRICE
!
    call dismoi('NOM_NUME_DDL', raildl, 'MATR_ASSE', repk=numddl)
    numddl(15:19) = '.NUME'
!
    call dismoi('NB_EQUA', numddl, 'NUME_DDL', repi=neq)
!
    nomcvn = '&&'//pgc//'.CONV.NUMDDL'
    numref(15:19) = '.NUME'
    call convnu(numddl, numref, nomcvn, 'V', neqr)
    call jeveuo(nomcvn, 'L', ltcvn)
!
! --- CREATION DU VECTEUR DE TRAVAIL
!
    call wkvect('&&'//pgc//'CHAMNO', 'V V R', neq, ltcham)
!
! --- BOUCLE DE CALCUL DES DEFORMEES
!
    do i = 1, nbdef
!
! ----- GENERATION DU NOM DU CHAMP RESULTAT
!
        call rsexch(' ', nomres, depl, inord, chamno, &
                    ier)
        if (ier .eq. 0) then
            call utmess('A', 'ALGORITH2_64', sk=chamno)
        else if (ier .eq. 100) then
            call vtcrea(chamno, crefe, 'G', 'R', neqr)
        else
            call utmess('F', 'ALGORITH2_65')
        end if
!
! ----- INITIALISATION DU SECOND MEMBRE
!
        call jeveuo(raildl(1:19)//'.&INT', 'E', lmat)
        call inivec(zr(ltcham), neq, lddl(1, i), nbfor)
!
! ----- RESOLUTION EN PLACE
!
        matpre = '&&OP0099.MATPRE'
        solveu = '&&OP0099.SOLVEUR'
        call resoud(raildl, matpre, solveu, ' ', 1, &
                    ' ', ' ', ' ', zr(ltcham), [cbid], &
                    ' ', .true._1, 0, iret)
!
! ----- CONVERSION NUMEROTATION
!
        chamno(20:24) = '.VALE'
        call jeveuo(chamno, 'E', iadval)
        call copcvn(neqr, zr(ltcham), zr(iadval), zi(ltcvn), 1.d0)
!
! ----- SAUVEGARDE DE LA SOLUTION
!
        call jelibe(chamno)
!
!        CHAMNO(20:24)  = '.DESC'
!
!
!        CHAMNO(20:24)  = '.REFE'
!
!
        call rsnoch(nomres, depl, inord)
!
! ----- STOCKAGE DES PSEUDO PARAMETRES
!
        call rsadpa(nomres, 'E', nbpabm, bmpara, inord, &
                    0, tjv=ldpar, styp=kbid)
        zi(ldpar(1)) = inord
        zr(ldpar(2)) = 0.d0
        zk24(ldpar(3)) = blanc
        zk16(ldpar(4)) = nocmp(i)
        zk16(ldpar(5)) = tydef
        zr(ldpar(6)) = 0.d0
        zr(ldpar(7)) = 0.d0
        zr(ldpar(8)) = 0.d0
        zk16(ldpar(9)) = 'MODE_STA'
!
! ----- INCREMENT DU NUMERO D'ORDRE
!
        inord = inord+1
!
    end do
!
    call jedetr(nomcvn)
    call jedetr('&&'//pgc//'CHAMNO')
!
999 continue
    call jedema()
end subroutine
