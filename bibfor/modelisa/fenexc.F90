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
subroutine fenexc(noma, nomnoa, long, nbn, nuno, &
                  diax, nbnfen, noefen, disfen)
    implicit none
!     EXTRACTION DE LA FENETRE EXCITEE SUR LA STRUCTURE
!     APPELANT : SPECFF
!-----------------------------------------------------------------------
! IN  : NOMA   : CHARACTER*8 , NOM DU CONCEPT MAILLAGE
! IN  : NOMNOA : CHARACTER*8 , NOM DU NOEUD D'APPLICATION (CENTRE DE
!                LA FENETRE EXCITEE)
! IN  : LONG   : REAL*8 , LONGUEUR EXCITEE
! IN  : NBN    : INTEGER , NOMBRE DE NOEUDS DU MAILLAGE
! IN  : NUNO   : INTEGER , VECTEUR DE DIM = NBRE DE NOEUDS DU MAILLAGE
!                LISTE DES NUMEROS DES NOEUDS DU MAILLAGE, REORDONNEE
!                PAR VALEURS CROISSANTES DE LA COORDONNEE D'ESPACE LE
!                LONG DE L'AXE DIRECTEUR DE LA STRUCTURE
! IN  : DIAX   : REAL*8 , VECTEUR DE DIM = NBRE DE NOEUDS DU MAILLAGE
!                LISTE CROISSANTE DES VALEURS DE LA COORDONNEE D'ESPACE
!                LE LONG DE L'AXE DIRECTEUR DE LA STRUCTURE
! OUT : NBNFEN : INTEGER , NOMBRE DE NOEUDS DU MAILLAGE APPARTENANT
!                A LA FENETRE EXCITEE
! OUT : NOEFEN : INTEGER , VECTEUR DE DIM = NBRE DE NOEUDS DU MAILLAGE
!                LISTE DES NUMEROS DES NOEUDS DU MAILLAGE APPARTENANT
!                A LA FENETRE EXCITEE, ORDONNEE PAR VALEURS CROISSANTES
!                DE LA COORDONNEE D'ESPACE LE LONG DE L'AXE DIRECTEUR
!                DE LA STRUCTURE
! OUT : DISFEN : REAL*8 , VECTEUR DE DIM = NBRE DE NOEUDS DU MAILLAGE
!                DISCRETISATION SPATIALE CORRESPONDANT A LA FENETRE
!                EXCITEE, ORDONNEE PAR VALEURS CROISSANTES ET TRANSLATEE
!                POUR ETRE SUPERPOSEE A LA DISCRETISATION DES FONCTIONS
!                DE FORME ASSOCIEES A L'EXCITATION
!
!
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/lexseg.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: noma, nomnoa
    integer(kind=8) :: nbn, nuno(*), nbnfen, noefen(*)
    real(kind=8) :: long, diax(*), disfen(*)
!
    character(len=24) :: connex, typmai
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    integer(kind=8) :: idiscp, idiscs, indap, ino, inunop, inunos, isensp
    integer(kind=8) :: isenss, lp, ls, nbnop, nbnos, nbrma, nunoap
    integer(kind=8) :: nunod, nunog
    real(kind=8) :: difx1, difx2, refx1, tol, x, x1
    real(kind=8) :: x2, xap, xdecal
!-----------------------------------------------------------------------
    call jemarq()
    tol = 100.d0*r8prem()
!
! --- 1.EXTRACTION DU NUMERO DU NOEUD D'APPLICATION DE L'EXCITATION
!
    nunoap = char8_to_int(nomnoa)
    if (nunoap .eq. 0) then
        call utmess('F', 'MODELISA4_53')
    end if
!
! --- 2.DETERMINATION DU RANG DU NOEUD D'APPLICATION DANS LA LISTE
! ---   ORDONNEE DES NOEUDS DE LA STRUCTURE
!
    indap = 1
    xap = 0.d0
    do ino = 1, nbn
        if (nuno(ino) .eq. nunoap) then
            indap = ino
            xap = diax(ino)
            goto 11
        end if
    end do
11  continue
!
! --- 3.TESTS DE RECOUVREMENT PAR RAPPORT AU DOMAINE DE DEFINITION
! ---   DU MAILLAGE
!
    if (indap .eq. 1 .or. indap .eq. nbn) then
        call utmess('F', 'MODELISA4_54')
    end if
!
    x1 = xap-long/2.d0
    x2 = xap+long/2.d0
    difx1 = diax(1)-x1
    difx2 = x2-diax(nbn)
    if (difx1 .gt. tol .or. difx2 .gt. tol) then
        call utmess('F', 'MODELISA4_55')
    end if
!
!.....DEDUCTION DE XDECAL UTILE PLUS LOIN POUR TRANSLATION DE
!.....LA DISCRETISATION
!
    xdecal = 0.d0
    refx1 = dble(abs(x1))
    if (refx1 .gt. tol) xdecal = x1
!
! --- 4.CREATION D'OBJETS DE TRAVAIL - ACCES AUX OBJETS UTILES
!
    lp = indap
    call wkvect('&&FENEXC.TEMP.NNOP', 'V V I', lp, inunop)
    call wkvect('&&FENEXC.TEMP.DISP', 'V V R', lp, idiscp)
!
    ls = nbn-indap+1
    call wkvect('&&FENEXC.TEMP.NNOS', 'V V I', ls, inunos)
    call wkvect('&&FENEXC.TEMP.DISS', 'V V R', ls, idiscs)
!
    connex = noma//'.CONNEX'
    typmai = noma//'.TYPMAIL'
    call jelira(typmai, 'LONMAX', nbrma)
!
! --- 5.DETERMINATION DE L'ENSEMBLE DES NOEUDS APPARTENANT A LA
! ---   DEMI-FENETRE EXCITEE EN AMONT DU NOEUD CENTRAL D'APPLICATION
!
    zi(inunop+lp-1) = nuno(indap)
    zr(idiscp+lp-1) = diax(indap)
    ino = 0
    nbnop = 1
!
!.....REPETER
!
20  continue
    ino = ino+1
    if (ino .gt. lp-1) goto 21
    nunog = nuno(indap-ino)
    nunod = zi(inunop+lp-nbnop)
    if (nbnop .eq. 1) then
        if (lexseg(connex, typmai, nbrma, nunog, nunod)) then
            isensp = 1
            nbnop = nbnop+1
            zi(inunop+lp-nbnop) = nunog
            x = diax(indap-ino)
            zr(idiscp+lp-nbnop) = x
            if (x-x1 .lt. tol) goto 21
        else if (lexseg(connex, typmai, nbrma, nunod, nunog)) then
            isensp = -1
            nbnop = nbnop+1
            zi(inunop+lp-nbnop) = nunog
            x = diax(indap-ino)
            zr(idiscp+lp-nbnop) = x
            if (x-x1 .lt. tol) goto 21
        end if
    else if (isensp .eq. 1) then
        if (lexseg(connex, typmai, nbrma, nunog, nunod)) then
            nbnop = nbnop+1
            zi(inunop+lp-nbnop) = nunog
            x = diax(indap-ino)
            zr(idiscp+lp-nbnop) = x
            if (x-x1 .lt. tol) goto 21
        end if
    else
        if (lexseg(connex, typmai, nbrma, nunod, nunog)) then
            nbnop = nbnop+1
            zi(inunop+lp-nbnop) = nunog
            x = diax(indap-ino)
            zr(idiscp+lp-nbnop) = x
            if (x-x1 .lt. tol) goto 21
        end if
    end if
    goto 20
!
!.....SORTIE EN ERREUR FATALE LE CAS ECHEANT
!
21  continue
    if (nbnop .eq. 1) then
        call utmess('F', 'MODELISA4_56')
    end if
!
    difx1 = zr(idiscp+lp-nbnop)-x1
    if (difx1 .gt. tol) then
        call utmess('F', 'MODELISA4_57')
    end if
!
! --- 6.DETERMINATION DE L'ENSEMBLE DES NOEUDS APPARTENANT A LA
! ---   DEMI-FENETRE EXCITEE EN AVAL DU NOEUD CENTRAL D'APPLICATION
!
    zi(inunos) = nuno(indap)
    zr(idiscs) = diax(indap)
    ino = 0
    nbnos = 1
!
!.....REPETER
!
30  continue
    ino = ino+1
    if (ino .gt. ls-1) goto 31
    nunog = zi(inunos+nbnos-1)
    nunod = nuno(indap+ino)
    if (nbnos .eq. 1) then
        if (lexseg(connex, typmai, nbrma, nunog, nunod)) then
            isenss = 1
            if (isenss .ne. isensp) then
                call utmess('F', 'MODELISA4_58')
            end if
            nbnos = nbnos+1
            zi(inunos+nbnos-1) = nunod
            x = diax(indap+ino)
            zr(idiscs+nbnos-1) = x
            if (x2-x .lt. tol) goto 31
        else if (lexseg(connex, typmai, nbrma, nunod, nunog)) then
            isenss = -1
            if (isenss .ne. isensp) then
                call utmess('F', 'MODELISA4_58')
            end if
            nbnos = nbnos+1
            zi(inunos+nbnos-1) = nunod
            x = diax(indap+ino)
            zr(idiscs+nbnos-1) = x
            if (x2-x .lt. tol) goto 31
        end if
    else if (isenss .eq. 1) then
        if (lexseg(connex, typmai, nbrma, nunog, nunod)) then
            nbnos = nbnos+1
            zi(inunos+nbnos-1) = nunod
            x = diax(indap+ino)
            zr(idiscs+nbnos-1) = x
            if (x2-x .lt. tol) goto 31
        end if
    else
        if (lexseg(connex, typmai, nbrma, nunod, nunog)) then
            nbnos = nbnos+1
            zi(inunos+nbnos-1) = nunod
            x = diax(indap+ino)
            zr(idiscs+nbnos-1) = x
            if (x2-x .lt. tol) goto 31
        end if
    end if
    goto 30
!
!.....SORTIE EN ERREUR FATALE LE CAS ECHEANT
!
31  continue
    if (nbnos .eq. 1) then
        call utmess('F', 'MODELISA4_59')
    end if
!
    difx2 = x2-zr(idiscs+nbnos-1)
    if (difx2 .gt. tol) then
        call utmess('F', 'MODELISA4_60')
    end if
!
! --- 7.AFFECTATION DES RESULTATS
!
!.....NOMBRE DE NOEUDS APPARTENANT A LA FENETRE EXCITEE
!
    nbnfen = nbnop+nbnos-1
!
!.....LISTE DES NUMEROS DES NOEUDS APPARTENANT A LA FENETRE EXCITEE
!.....DISCRETISATION SPATIALE CORRESPONDANTE (TRANSLATEE)
!
    do ino = 1, nbnop
        noefen(ino) = zi(inunop+lp-nbnop+ino-1)
        disfen(ino) = zr(idiscp+lp-nbnop+ino-1)-xdecal
    end do
    do ino = 2, nbnos
        noefen(nbnop+ino-1) = zi(inunos+ino-1)
        disfen(nbnop+ino-1) = zr(idiscs+ino-1)-xdecal
    end do
!
! --- MENAGE
!
    call jedetr('&&FENEXC.TEMP.NNOP')
    call jedetr('&&FENEXC.TEMP.DISP')
    call jedetr('&&FENEXC.TEMP.NNOS')
    call jedetr('&&FENEXC.TEMP.DISS')
!
    call jedema()
end subroutine
