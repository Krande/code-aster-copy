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

subroutine preml1(neq, n2, diag, delg, col, &
                  xadj, adjncy, parent, adress, supnd, &
                  nnz, qsize, llist, suiv, p, &
                  q, invp, perm, lgind, ddlmoy, &
                  nbsn, optnum, lgadjn, nrl, deb, &
                  vois, suit, ier, nec, prno, &
                  deeq, noeud, ddl, invpnd, permnd, &
                  spndnd, xadjd, matgen)
! person_in_charge: olivier.boiteau at edf.fr
! aslint: disable=W1504
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/ismaem.h"
#include "asterfort/amdapt.h"
#include "asterfort/amdbar.h"
#include "asterfort/genmmd.h"
#include "asterfort/infniv.h"
#include "asterfort/prmadj.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: neq, diag(0:neq), lgind, lgadjn
    integer(kind=8) :: col(*), xadj(neq+1), adjncy(lgadjn)
    integer(kind=8) :: delg(neq), nbsn, adress(neq), parent(neq)
    integer(kind=8) :: supnd(neq), optnum
    integer(kind=8) :: invp(neq), perm(neq), ier
    integer(kind=8) :: deb(*), vois(*), suit(*)
    integer(kind=8) :: noeud(*), ddl(*), permnd(*), invpnd(*), spndnd(*), xadjd(*)
!     VARIABLES LOCALES
    integer(kind=8) :: i, j, idelta, maxint, nadj
    integer(kind=8) :: n2, fctnzs, innz, numi, numj, numl, num, ii
    real(kind=8) :: nbops
    integer(kind=8) :: nnz(1:neq), qsize(neq), llist(neq), suiv(neq)
    integer(kind=8) :: libre, iovflo, ncmpa, ifm, niv, p(neq), q(n2), nrl
    integer(kind=8) :: it
    integer(kind=8) :: nec, prno(*), deeq(*), ino, nbcmp
    aster_logical :: matgen, liaiso
!--------------------------------------------------------------
!
!     VERSION RENUMEROTATION PAR NOEUD
    integer(kind=8) :: nbnd, nd, nbnd1, ddlmoy, renum
    integer(kind=8) :: pas, k, ndanc, iddl, sni, ind, nddl
    integer(kind=8) :: vali(2)
    character(len=24) :: noxadj, noadjn, noperm, noinvp, nopare, nospnd
!
    data noxadj/'&&PREML1.NOMXADJ_PROV   '/
    data noadjn/'&&PREML1.NOMADJN_PROV   '/
    data noperm/'&&PREML1.NOPERM_PROV    '/
    data noinvp/'&&PREML1.NOINVP_PROV    '/
    data nopare/'&&PREML1.NOPARE_PROV    '/
    data nospnd/'&&PREML1.NOSPND_PROV    '/
!****************************************************************
!-----RECUPERATION DU NIVEAU D'IMPRESSION
!
!
    call infniv(ifm, niv)
!--------------------------------------------------------------
!     LA RENUMEROTATION VA ETRE FAITE AVEC LA CONNECTIVITE NODALE
!     ET NON PLUS LA CONNECTIVITE PAR DDL COMME AVANT
!
!     1) ON CALCULE   NOEUD ET DDL A PARTIR DE PRNO ET DEEQ
!     NOEUD(1:NEQ) SURJECTION VERS (1:NBND) OU NBND EST LE NOMBRE
!     DE NOEUDS AU SENS DE LA DISCRETISATION
!     LA FONCTION INVERSE EST DDL(1:NBND)
!     EN CAS DE MATRICE GENERALISEE PRNO ET DEEQ N ONT PAS LE MEME SENS
!     ON NE LES UTILISE PAS MAIS ON SAIT QUE NOEUDS ET DDL SONT
!     IDENTIQUES, LES TABLEAUX DDL ET NOEUD SONT ALORS IMMEDIATS
!     SI ON N'AVAIT PLUS UN DDL PAR NOEUD EN MATRICE GENERALISEE
!     CECI NE SERAIT PLUS VALABLE
    renum = optnum
    if (matgen) then
        nbnd = neq
!        DANS LE CAS MATGEN, ON ASSIMILE LES NOEUDS AUX DDL
!        ON SUPPOSE AUSSI QUE NEQ == N2 IE PAS DE LAGRANGE
!        SI ON AVAIT N2< NEQ ( CERTAINS DELG NON NULS)
!        IL FAUDRAIT P.E DEPLACER LA BOUCLE 3 APRES LA 350
!        AVEC N2 EN PLACE DE NEQ
        do i = 1, neq+1
            ddl(i) = i
            noeud(i) = i
        end do
!     RECHERCHE DE PRESENCE DE LIAISON NON DETECTEE PAR PREML0(NRL=0)
!     (NON SIGNALEE PAR DES VALEURS DE DELG NEGATIVES)
        liaiso = .false.
        do i = 1, neq
            if (deeq(2*i) .lt. 0) then
                liaiso = .true.
            end if
        end do
        if (liaiso) then
!          EN CAS DE LIAISON DANS UNE MATRICE GENERALISEE
!          LES DDL SONT CORRECTEMENT ORDONNES A PRIORI
!          ON NE FAIT PAS DE RENUMEROTATION
            call utmess('A', 'ALGELINE3_35')
            renum = 3
        end if
    else
        do i = 1, neq
            noeud(i) = 0
        end do
        nbnd = 0
        i = 1
        ddl(i) = 1
101     continue
!     DO WHILE(I.LE.NEQ)
        if (i .le. neq) then
            if (deeq(2*i) .gt. 0) then
                ino = deeq(2*i-1)
                nbnd = nbnd+1
                nbcmp = prno((ino-1)*(2+nec)+2)
                do j = i, i+nbcmp-1
                    noeud(j) = nbnd
                end do
                ddl(nbnd+1) = ddl(nbnd)+nbcmp
                pas = nbcmp
            else
                pas = 1
            end if
            i = i+pas
            goto 101
        end if
    end if
    ddlmoy = (ddl(nbnd+1)-1)/nbnd
!--------------------------------------------------------------------
!     2) CALCUL DE (ADJNCY, XADJ) EN DDL DANS LA NUMEROTATION DE 1 À N2
!     COMME DANS LA VERSION INITIALE
!     INITIALISATION DE NNZ : NBRE DE TERMES A AJOUTER
!     POUR CHAQUE LIGNE
    num = 0
    do i = 1, neq
        if (delg(i) .eq. 0) then
            num = num+1
!     ON COMPTE LES NON-ZEROS
            innz = 0
            do j = diag(i-1)+1, diag(i)-1
                if (delg(col(j)) .eq. 0) then
!     PARTIE TRIANGULAIRE INFERIEURE
                    innz = innz+1
!     PARTIE TRIANGULAIRE SUPERIEURE
                    nnz(p(col(j))) = nnz(p(col(j)))+1
                end if
            end do
            nnz(num) = innz
        end if
    end do
    if (nrl .ne. 0) then
        do j = 1, neq
            it = deb(j)
219         continue
            if (it .gt. 0) then
                nnz(p(j)) = nnz(p(j))+1
                it = suit(it)
                goto 219
            end if
        end do
!     VERIFICATION
!
        do j = 1, neq
!     TERMES A AJOUTER PARTIE SUPERIEURE
            it = deb(j)
324         continue
            if (it .gt. 0) then
                nnz(p(vois(it))) = nnz(p(vois(it)))+1
                it = suit(it)
                goto 324
            end if
        end do
    end if
!
    xadj(1) = 1
    do j = 1, n2
        xadj(j+1) = xadj(j)+nnz(j)
        nnz(j) = 0
    end do
    if ((xadj(neq+1)-1) .gt. lgadjn) then
!     TEST D'ESPACE SUFFISANT DANS ADJNCY
        vali(1) = lgadjn
        vali(2) = xadj(neq+1)-1
        call utmess('F', 'ALGELINE4_4', ni=2, vali=vali)
    end if
!
    do j = 1, neq
        if (delg(j) .eq. 0) then
            numj = p(j)
            do ii = diag(j-1)+1, diag(j)-1
                i = col(ii)
                if (delg(i) .eq. 0) then
                    numi = p(i)
                    adjncy(xadj(numj)+nnz(numj)) = numi
                    nnz(numj) = nnz(numj)+1
                    adjncy(xadj(numi)+nnz(numi)) = numj
                    nnz(numi) = nnz(numi)+1
                end if
            end do
            if (nrl .ne. 0) then
                it = deb(j)
344             continue
                if (it .gt. 0) then
                    numl = p(vois(it))
                    adjncy(xadj(numj)+nnz(numj)) = numl
                    nnz(numj) = nnz(numj)+1
                    adjncy(xadj(numl)+nnz(numl)) = numj
                    nnz(numl) = nnz(numl)+1
                    it = suit(it)
                    goto 344
                end if
            end if
        end if
    end do
    nbnd1 = nbnd+1
    libre = xadj(nbnd1)
    nadj = libre-1
!-----------------------------------------------------------
!     3) MODIFICATION DE (ADJNCY, XADJ) VOISINAGE PAR DDL EN
!     (ADJNCY,XADJD) VOISINAGE PAR  NOEUD POUR LA RENUMEROTATION
    call prmadj(nbnd, neq, n2, adjncy, xadj, &
                xadjd, llist, q, noeud)
!-----------------------------------------------------------
    nbnd1 = nbnd+1
    libre = xadjd(nbnd1)
    nadj = libre-1
    if (renum .eq. 0) then
!----------------------------------MINIMUM DEGRE : GENMMD
        idelta = 0
        maxint = 2*nbnd
!
        call genmmd(nbnd, nbnd1, nadj, xadjd, adjncy, &
                    maxint, idelta, invpnd, permnd, nbsn, &
                    spndnd, adress, parent, lgind, fctnzs, &
                    nbops, nnz, qsize, llist, suiv)
    else if (renum .eq. 1) then
!----------------------------------MINIMUM DEGRE : APPROXIMATE MIN DEG
        iovflo = ismaem()
        do i = 1, n2
            qsize(i) = xadjd(i+1)-xadjd(i)
        end do
        call amdbar(nbnd, xadjd, adjncy, qsize, lgadjn, &
                    libre, suiv, llist, permnd, nnz, &
                    invpnd, parent, ncmpa, adress, iovflo)
        call amdapt(neq, nbnd, nbsn, xadjd, suiv, &
                    invpnd, parent, spndnd, adress, lgind, &
                    fctnzs, nbops, llist, qsize)
    else if (renum .eq. 3) then
!-----------------MATRICE GENERALISEE PAS DE RENUMEROTATION
!     ON L'EMULE EN CREANT UN SEUL SUPER NOEUD
!
        nbsn = 1
        do i = 1, nbnd
            invpnd(i) = i
            permnd(i) = i
        end do
        parent(nbsn) = 0
        spndnd(1) = 1
        spndnd(nbsn+1) = nbnd+1
        lgind = 10
        nbops = 0.d0
        fctnzs = nbnd*(nbnd+1)/2
    end if
!****************************************************************
!****************************************************************
!.....................................................................
    lgind = lgind*ddlmoy
!     4) TRAITEMENT POUR LE PASSAGE DE LA RENUMEROTATION PAR NOEUD
!     A CELLE PAR DDL
!     PERM, INVP, SUPND SONT RECONSTITUES
    ind = 0
    do nd = 1, nbnd
        ndanc = permnd(nd)
        nddl = ddl(ndanc+1)-ddl(ndanc)
        do k = 1, nddl
            ind = ind+1
            perm(ind) = ddl(ndanc)+k-1
        end do
    end do
    if (niv .eq. 2 .and. ind .ne. n2) then
        vali(1) = n2
        vali(2) = ind
        call utmess('F', 'ALGELINE4_60', ni=2, vali=vali)
    end if
    do iddl = 1, n2
        invp(perm(iddl)) = iddl
    end do
!     SUPND
    supnd(nbsn+1) = n2+1
    do sni = 1, nbsn
!        ND 1ER NOEUD DU SNI
        nd = spndnd(sni)
!        ND : ANCIEN NUMERO DE ND
        nd = permnd(nd)
!        IDDL : 1ER DDL DE CE NOEUD
        iddl = ddl(nd)
        supnd(sni) = invp(iddl)
!        INVP(IDDL) 1ER DDL DU SNI
    end do
!.....................................................................
    if (niv .ge. 1) then
!     WRITE(IFM,*)'RENUMEROTATION PAR MINIMUM DEGRE'//
!     *  'TEMPS CPU',TEMPS(3),
!     *  ' + TEMPS CPU SYSTEME ',TEMPS(6)
    end if
    fctnzs = fctnzs+neq
!
    if (niv .ge. 2) then
        write (ifm, *) '--- RESULTATS DE LA RENUMEROTATION : '
        write (ifm, *) '   --- NOMBRE DE NOEUDS ', nbnd
        write (ifm, *) '   --- LONGUEUR DE LA MATRICE INITIALE ', diag(neq)
        write (ifm, *) '   --- NOMBRE DE SUPERNOEUDS ', nbsn
        if (renum .eq. 2) then
            write (ifm, *) '   --- NOMBRE D''OP. FLOTTANTES ', nbops
        end if
    end if
    ier = 0
end subroutine
