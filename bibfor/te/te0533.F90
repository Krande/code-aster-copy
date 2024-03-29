! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
! person_in_charge: samuel.geniaut at edf.fr
!
subroutine te0533(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elelin.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/xdocon.h"
#include "asterfort/xmcont.h"
#include "asterfort/xmfrot.h"
#include "asterfort/xmprep.h"
#include "asterfort/xmulco.h"
#include "asterfort/xteddl.h"
#include "asterfort/xteini.h"
#include "asterfort/xkamat.h"
!
    character(len=16) :: option, nomte
!
!         CALCUL DES MATRICES DE CONTACT FROTTEMENT POUR X-FEM
!                       (METHODE CONTINUE)
!
!
!  OPTION : 'RIGI_CONT'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!......................................................................
!
    integer :: i, j, ij, ifa, ipgf, isspg
    integer :: jindco, jdonco, jlsn, ipoids, ivf, idfde, jgano, igeom
    integer :: idepm, idepd, imatt, jlst, jptint, jaint, jcface, jlonch
    integer :: ivff, iadzi, iazk24, ibid, jbasec, jseuil
    integer :: ndim, nfh, ddlc, ddls, nddl, nno, nnos, nnom, nnof, ddlm
    integer :: npg, npgf, algocr, algofr, vstnc(32)
    integer :: indco, ninter, nface, cface(30, 6)
    integer :: nfe, singu, jstno, nvit, jcoheo, ncompv, jheavn
    integer :: nnol, pla(27), lact(8), nlact, nptf
    integer :: contac, nfiss, jfisno, jmate, jcohes, nbspg, nspfis
    real(kind=8) :: ffp(27), ffc(8), coefcp, coefcr, coeffp
    real(kind=8) :: mmat(216, 216), jac, mu, coeffr, rela
    real(kind=8) :: tau1(3), tau2(3)
    real(kind=8) :: nd(3), seuil, cohes(3), coheo(3)
    real(kind=8) :: rr
    real(kind=8) :: ka, fk(27, 3, 3), mu2
    integer :: jheano, ifiss, jheafa, ncomph, jta2(3)
    integer :: jtab(7), iret, ncompd, ncompp, ncompa, ncompb, ncompc, ncompn
    integer :: jbaslo
    aster_logical :: matsym, lelim
    aster_logical :: axi
    character(len=8) :: elref, elrefc, elc, fpg
    character(len=24) :: typma
!......................................................................
!
    rela = 0.d0
!
! --- INITIALISATIONS
!
    do i = 1, 8
        lact(i) = 0
    end do
    ffp(:) = 0.d0
    tau2(:) = 0.d0
    lelim = .false.
    nbspg = 0
! INTIALISATION JMATE POUR DETECTER EVENTUELLES ERREURS JEVEUX
    jmate = 1
    call elref1(elref)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                nnom, ddls, nddl, ddlm, nfiss, &
                contac)
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3)
!
!     INITIALISATION DE LA MATRICE DE TRAVAIL
    mmat(:, :) = 0.d0
!
! --- ROUTINE SPECIFIQUE P2P1, A CONSERVER
!
    call elelin(contac, elref, elrefc, ibid, ibid)
!
! --- RECUPERATION DES ENTREES / SORTIE
!
    call jevech('PGEOMER', 'L', igeom)
!     DEPMOI
    call jevech('PDEPL_M', 'L', idepm)
!     DEPDEL
    call jevech('PDEPL_P', 'L', idepd)
    call jevech('PINDCOI', 'L', jindco)
    call jevech('PDONCO', 'L', jdonco)
    call jevech('PSEUIL', 'L', jseuil)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
    call jevech('PPINTER', 'L', jptint)
    call jevech('PAINTER', 'L', jaint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PLONGCO', 'L', jlonch)
    call jevech('PBASECO', 'L', jbasec)
    if (nfh .gt. 0) then
        call jevech('PHEA_NO', 'L', jheavn)
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
    end if
    if (nfiss .gt. 1) then
        call jevech('PFISNO', 'L', jfisno)
        call jevech('PHEAVNO', 'L', jheano)
        call jevech('PHEA_FA', 'L', jheafa)
        call tecach('OOO', 'PHEA_FA', 'L', iret, nval=2, &
                    itab=jtab)
        ncomph = jtab(2)
    end if
!     NB COMPOSANTES DES MODES LOCAUX
!     ASSOCIES AUX CHAMPS DANS LE CATALOGUE
    call tecach('OOO', 'PDONCO', 'L', iret, nval=2, itab=jtab)
    ncompd = jtab(2)
    call tecach('OOO', 'PPINTER', 'L', iret, nval=2, itab=jtab)
    ncompp = jtab(2)
    call tecach('OOO', 'PAINTER', 'L', iret, nval=2, itab=jtab)
    ncompa = jtab(2)
    call tecach('OOO', 'PBASECO', 'L', iret, nval=2, itab=jtab)
    ncompb = jtab(2)
    call tecach('OOO', 'PCFACE', 'L', iret, nval=2, itab=jtab)
    ncompc = jtab(2)
    if (nfe .gt. 0) then
        call jevech('PMATERC', 'L', jmate)
        call jevech('PBASLOR', 'L', jbaslo)
        call jevech('PSTANO', 'L', jstno)
        axi = lteatt('AXIS', 'OUI')
        call xkamat(zi(jmate), ndim, axi, ka, mu2)
    else
        ka = 3.d0
        mu2 = 1.d0
        axi = .false._1
        jmate = 0
        jbaslo = 0
        jlsn = 0
        jstno = 0
    end if
!
!     STATUT POUR L'ÉLIMINATION DES DDLS DE CONTACT
    do i = 1, max(1, nfh)*nnos
        vstnc(i) = 1
    end do
!
! - Loop on cracks
!
    do ifiss = 1, nfiss
!
! --- RECUPERATION DIVERSES DONNEES CONTACT
!
        call xdocon(algocr, algofr, cface, contac, coefcp, &
                    coeffp, coefcr, coeffr, elc, fpg, &
                    ifiss, ivff, jcface, jdonco, jlonch, &
                    mu, nspfis, ncompd, ndim, nface, &
                    ninter, nnof, nomte, npgf, nptf, &
                    rela)
        if (ninter .eq. 0) goto 91
!
! --- RECUPERATION MATERIAU ET VARIABLES INTERNES COHESIF
!
        if (algocr .eq. 3) then
            call jevech('PMATERC', 'L', jmate)
            call jevech('PCOHES', 'L', jcohes)
            call tecach('OOO', 'PCOHES', 'L', iret, nval=3, &
                        itab=jta2)
!
!           CAS COLLOCATION AUX POINTS DE GAUSS
            if (contac .eq. 1 .or. contac .eq. 3) ncompv = jta2(2)
!
!           CAS CHAMP ELNO
            if (contac .eq. 2) ncompv = jta2(2)/jta2(3)
            call jevech('PCOHESO', 'E', jcoheo)
        end if
!
! --- RECUP MULTIPLICATEURS ACTIFS ET LEURS INDICES
!
        call xmulco(contac, ddls, ddlc, ddlm, jaint, ifiss, &
                    jheano, vstnc, lact, .true._1, lelim, &
                    ndim, nfh, nfiss, ninter, &
                    nlact, nno, nnol, nnom, nnos, &
                    pla, typma)
!
! --- BOUCLE SUR LES FACETTES
!
        do ifa = 1, nface
!
! --- BOUCLE SUR LES POINTS DE GAUSS DES FACETTES
!
            do ipgf = 1, npgf
!
! --- RECUPERATION DES STATUTS POUR LE POINT DE GAUSS
!
                isspg = npgf*(ifa-1)+ipgf
                indco = zi(jindco-1+nbspg+isspg)
                if (algofr .ne. 0) seuil = zr(jseuil-1+nbspg+isspg)
                if (algocr .eq. 3 .and. contac .ne. 2) then
                    do i = 1, ncompv
                        cohes(i) = zr(jcohes+ncompv*(nbspg+isspg-1)-1+i)
                    end do
                end if
!
! --- PREPARATION DU CALCUL
!
                call xmprep(cface, contac, elref, elrefc, elc, &
                            ffc, ffp, fpg, jaint, jbasec, &
                            jptint, ifa, igeom, ipgf, jac, &
                            jlst, lact, nd, ndim, ninter, &
                            nlact, nno, nnos, nptf, nvit, &
                            rr, singu, tau1, tau2, ka, mu2, &
                            jbaslo, jstno, jlsn, fk)
!
! --- CALCUL DES MATRICES DE CONTACT
!     ..............................
!
                call xmcont(algocr, coefcr, coefcp, cohes, coheo, &
                            jcohes, jcoheo, ncompv, &
                            ddlm, ddls, ffc, ffp, idepd, &
                            idepm, ifa, ifiss, jmate, indco, &
                            ipgf, jac, jheavn, ncompn, jheafa, mmat, &
                            lact, ncomph, nd, nddl, ndim, &
                            nfh, nfiss, nno, nnol, nnos, &
                            nvit, pla, rela, singu, fk, &
                            tau1, tau2)

!
! --- SI COHESIF CLASSIQUE ON ACTUALISE LA VARIABLE INTERNE
!
                if (algocr .eq. 3 .and. (contac .eq. 1 .or. contac .eq. 3)) then
                    do i = 1, ncompv
                        zr(jcoheo+ncompv*(nbspg+isspg-1)-1+i) = coheo(i)
                    end do
                end if
                if (rela .eq. 0.d0 .or. rela .eq. 1.d0 .or. rela .eq. 2.d0) then
                    call xmfrot(algofr, coeffr, coeffp, ddlm, ddls, &
                                ffc, ffp, idepd, idepm, indco, &
                                jac, lact, mmat, mu, nd, &
                                ndim, nfh, nfiss, nno, nnol, &
                                nnos, nvit, pla, seuil, &
                                singu, fk, tau1, tau2)
                end if

!
! --- FIN DE BOUCLE SUR LES POINTS DE GAUSS
            end do
!
! --- FIN DE BOUCLE SUR LES FACETTES
        end do
! --- FIN BOUCLE SUR LES FISSURES : DECALAGE D INDICES
        nbspg = nbspg+nspfis
91      continue
        jbasec = jbasec+ncompb
        jptint = jptint+ncompp
        jaint = jaint+ncompa
        jcface = jcface+ncompc
    end do
!
!-----------------------------------------------------------------------
!     COPIE DES CHAMPS DE SORTIES ET FIN
!-----------------------------------------------------------------------
!
    if (algocr .eq. 2 .or. algofr .eq. 2 .or. algocr .eq. 3 &
        .and. rela .ne. 5.d0 .and. rela .ne. 4.d0) then
! --- RECUPERATION DE LA MATRICE 'OUT' NON SYMETRIQUE
        matsym = .false.
        call jevech('PMATUNS', 'E', imatt)
        do j = 1, nddl
            do i = 1, nddl
                ij = j+nddl*(i-1)
                zr(imatt+ij-1) = mmat(i, j)
            end do
        end do
    else
! --- RECUPERATION DE LA MATRICE 'OUT' SYMETRIQUE
        matsym = .true.
        call jevech('PMATUUR', 'E', imatt)
        do j = 1, nddl
            do i = 1, j
                ij = (j-1)*j/2+i
                zr(imatt+ij-1) = mmat(i, j)
            end do
        end do
    end if
! --- SUPPRESSION DES DDLS DE DEPLACEMENT SEULEMENT POUR LES XHTC
    if (nfh .ne. 0) then
        call jevech('PSTANO', 'L', jstno)
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, zi(jstno), .false._1, matsym, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    mat=zr(imatt))
    end if
! --- SUPPRESSION DES DDLS DE CONTACT
    if (lelim) then
        call xteddl(ndim, nfh, nfe, ddls, nddl, &
                    nno, nnos, vstnc, .true._1, matsym, &
                    option, nomte, ddlm, nfiss, jfisno, &
                    mat=zr(imatt))
    end if
!
end subroutine
