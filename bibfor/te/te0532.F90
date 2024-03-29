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
!
subroutine te0532(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elelin.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/xdocon.h"
#include "asterfort/xkamat.h"
#include "asterfort/xmmsa2.h"
#include "asterfort/xmmsa3.h"
#include "asterfort/xmmsa5.h"
#include "asterfort/xmmsa6.h"
#include "asterfort/xmprep.h"
#include "asterfort/xmulco.h"
#include "asterfort/xteini.h"
#include "asterfort/xxlag2.h"
#include "asterfort/xxlagm.h"
#include "asterfort/xxlan5.h"
#include "asterfort/xxmmvd.h"
!
    character(len=16) :: option, nomte
!
! person_in_charge: samuel.geniaut at edf.fr
!.......................................................................
!
!                CONTACT X-FEM MÉTHODE CONTINUE :
!         CONVERGENCE DE LA BOUCLE SUR LES CONTRAINTES ACTIVES
!
!
!  OPTION : 'XCVBCA' (X-FEM CONVERGENCE BOUCLE CONTRAINTES ACTIVES )
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!......................................................................
!
    integer :: algocr, i, j, ifa, ipgf, isspg, ino
    integer :: jindco, jdonco, jlst, ipoids, ivf, idfde, jgano, igeom
    integer :: idepl, jptint, jaint, jcface, jlonch, jgliss
    integer :: ivff, iadzi, iazk24, ibid, jout1, jout2
    integer :: jout3, jmemco, ndim, nfh, ddlc, ddls, ddlm
    integer :: npg, npgf, incoca, nfe, ninter, nnof, vstnc(1)
    integer :: indco, gliss, memco, nface, cface(30, 6)
    integer :: nno, nnos, nnom, nnol, pla(27), lact(8), nlact, nvec
    integer :: contac, jbasec, nddl, nfiss, jfisno
    integer :: jmate, singu, jcohes, jcoheo, jheano, ifiss, jheafa, jheavn, ncomph
    integer :: jtab(7), iret, ncompd, ncompp, ncompa, ncompb, ncompc, ncompn
    integer :: nbspg, nspfis, nvit, ncompv, jta2(3)
    integer :: nptf
    character(len=8) :: elref, elrefc, job, elc, fpg, champ
    character(len=24) :: typma
    real(kind=8) :: ffpc(27), rela, eps, rhon
    real(kind=8) :: reac, ffp(27), ffc(8), r3bid(3), jac
    real(kind=8) :: prec, nd(3), dn, saut(3), rr, rbid, r3bd(3)
    real(kind=8) :: coefcp, coeffp, coefcr, coeffr, r6bid(6), wsaut(3)
    real(kind=8) :: p(3, 3), pp(3, 3), dsidep(6, 6), tau1(3), lamb(3)
    real(kind=8) :: tau2(3), alpha(3), dnor(3), dtang(3), am3(3), sigma(6)
    real(kind=8) :: cohes(3), mat3bd(3, 3), mat6bd(6, 6)
    real(kind=8) :: ka, mu, fk(27, 3, 3)
    integer :: jbaslo, jlsn, jstno
    parameter(prec=1.d-16)
    aster_logical :: imprim, lbid
    aster_logical :: axi
    integer :: zxain
!......................................................................
!
!
    imprim = .false.
    incoca = 0
    nbspg = 0
    zxain = xxmmvd('ZXAIN')
    ffpc(:) = 0.d0
    ffp(:) = 0.d0
    ffc(:) = 0.d0
    do i = 1, 8
        lact(i) = 0
    end do
!
    call elref1(elref)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                nnom, ddls, nddl, ddlm, nfiss, &
                contac)
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3)
!
! --- ROUTINE SPECIFIQUE P2P1
!.
    call elelin(contac, elref, elrefc, ibid, ibid)
!
!     RECUPERATION DES ENTRÉES / SORTIE
    call jevech('PGEOMER', 'L', igeom)
!     DEPLACEMENT TOTAL COURANT (DEPPLU) : 'PDEPL_P'
    call jevech('PDEPL_P', 'L', idepl)
    call jevech('PINDCOI', 'L', jindco)
    call jevech('PDONCO', 'L', jdonco)
    call tecach('OOO', 'PDONCO', 'L', ibid, nval=2, &
                itab=jtab)
    ncompd = jtab(2)
    call jevech('PGLISS', 'L', jgliss)
    call jevech('PLST', 'L', jlst)
    call jevech('PPINTER', 'L', jptint)
    call jevech('PAINTER', 'L', jaint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PLONGCO', 'L', jlonch)
    call jevech('PMEMCON', 'L', jmemco)
    call jevech('PBASECO', 'L', jbasec)
    if (nfh .gt. 0 .and. option .ne. 'XCVBCA_MORTAR') then
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
!     DIMENSSION DES GRANDEURS DANS LA CARTE
    call tecach('OOO', 'PDONCO', 'L', iret, nval=2, &
                itab=jtab)
    ncompd = jtab(2)
    call tecach('OOO', 'PPINTER', 'L', iret, nval=2, &
                itab=jtab)
    ncompp = jtab(2)
    call tecach('OOO', 'PAINTER', 'L', iret, nval=2, &
                itab=jtab)
    ncompa = jtab(2)
    call tecach('OOO', 'PBASECO', 'L', iret, nval=2, &
                itab=jtab)
    ncompb = jtab(2)
    call tecach('OOO', 'PCFACE', 'L', iret, nval=2, &
                itab=jtab)
    ncompc = jtab(2)
    if (nfe .gt. 0) then
        call jevech('PMATERC', 'L', jmate)
        call jevech('PBASLOR', 'L', jbaslo)
        call jevech('PLSN', 'L', jlsn)
        call jevech('PSTANO', 'L', jstno)
        axi = lteatt('AXIS', 'OUI')
        call xkamat(zi(jmate), ndim, axi, ka, mu)
    else
        ka = 3.d0
        mu = 1.d0
        axi = .false._1
        jmate = 0
        jbaslo = 0
        jlsn = 0
        jstno = 0
    end if
!
    call jevech('PINCOCA', 'E', jout1)
    call jevech('PINDCOO', 'E', jout2)
    call jevech('PINDMEM', 'E', jout3)
!
! --- BOUCLE SUR LES FISSURES
!
    do ifiss = 1, nfiss
!
! --- RECUPERATION DIVERSES DONNEES CONTACT
!
        call xdocon(algocr, ibid, cface, contac, coefcp, &
                    coeffp, coefcr, coeffr, elc, fpg, &
                    ifiss, ivff, jcface, jdonco, jlonch, &
                    rbid, nspfis, ncompd, ndim, nface, &
                    ninter, nnof, nomte, npgf, nptf, &
                    rela)
        if (ninter .eq. 0 .and. (contac .eq. 1 .or. contac .eq. 3)) then
            jbasec = jbasec+ncompb
            jptint = jptint+ncompp
            jaint = jaint+ncompa
            jcface = jcface+ncompc
        end if
        if (ninter .eq. 0) goto 90
!
! --- RECUPERATION MATERIAU ET VARIABLES INTERNES COHESIF
!
        if (algocr .eq. 3) then
            call jevech('PMATERC', 'L', jmate)
            call jevech('PCOHES', 'L', jcohes)
            call jevech('PCOHESO', 'E', jcoheo)
            call tecach('OOO', 'PCOHES', 'L', iret, nval=3, &
                        itab=jta2)
            if (contac .eq. 2) ncompv = jta2(2)/jta2(3)
            if (contac .eq. 1 .or. contac .eq. 3) ncompv = jta2(2)
        end if
!
!       IMPRESSION (1ERE PARTIE)
!
        if (imprim) then
            write (6, *) ' '
            write (6, 697)
697         format(4x, 'I_IN', 7x, 'DN', 11x, 'REAC', 7x, 'I_OUT')
        end if
!
! --- RECUP MULTIPLICATEURS ACTIFS ET LEURS INDICES
!
        call xmulco(contac, ddls, ddlc, ddlm, jaint, &
                    ifiss, jheano, vstnc, lact, .false._1, &
                    lbid, ndim, nfh, nfiss, ninter, &
                    nlact, nno, nnol, nnom, nnos, &
                    pla, typma)
!
!       SI CONTACT "MORTAR"
!
        if (contac .eq. 2) then
!
            call xmprep(cface, contac, elref, elrefc, elc, &
                        ffc, ffp, fpg, jaint, jbasec, &
                        jptint, 1, igeom, 1, jac, &
                        jlst, lact, nd, ndim, ninter, &
                        nlact, nno, nnos, nptf, nvit, &
                        rr, singu, tau1, tau2)
            do ino = 1, nnol
                do i = 1, ncompv
                    cohes(i) = zr(jcohes+ncompv*(ino-1)-1+i)
                end do
                nvec = 1
                champ = 'W'
                call xxlan5(ino, idepl, ibid, ibid, lact, &
                            ndim, pla, wsaut, nvec, champ)
                nvec = 1
                champ = 'LAMBDA'
                call xxlan5(ino, idepl, ibid, ibid, lact, &
                            ndim, pla, lamb, nvec, champ)
                job = 'ACTU_VI'
                call xmmsa6(ndim, ipgf, zi(jmate), lamb, wsaut, &
                            nd, tau1, tau2, cohes, job, &
                            rela, alpha, dsidep, sigma, p, &
                            am3, rbid)
                do i = 1, ncompv
                    zr(jcoheo+ncompv*(ino-1)-1+i) = alpha(i)
                end do
                eps = r8prem()
                ASSERT((alpha(1)+eps) .ge. cohes(1))
!
            end do
!
!        SI CONTACT CLASSIQUE
!
        else if (contac .eq. 3 .or. contac .eq. 1) then
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
                    dn = 0.d0
                    if (algocr .eq. 3) then
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
                                rr, singu, tau1, tau2, ka, &
                                mu, jbaslo, jstno, jlsn, fk)
!
!            CALCUL COMPOSANTE NORMALE SAUT DE DEPLACEMENT
!
                    nvec = 1
                    call xmmsa3(ndim, nno, nnos, ffp, nddl, &
                                nvec, zr(idepl), zr(idepl), zr(idepl), nfh, &
                                singu, fk, ddls, ddlm, jheavn, &
                                ncompn, nfiss, ifiss, jheafa, ncomph, &
                                ifa, saut)
!
                    if (algocr .eq. 1 .or. algocr .eq. 2) then
                        gliss = zi(jgliss-1+nbspg+isspg)
                        memco = zi(jmemco-1+nbspg+isspg)
                        do j = 1, ndim
                            dn = dn+saut(j)*nd(j)
                        end do
                        nvec = 1
                        call xxlagm(ffc, idepl, ibid, lact, ndim, &
                                    nnol, pla, reac, r3bid, tau1, &
                                    tau2, nvec)
!
                        if (indco .eq. 0) then
!
!              ON REGARDE LA DISTANCE DN DES POINTS SUPPOSÉS
!              NON CONTACTANTS :
!              INTERPÉNÉPRATION EQUIVAUT À DN > 0 (ICI DN > 1E-16 )
!
                            if (dn .gt. prec) then
                                zi(jout2-1+nbspg+isspg) = 1
                                zi(jout3-1+nbspg+isspg) = 1
                                incoca = 1
                            else
                                zi(jout2-1+nbspg+isspg) = 0
                            end if
!
!                ON REGARDE LA REACTION POUR LES POINTS
!                SUPPOSES CONTACTANT :
                        else if (indco .eq. 1) then
                            if (coefcr .eq. 0.d0) rhon = 100.d0
                            if (coefcr .ne. 0.d0) rhon = coefcr
                            if ((reac-rhon*dn) .gt. r8prem()) then
!                  SI GLISSIERE=OUI ET IL Y A EU DU CONTACT DEJA SUR CE
!                  POINT (MEMCON=1), ALORS ON FORCE LE CONTACT
                                if ((gliss .eq. 1) .and. (memco .eq. 1)) then
                                    zi(jout2-1+nbspg+isspg) = 1
                                    zi(jout3-1+nbspg+isspg) = 1
                                else if (gliss .eq. 0) then
                                    zi(jout2-1+nbspg+isspg) = 0
                                    incoca = 1
                                end if
                            else
                                zi(jout2-1+nbspg+isspg) = 1
                                zi(jout3-1+nbspg+isspg) = 1
                            end if
                        else
!                SI INDCO N'EST NI ÉGAL À 0 NI ÉGAL À 1:
!                PROBLEME DE STATUT DE CONTACT.
                            ASSERT(indco .eq. 0 .or. indco .eq. 1)
                        end if
!
!           IMPRESSION (2EME PARTIE)
                        if (imprim) then
                            write (6, 698) indco, dn, reac, zi(jout2-1+nbspg+ &
                                                               isspg)
698                         format(5x, i1, 4x, 1pe12.5, 4x, 1pe12.5, 4x, i1)
                        end if
!
                    else if (algocr .eq. 3) then
!
!              CALCUL SAUT DE DEPLACEMENT EQUIVALENT
                        if (rela .eq. 1.d0 .or. rela .eq. 2.d0) then
                            job = 'ACTU_VI'
                            call xmmsa2(ndim, ipgf, zi(jmate), saut, nd, &
                                        tau1, tau2, cohes, job, rela, &
                                        alpha, dsidep, sigma, pp, dnor, &
                                        dtang, p, am3)
                        else if (rela .eq. 3.d0 .or. rela .eq. 4.d0) then
! NOUVEAUTE, IL FAUT RENSEIGNER LAMBDA
                            nvec = 1
                            call xxlag2(ffc, idepl, ibid, lact, ndim, &
                                        nnol, pla, lamb, nvec)
                            job = 'ACTU_VI'
                            call xmmsa5(ndim, ipgf, zi(jmate), saut, lamb, &
                                        nd, tau1, tau2, cohes, job, &
                                        rela, alpha, mat6bd, r6bid, mat3bd, &
                                        r3bd, rbid)
                        end if
!
! --- ACTUALISATION VARIABLE INTERNE
!
                        if (algocr .eq. 3) then
                            do i = 1, ncompv
                                zr(jcoheo+ncompv*(nbspg+isspg-1)-1+i) = &
                                    alpha(i)
                            end do
                            eps = r8prem()
                            ASSERT((alpha(1)+eps) .ge. cohes(1))
                        end if
!
! CHAMPS LIES AU CONTACT INUTILES POUR LE COHESIF
!
                        zi(jout2-1+nbspg+isspg) = 0
                        zi(jout3-1+nbspg+isspg) = 0
                        incoca = 0
                    end if
!
                end do
            end do
            nbspg = nbspg+nspfis
            jbasec = jbasec+ncompb
            jptint = jptint+ncompp
            jaint = jaint+ncompa
            jcface = jcface+ncompc
        end if
90      continue
    end do
!
!     ENREGISTREMENT DES CHAMPS DE SORTIE
    zi(jout1) = incoca
!
end subroutine
