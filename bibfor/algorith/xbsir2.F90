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

subroutine xbsir2(elref, contac, ddlc, ddlm, ddls, &
                  igeom, jheavn, jlst, ivectu, singu, &
                  nddl, ndim, nfh, nfiss, &
                  nno, nnom, nnos, depref, sigref, &
                  jbaslo, jstno, jlsn)
!
! aslint: disable=W1504
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elelin.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/xmprep.h"
#include "asterfort/xmulco.h"
#include "asterfort/xmvco3.h"
! IN ELREF  : ELEMENT DE REF PARENT
! IN CONTAC : DISCRETISATION, 1 POUR P1P1, 3 POUR P2P1
! IN DDLC   : NB DDL DE CONTACT PAR NOEUD SOMMET
! IN DDLM   : NB DDL PAR NOEUD MILIEU
! IN DDLS   : NB DDL TOT PAR NOEUD SOMMET
! IN GEOM   : ADRESSE POUR COORDONNEES NOEUD PARENT
! IN JLST   : ADRESSE LST
! IN/OUT IVECTU : VECTEUR RESIDUS DE REF
! IN SINGU
! IN NDDL   : NB TOTAL DDL ELEMENT
! NDIM      : DIMENSION DU MODELE
! NFH       : IDEM HEAVISIDE
! NFISS     : NB FISSURES
! NNO       : NB NOEUD ELEM PARENT
! NNOM      : DONT NB NOEUDS MILIEUX
! NNOS      : DONT NB NOEUDS SOMMETS
! DEPREF    : DEPLACEMENT DE REFERENCE
! SIGREF    : CONTRAINTE DE REFERENCE
! NOMTE     : TYPE D ELEMENT
! -------------------
! CALCUL RESIDU DE REFERENCE ELEMENTS COHESIF MIXTE XFEM
! TERMES D INTERFACE
! -------------------
    integer :: cface(30, 6), contac, ddlc, ddlm, ddls
    integer :: i, iadzi, iazk24, vstnc(1), ibid, ifa, ifiss, igeom, ipgf
    integer :: iret, jaint, jbasec, jcface
    integer :: jheavn, jheafa, jheano, jlonch, jlst, jptint, jtab(7)
    integer :: ivectu, lact(8), singu, jbaslo, jstno, jlsn
    integer :: nbspg, ncompa, ncompb, ncompc, ncomph, ncompp, ncompn
    integer :: nddl, ndim, nface, nfh, nfiss, ninter, nlact
    integer :: nno, nnol, nnom, nnos, npgf, nptf, nspfis, pla(27)
    integer :: idfdef, ipoidf, ivff, j, nnof
    real(kind=8) :: depref, ffc(8), ffp(27), jac
    real(kind=8) :: r3bid(3), rr, sigref, vtmp(400)
    real(kind=8) :: fk(27, 3, 3)
    aster_logical :: lbid
    character(len=8) :: elc, elref, elrefc, fpg, typma
!
! --- INITIALISATIONS
!
    do i = 1, 8
        lact(i) = 0
    end do
    ffp(:) = 0.d0
    vtmp(:) = 0.d0
    rr = 0.d0
    ncomph = 0
    nbspg = 0
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
!
! --- ROUTINE SPECIFIQUE P2P1
!
    call elelin(contac, elref, elrefc, ibid, ibid)
!
! --- ARGUMENTS SUPPLEMENTAIRES NECESSAIRES PAR RAPPORT
! --- AUX ELEMENTS VOLUMIQUES
!
    call jevech('PPINTER', 'L', jptint)
    call jevech('PAINTER', 'L', jaint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PBASECO', 'L', jbasec)
    call jevech('PLONFA', 'L', jlonch)
!
!   RECUPERATION DE LA DEFINITION DES FONCTIONS HEAVISIDES
    if (nfh .gt. 0) then
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
    end if
    if (nfiss .gt. 1) then
        call jevech('PHEAVNO', 'L', jheano)
        call jevech('PHEA_FA', 'L', jheafa)
        call tecach('OOO', 'PHEA_FA', 'L', iret, nval=2, &
                    itab=jtab)
        ncomph = jtab(2)
    end if
!     DIMENSIONS DES GRANDEURS DANS LA CARTE
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
!
! --- BOUCLE SUR LES FISSURES
!
    do ifiss = 1, nfiss
!
! --- RECUPERATION DIVERSES DONNEES CONTACT
!
        ninter = zi(jlonch+3*(ifiss-1)-1+1)
        if (ninter .eq. 0) goto 90
!
        fpg = 'FPG2'
! SCHEMA EN DUR POUR LE MOMENT
        if (ndim .eq. 3) then
            elc = 'TR3'
        else if (ndim .eq. 2) then
            if (contac .le. 2) then
                elc = 'SE2'
            else
                elc = 'SE3'
            end if
        end if
!
        call elrefe_info(elrefe=elc, fami=fpg, nno=nnof, npg=npgf, jpoids=ipoidf, &
                         jvf=ivff, jdfde=idfdef)
        nface = zi(jlonch+3*(ifiss-1)-1+2)
        nptf = zi(jlonch+3*(ifiss-1)-1+3)
        do i = 1, nface
            do j = 1, nptf
                cface(i, j) = zi(jcface-1+nptf*(i-1)+j)
            end do
        end do
!
        nspfis = npgf*nface
!
!
! --- RECUP MULTIPLICATEURS ACTIFS ET LEURS INDICES
!
        call xmulco(contac, ddls, ddlc, ddlm, jaint, ifiss, &
                    jheano, vstnc, lact, .false._1, lbid, &
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
! --- PREPARATION DU CALCUL
!
                call xmprep(cface, contac, elref, elrefc, elc, &
                            ffc, ffp, fpg, jaint, jbasec, &
                            jptint, ifa, igeom, ipgf, jac, &
                            jlst, lact, r3bid, ndim, ninter, &
                            nlact, nno, nnos, nptf, ibid, &
                            rr, singu, r3bid, r3bid, 3.d0, sigref, &
                            jbaslo, jstno, jlsn, fk)
!
! --- CALCUL VECTEURS DE REFERENCE POUR LA LOI D INTERFACE
!
                call xmvco3(sigref, depref, ndim, nno, nnol, &
                            nnos, pla, lact, nfh, ddls, &
                            ddlm, nfiss, ifiss, jheafa, ifa, &
                            ncomph, jheavn, ncompn, jac, ffc, ffp, &
                            singu, fk, vtmp)
! --- FIN DE BOUCLE SUR LES POINTS DE GAUSS
            end do
!
! --- FIN DE BOUCLE SUR LES FACETTES
        end do
! --- FIN BOUCLE SUR LES FISSURES
        nbspg = nbspg+nspfis
        jbasec = jbasec+ncompb
        jptint = jptint+ncompp
        jaint = jaint+ncompa
        jcface = jcface+ncompc
90      continue
    end do
!
!-----------------------------------------------------------------------
!     COPIE DES CHAMPS DE SORTIES ET FIN
!-----------------------------------------------------------------------
!
    do i = 1, nddl
        zr(ivectu-1+i) = zr(ivectu-1+i)+vtmp(i)
    end do
!
end subroutine
