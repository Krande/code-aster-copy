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
! person_in_charge: daniele.colombo at ifpen.fr
!
subroutine te0557(option, nomte)
!
    use THM_type
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
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/teattr.h"
#include "asterfort/xasshv_frac.h"
#include "asterfort/xhmini.h"
#include "asterfort/xhmddl.h"
#include "asterfort/xmulhm.h"
#include "asterfort/xminte.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!         CALCUL DES SECONDS MEMBRES DE CONTACT POUR MODELE HM-XFEM
!
!  OPTION : 'CHAR_MECA_CONT' (CALCUL DES SECONDS MEMBRES DE CONTACT)
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer :: i, j
    integer :: nfh, ddld, ddlm, ddlp, ddlc, nddls, nddlm
    integer :: nnop, nnops, nnopm, dimuel, jheavn, ncompn, pos(16)
    integer :: lact(16), vstnc(32), nlact(2), algocr, jheafa, ncomph
    integer :: ndim, nbspg, contac, ibid, iadzi, iazk24, ifisc, nfisc
    integer :: igeom, idepm, idepd, jlsn, jfisco
    integer :: jptint, jaint, jcface, jlonch, jbasec, jdonco
    integer :: iinstm, iinstp, icarcr, jmate, jcohes
    integer :: jtab(7), iret, ncompv, nnol, pla(27), jv_cont, ncompp
    integer :: npgf, ipoidf, ivff, idfdef, jstno, jfisno, jheano
    integer :: ninter, nface, nptf, nnof, cface(30, 6), ncompd
    integer :: nint, ninteg, ifiss, nfiss, ncompa, ncompb, ncompc
    integer :: jfiss, nfisc2, kfiss
    aster_logical :: lelim
    real(kind=8) :: rela, vcont(560), mat(1)
    character(len=8) :: elrefp, elrefc, elc, fpg, typma
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: enr
!
    integer :: nfimax
    parameter(nfimax=10)
    integer :: fisc(2*nfimax), fisco(2*nfimax)
    type(THM_DS) :: ds_thm
!
! --------------------------------------------------------------------------------------------------
!

!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm)
    rela = 0.d0
    nbspg = 0

!   INTIALISATION JMATE POUR DETECTER EVENTUELLES ERREURS JEVEUX
    jmate = 1
!
!   INITIALISATION DU NOMBRE DE DDL PAR NEOUD, DU TYPE DE CONTACT ET
!   DES ADRESSES POUR LES DIFFERENTS TERMES DE L'OPERATEUR TANGENT
!   (CAS DE LA FRACTURE UNIQUEMENT)
    vcont(:) = 0.d0
    call xhmini(nomte, nfh, ddld, ddlm, ddlp, nfiss, ddlc, contac)
!
!   INITIALISATION DE LA DIMENSION DE L'ELEMENT PRINCIPAL, DU NOMBRE DE
!   NOEUDS PARENTS (SOMMET + MILIEU)
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop, nnos=nnops)
!
    nddls = ddld+ddlp+ddlc
    nddlm = ddlm
    nnopm = nnop-nnops
    dimuel = nnops*nddls+nnopm*nddlm
!
!   RECUPERATION DU TYPE DE MAILLE POUR LA SELECTION DES LAGRANGES ACTIFS
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
!
    if (nfh .gt. 0) then
        call jevech('PHEA_NO', 'L', jheavn)
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
    end if
    jfisno = 1
    jheano = 1
    jfisco = 1
    jheafa = 1
    jlsn = 1
    ncomph = 1
!
    do i = 1, 2*nfimax
        fisco(i) = 0
        fisc(i) = 0
    end do
!
    if (nfiss .gt. 1) then
        call jevech('PFISNO', 'L', jfisno)
        call jevech('PHEAVNO', 'L', jheano)
        call jevech('PHEA_FA', 'L', jheafa)
        call tecach('OOO', 'PHEA_FA', 'L', iret, nval=2, &
                    itab=jtab)
        ncomph = jtab(2)
        call jevech('PFISCO', 'L', jfisco)
        do i = 1, 2*nfiss
            fisco(i) = zi(jfisco-1+i)
        end do
        call jevech('PLSN', 'L', jlsn)
    end if
!
!   RECUPERATION DES ADRESSES DES CATALOGUES POUR LES ELEMENTS PRINCIPAUX
!   HM-XFEM -- PARAMETRES EN ENTREE
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPL_M', 'L', idepm)
    call jevech('PDEPL_P', 'L', idepd)
    call jevech('PDONCO', 'L', jdonco)
    call jevech('PPINTER', 'L', jptint)
    call jevech('PAINTER', 'L', jaint)
    call jevech('PLONGCO', 'L', jlonch)
    call jevech('PBASECO', 'L', jbasec)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)

    call jevech('PSTANO', 'L', jstno)
!
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
!
!   ON RECUPERE L'ELEMENT PRINCIPAL ...
    call elref1(elrefp)
!   ... AINSI QUE L'ELEMENT LINEAIRE ASSOCIE, POUR LE CONTACT P2P1
    if (contac .ge. 2) then
        call elelin(contac, elrefp, elrefc, ibid, ibid)
    else
        ASSERT(.false.)
    end if
!
!   STATUT POUR L'ELIMINATION DES PRE_FLU, LAG_FLI, LAG_FLS, LAG1_HM ET
!   LAG2_HM PORTES PAR UNE ARETE NON VITALE
    do i = 1, 32
        vstnc(i) = 1
    end do
!
!   LOGICAL POUR L'ELIMINATION DES LAGRANGES EN TROP
    lelim = .false.
!
    do ifiss = 1, nfiss
!
        do i = 1, 2*nfiss
            fisc(i) = 0
        end do
        ifisc = ifiss
        nfisc = 0
80      continue
        if (fisco(2*ifisc-1) .gt. 0) then
!   STOCKAGE DES FISSURES SUR LESQUELLES IFISS SE BRANCHE
            nfisc = nfisc+1
            fisc(2*(nfisc-1)+2) = fisco(2*ifisc)
            ifisc = fisco(2*ifisc-1)
            fisc(2*(nfisc-1)+1) = ifisc
            goto 80
        end if
!
        nfisc2 = 0
        do jfiss = ifiss+1, nfiss
!   STOCKAGE DES FISSURES QUI SE BRANCHENT SUR IFISS
            kfiss = fisco(2*jfiss-1)
            do i = nfisc+1, nfisc+nfisc2
                if (fisc(2*(i-1)+1) .eq. kfiss) then
                    nfisc2 = nfisc2+1
                    fisc(2*(nfisc+nfisc2-1)+1) = jfiss
                    fisc(2*(nfisc+nfisc2)) = fisco(2*jfiss)
                end if
            end do
            if (kfiss .eq. ifiss) then
                nfisc2 = nfisc2+1
                fisc(2*(nfisc+nfisc2-1)+1) = jfiss
                fisc(2*(nfisc+nfisc2)) = fisco(2*jfiss)
            end if
        end do
    end do
!
    do ifiss = 1, nfiss
!
!   INITIALISATION DU VECTEUR DES PRE_FLU, LAG_FLI, LAG_FLS, LAG1_HM ET
!   LAG2_HM ACTIFS PORTES PAR UNE ARETE VITALE
        do i = 1, 16
            lact(i) = 0
            pos(i) = 0
        end do
        do i = 1, 27
            pla(i) = 0
        end do
!
!   DEFINITION DU NOMBRE DE POINTS D'INTERSECTION, DU NOMBRE FACETTE ET
!   DU NOMBRE DE POINT PAR FACETTE
        ninter = zi(jlonch+3*(ifiss-1)-1+1)
        nface = zi(jlonch+3*(ifiss-1)-1+2)
        nptf = zi(jlonch+3*(ifiss-1)-1+3)
!   SELECTION DES DDLS PRE_FLU, LAG_FLI, LAG_FLS, LAG1_HM ET
!   LAG2_HM ACTIFS
!
        call xmulhm(contac, nddls, ddlc, nddlm, jaint, &
                    ifiss, jheano, vstnc, lact, .true._1, lelim, &
                    nfh, nfiss, ninter, nlact, nnop, &
                    nnol, nnopm, nnops, pla, pos, typma, jstno)
!
!   SI IL N'Y A PAS DE FACETTES POUR LA FRACTURE ON SORT
        if (ninter .eq. 0) goto 200
!
!   RECUPERATION DU TYPE D'ELEMENT POUR LA FACETTE DE CONTACT ET
!   DE LA FAMILLE DE POINT D'INTEGRATION (P2P1 UNIQUEMENT)
        if ((ndim .eq. 2) .and. (contac .ge. 2)) then
            elc = 'SE3'
            ninteg = nint(zr(jdonco-1+(ifiss-1)*ncompd+4))
            call xminte(ndim, ninteg, fpg)
        elseif ((ndim .eq. 3) .and. (contac .ge. 2)) then
            elc = 'TR6'
            ninteg = nint(zr(jdonco-1+(ifiss-1)*ncompd+4))
            call xminte(ndim, ninteg, fpg)
        end if
!
!   RECUPERATION DES DONNEES RELATIVES AU CONTACT AVEC LOI COHESIVE
!
        algocr = nint(zr(jdonco-1+(ifiss-1)*ncompd+6))
!
        if (algocr .eq. 3) then
            call teattr('S', 'XFEM', enr, ibid, typel=nomte)
            ASSERT(enr(3:3) .eq. 'C' .or. enr(4:4) .eq. 'C')
            rela = zr(jdonco-1+(ifiss-1)*ncompd+10)
            if (contac .eq. 3) then
                ASSERT((rela .eq. 3.d0) .or. (rela .eq. 4.d0))
            elseif (contac .eq. 2) then
                ASSERT(rela .eq. 5.d0)
            end if
            call jevech('PMATERC', 'L', jmate)
            call jevech('PCOHES', 'L', jcohes)
            call tecach('OOO', 'PCOHES', 'L', iret, nval=3, &
                        itab=jtab)
            if (contac .eq. 2) ncompv = jtab(2)/jtab(3)
            if (contac .eq. 3) ncompv = jtab(2)
        else
            ASSERT(.false.)
        end if

!   RECUPERATION DES DIFFERENTES ADRESSES POUR L'INTEGRATION SUR LES
!   FACETTES DE CONTACT

        call elrefe_info(elrefe=elc, fami=fpg, nno=nnof, &
                         npg=npgf, jpoids=ipoidf, jvf=ivff, jdfde=idfdef)

!   DEFINTION DE LA CONNECTIVITE DES FACETTES DE CONTACT
!
        do i = 1, 30
            do j = 1, 6
                cface(i, j) = 0
            end do
        end do
!
        do i = 1, nface
            do j = 1, nptf
                cface(i, j) = zi(jcface-1+nptf*(i-1)+j)
            end do
        end do
!
!   CALCUL DES SECONDS MEMBRES POUR LA FRACTURE
!
        call xasshv_frac(ds_thm, &
                         nddls, nddlm, nnop, nnops, &
                         lact, elrefp, elrefc, elc, contac, &
                         dimuel, nface, npgf, nbspg, nptf, &
                         jcohes, jptint, igeom, jbasec, &
                         nlact, cface, zr(iinstp), &
                         zr(iinstm), zr(icarcr), fpg, ncompv, &
                         vcont, compor, jmate, ndim, idepm, &
                         idepd, pla, algocr, rela, jheavn, ncompn, &
                         ifiss, nfiss, nfh, jheafa, ncomph, pos)
!
        nbspg = nbspg+npgf*nface
!
200     continue
        jbasec = jbasec+ncompb
        jptint = jptint+ncompp
        jaint = jaint+ncompa
        jcface = jcface+ncompc
    end do
!
    call jevech('PVECTCR', 'E', jv_cont)
    do i = 1, dimuel
        zr(jv_cont-1+i) = vcont(i)
    end do
!
!   SUPPRESSION DES DDLS DE DEPLACEMENT ET DE PRESSION HPRE1

    if (nfh .ne. 0) then
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    zi(jstno), .false._1, option, nomte, &
                    mat, zr(jv_cont), nddlm, nfiss, jfisno, .false._1, contac)
    end if

!   SUPPRESSION DES DDLS DE CONTACT

    if (lelim) then
        call xhmddl(ndim, nfh, nddls, dimuel, nnop, nnops, &
                    vstnc, .false._1, option, nomte, &
                    mat, zr(jv_cont), nddlm, nfiss, jfisno, .true._1, contac)
    end if
!
end subroutine
