! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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

subroutine te0363(option, nomte)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elelin.h"
#include "asterfort/elrfvf.h"
#include "asterfort/jevech.h"
#include "asterfort/mmnorm.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/xlacti.h"
#include "asterfort/xmelet.h"
#include "asterfort/xmmjec.h"
#include "asterfort/xmoffc.h"
#include "asterfort/xtform.h"
#include "asterfort/xtlagc.h"
#include "asterfort/tecach.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xkamat.h"
    character(len=16) :: option, nomte
!
! ----------------------------------------------------------------------
!  CONTACT XFEM GRANDS GLISSEMENTS
!  REACTUALISATION DU STATUT DE CONTACT
!
!  OPTION : 'XCVBCA' (X-FEM MISE À JOUR DU STATUT DE CONTACT)
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!
!
!
    character(len=8) :: typmai, elrees, elrema, elreco, typmec
    integer :: ndim, nddl, nne(3), nnm(3), nnc
    integer :: nsinge, nsingm, lact(8), nlact, ninter
    integer :: jpcpo, jpcai, jheafa, jheano
    integer :: jgeom, jdepde, jheavn, ncompn, jtab(7), iret
    integer :: indco, memco, indnor, igliss, nfaes
    integer :: jout
    integer :: incoca, nfhe, nfhm
    real(kind=8) :: jeuca, tau1(3), tau2(3), norm(3)
    real(kind=8) :: coore(3), coorm(3), coorc(2)
    real(kind=8) :: dlagrc, rhon
    real(kind=8) :: ffc(9), ffe(20), ffm(20), dffc(3, 9)
    real(kind=8) :: prec, noor
    real(kind=8) :: rre, rrm, ffec(8)
    parameter    (prec=1.d-16)
    integer :: contac, ddle(2), ddlm(2), ibid, ndeple
    integer :: imate, jbaslo, jlsn, jstno
    real(kind=8) :: fk_escl(27,3,3), fk_mait(27,3,3), ka, mu
    aster_logical :: lmulti
!
! ----------------------------------------------------------------------
!
!
! --- INFOS SUR LA MAILLE DE CONTACT
!
    call xmelet(nomte, typmai, elrees, elrema, elreco,&
                ndim, nddl, nne, nnm, nnc,&
                ddle, ddlm, contac, ndeple, nsinge,&
                nsingm, nfhe, nfhm)
!
    ASSERT(nddl.le.336)
    lmulti = .false.
    if (nfhe .gt. 1 .or. nfhm .gt. 1) lmulti = .true.
!
! --- INITIALISATIONS
!
    dlagrc = 0.d0
    fk_escl(:,:,:) = 0.d0
    fk_mait(:,:,:) = 0.d0
! --- INITIALISATION DE LA VARIABLE DE TRAVAIL
    incoca = 0
    ASSERT(option.eq.'XCVBCA')
!
! --- RECUPERATION DES DONNEES DE LA CARTE CONTACT 'POINT' (VOIR XMCART)
!
    call jevech('PCAR_PT', 'L', jpcpo)
! --- LES COORDONNEES ESCLAVE DANS L'ELEMENT DE CONTACT
    coorc(1) = zr(jpcpo-1+1)
    coorc(2) = zr(jpcpo-1+10)
    tau1(1) = zr(jpcpo-1+4)
    tau1(2) = zr(jpcpo-1+5)
    tau1(3) = zr(jpcpo-1+6)
    tau2(1) = zr(jpcpo-1+7)
    tau2(2) = zr(jpcpo-1+8)
    tau2(3) = zr(jpcpo-1+9)
    rhon = zr(jpcpo-1+13)
    indco = nint(zr(jpcpo-1+11))
    ninter = nint(zr(jpcpo-1+31))
    indnor = nint(zr(jpcpo-1+17))
    igliss = nint(zr(jpcpo-1+20))
    memco = nint(zr(jpcpo-1+21))
    nfaes = nint(zr(jpcpo-1+22))
! --- LES COORDONNEES ESCLAVE ET MAITRES DANS L'ELEMENT PARENT
    coore(1) = zr(jpcpo-1+24)
    coore(2) = zr(jpcpo-1+25)
    coore(3) = zr(jpcpo-1+26)
    coorm(1) = zr(jpcpo-1+27)
    coorm(2) = zr(jpcpo-1+28)
    coorm(3) = zr(jpcpo-1+29)
! --- SQRT LSN PT MAITRE, ESCLAVE
    rre = zr(jpcpo-1+18)
    rrm = zr(jpcpo-1+23)
    if (nnm(1) .eq. 0) rre = 2*rre
!
! --- RECUPERATION DES DONNEES DE LA CARTE CONTACT AINTER (VOIR XMCART)
!
    call jevech('PCAR_AI', 'L', jpcai)
!
! --- RECUPERATION DE LA GEOMETRIE ET DES CHAMPS DE DEPLACEMENT
!
    call jevech('PGEOMER', 'L', jgeom)
    call jevech('PDEPL_P', 'L', jdepde)
!
    if (nfhe .gt. 0 .or. nfhm .gt. 0) then
      call jevech('PHEA_NO', 'L', jheavn)
      call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7,&
                itab=jtab)
      ncompn = jtab(2)/jtab(3)
    endif
!
    if (lmulti) then
!
! --- RECUPERATION DES FONCTION HEAVISIDES SUR LES FACETTES
!
        call jevech('PHEA_FA', 'L', jheafa)
!
! --- RECUPERATION DE LA PLACE DES LAGRANGES
!
        call jevech('PHEAVNO', 'L', jheano)
    else
        jheafa=1
        jheano=1
    endif
!
! --- RECUPERATION DES CHAMPS DE SORTIE
!
    call jevech('PINDCOO', 'E', jout)
!
! --- FONCTIONS DE FORMES
!
    call xtform(elrees, elrema, elreco,&
                nnm(1),coore, coorm, coorc,&
                ffe, ffm, dffc)
!
! --- CALCUL DES FCTS SINGULIERES
!
    call jevech('PSTANO', 'L', jstno)
    if (nsinge.eq.1 .and. nne(1).gt. 0) then
      call jevech('PLSNGG', 'L', jlsn)
      call jevech('PBASLOC', 'L', jbaslo)
      call jevech('PMATERC', 'L', imate)
      call xkamat(zi(imate), ndim, .false._1, ka, mu)
      call xcalfev_wrap(ndim, nne(1), zr(jbaslo), zi(jstno), -1.d0,&
                   zr(jlsn), zr(jlsn), zr(jgeom), ka, mu, ffe, fk_escl, face='ESCL')
    endif
!
! --- BRICOLAGES POUR RESPECTER LES ANCIENNES CONVENTIONS DE SIGNE
    fk_escl=-1.d0*fk_escl
    if (nnm(1) .eq. 0) fk_escl=2.d0*fk_escl
    if (nsingm.eq.1 .and. nnm(1).gt.0) then
      call jevech('PLSNGG', 'L', jlsn)
      call jevech('PBASLOC', 'L', jbaslo)
      call jevech('PMATERC', 'L', imate)
      call xkamat(zi(imate), ndim, .false._1, ka, mu)
      call xcalfev_wrap(ndim, nnm(1), zr(jbaslo+nne(1)*3*ndim), zi(jstno+nne(1)), +1.d0,&
                   zr(jlsn+nne(1)), zr(jlsn+nne(1)),zr(jgeom+nne(1)*ndim), &
                   ka, mu, ffm, fk_mait, face='MAIT')
    endif
!
! --- FONCTION DE FORMES POUR LES LAGRANGIENS
! --- SI ON EST EN LINEAIRE, ON IMPOSE QUE LE NB DE NOEUDS DE CONTACTS
! --- ET LES FFS LAGRANGES DE CONTACT SONT IDENTIQUES A CEUX
! --- DES DEPLACEMENTS DANS LA MAILLE ESCLAVE POUR LE CALCUL DES CONTRIB
!
    if (contac .eq. 1) then
        nnc = nne(2)
        call xlacti(typmai, ninter, jpcai, lact, nlact)
        call xmoffc(lact, nlact, nnc, ffe, ffc)
    else if (contac.eq.3) then
        nnc = nne(2)
        call elelin(contac, elrees, typmec, ibid, ibid)
        call elrfvf(typmec, coore, ffec)
        call xlacti(typmai, ninter, jpcai, lact, nlact)
        call xmoffc(lact, nlact, nnc, ffec, ffc)
    else
        ASSERT(contac.eq.0)
    endif
!
! --- CALCUL DE LA NORMALE
!
    if (ndim .eq. 2) then
        call mmnorm(ndim, tau1, tau2, norm, noor)
    else if (ndim.eq.3) then
        call provec(tau1, tau2, norm)
        call normev(norm, noor)
    endif
    if (noor .le. r8prem()) then
        ASSERT(.false.)
    endif
!
! --- CALCUL DE L'INCREMENT DE REACTION DE CONTACT
!
    call xtlagc(ndim, nnc, nne, ddle(1),&
                jdepde, ffc,&
                nfhe, lmulti, zi(jheano), dlagrc)
!
! --- EVALUATION DES JEUX - CAS DU CONTACT
!
    call xmmjec(ndim, nnm, nne, ndeple, nsinge,&
                nsingm, ffe, ffm, norm, jgeom,&
                jdepde, fk_escl, fk_mait, ddle, ddlm,&
                nfhe, nfhm, lmulti, zi(jheavn), zi(jheafa), jeuca)
!
!
! --- NOEUDS EXCLUS PAR PROJECTION HORS ZONE
!
    if (indnor .eq. 1) then
        if ((igliss.eq.0) .or. (memco.eq.0)) then
            zi(jout-1+2) = 0
            goto 999
        endif
    endif
!
! --- SI LE CONTACT A ETE POSTULE, ON TESTE LA VALEUR DE LA PRESSION
! --- DE CONTACT
!
    if (indco .eq. 1) then
        if ((dlagrc-rhon*jeuca) .gt. r8prem()) then
            if ((igliss.eq.1) .and. (memco.eq.1)) then
                zi(jout-1+2) = 1
                zi(jout-1+3) = 1
            else if (igliss.eq.0) then
                zi(jout-1+2) = 0
                incoca = 1
            endif
        else
            zi(jout-1+2) = 1
            zi(jout-1+3) = 1
        endif
!
! --- SI LE NON-CONTACT A ETE POSTULÉ, ON TESTE LA VALEUR DU JEU
!
    else if (indco .eq. 0) then
        if (jeuca .gt. prec) then
            zi(jout-1+2) = 1
            zi(jout-1+3) = 1
            incoca = 1
        else
            zi(jout-1+2) = 0
        endif
!
    else
        ASSERT(.false.)
    endif
!
999 continue
!
! --- ENREGISTREMENT DU CHAMP DE SORTIE
!
    zi(jout-1+1)=incoca
!
end subroutine
