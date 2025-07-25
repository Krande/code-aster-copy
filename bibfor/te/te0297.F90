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
subroutine te0297(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/cgverho.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/iselli.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/ltequa.h"
#include "asterfort/rcvad2.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/xcgfvo.h"
#include "asterfort/xsifel.h"
#include "asterfort/xsifle.h"
#include "asterfort/xteini.h"
!
    character(len=16) :: option, nomte
!
! person_in_charge: samuel.geniaut at edf.fr
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE POST-TRAITEMENT
!                          EN MÉCANIQUE DE LA RUPTURE
!                          POUR LES ÉLÉMENTS X-FEM
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
    integer(kind=8) :: ndim, nno, nnop, npg, ier
    integer(kind=8) :: nfh, nfe, ddlc, nse, ise, in, ino, ninter
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jbaslo, igeom, idepl
    integer(kind=8) :: ipres, ipref, itemps, jptint, jcface, jlongc, imate
    integer(kind=8) :: ithet, i, j, compt, igthet, ibid, jlsn, jlst, idecpg, icode
    integer(kind=8) :: nface, cface(30, 6), ifa, singu, jpmilt, ipuls, iret, jtab(7)
    integer(kind=8) :: irese, ddlm, jbasec, nptf, nfiss, jheavn, jstno
    integer(kind=8) :: contac
    real(kind=8) :: thet, valres(3), devres(3), presn(27), valpar(4)
    real(kind=8) :: pres, fno(81), coorse(81), puls
    integer(kind=8) :: icodre(3)
    character(len=8) :: elrefp, elrese(6), fami(6), nompar(4), enr
    character(len=16) :: nomres(3)
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
    data nomres/'E', 'NU', 'ALPHA'/
!
!
    call elref1(elrefp)
    call jevech('PTHETAR', 'L', ithet)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
!
!   SI LA VALEUR DE THETA EST NULLE SUR L'ÉLÉMENT, ON SORT
    compt = 0
    do i = 1, nnop
        thet = 0.d0
        do j = 1, ndim
            thet = thet+abs(zr(ithet+ndim*(i-1)+j-1))
        end do
        if (thet .lt. r8prem()) compt = compt+1
    end do
    if (compt .eq. nnop) goto 999
!
!   SOUS-ELEMENT DE REFERENCE : RECUP DE NNO, NPG ET IVF
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, npg=npg)
!
!   INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                ibid, ibid, ibid, ddlm, nfiss, &
                contac)
!
!     ------------------------------------------------------------------
!              CALCUL DE G, K1, K2, K3 SUR L'ELEMENT MASSIF
!     ------------------------------------------------------------------
!
!     PARAMÈTRES PROPRES À X-FEM
    if (option .ne. 'CALC_K_G_COHE') then
        call jevech('PPINTTO', 'L', jpintt)
        call jevech('PCNSETO', 'L', jcnset)
        call jevech('PHEAVTO', 'L', jheavt)
        call jevech('PLONCHA', 'L', jlonch)
    end if
    call jevech('PBASLOR', 'L', jbaslo)
    call jevech('PLST', 'L', jlst)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PMATERC', 'L', imate)
    call jevech('PGTHETA', 'E', igthet)
    call teattr('S', 'XFEM', enr, ier)
    if (nfh .gt. 0) call jevech('PHEA_NO', 'L', jheavn)
!   Propre aux elements 1d et 2d (quadratiques)
    if ((ier .eq. 0) .and. ltequa(elrefp, enr)) call jevech('PPMILTO', 'L', jpmilt)
    if (nfe .gt. 0) call jevech('PSTANO', 'L', jstno)
    if (option .eq. 'CALC_K_G_COHE') goto 98
    call jevech('PLSN', 'L', jlsn)
!
!   VERIFS DE COHERENCE RHO <-> PESANTEUR, ROTATION, PULSATION
    if (.not. cgverho(imate)) call utmess('F', 'RUPTURE1_26')
!
!   CALCUL DES FORCES NODALES CORRESPONDANT AUX CHARGES VOLUMIQUES
    call xcgfvo(option, ndim, nnop, fno)
!
! --- RECUPERATION DE LA PULSATION
!
    call tecach('ONO', 'PPULPRO', 'L', iret, nval=7, &
                itab=jtab)
    ipuls = jtab(1)
    if (iret .eq. 0) then
        puls = zr(ipuls)
    else
        puls = 0.d0
    end if
!
!   Recuperation de la subdivision de l'element en nse sous element
    nse = zi(jlonch-1+1)
!
!   Boucle sur les nse sous-elements
    do ise = 1, nse
!
!       Boucle sur les sommets du sous-tria (du sous-seg)
        do in = 1, nno
            ino = zi(jcnset-1+nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000- &
                                                              1)+j)
                else if (ino .gt. 2000 .and. ino .lt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-2000- &
                                                              1)+j)
                else if (ino .gt. 3000) then
                    coorse(ndim*(in-1)+j) = zr(jpmilt-1+ndim*(ino-3000- &
                                                              1)+j)
                end if
            end do
        end do
!
        idecpg = npg*(ise-1)
!
        call xsifel(elrefp, ndim, coorse, igeom, jheavt, &
                    ise, nfh, ddlc, ddlm, nfe, &
                    puls, zr(jbaslo), nnop, idepl, zr(jlsn), &
                    zr(jlst), idecpg, igthet, fno, nfiss, &
                    jheavn, jstno)
!
    end do
!
!     ------------------------------------------------------------------
!              CALCUL DE G, K1, K2, K3 SUR LES LEVRES
!     ------------------------------------------------------------------
!
    if (option .eq. 'CALC_K_G_XFEM') then
!       SI LA PRESSION N'EST CONNUE SUR AUCUN NOEUD, ON LA PREND=0.
        call jevecd('PPRESSR', ipres, 0.d0)
    else if (option .eq. 'CALC_K_G_XFEM_F') then
        call jevech('PPRESSF', 'L', ipref)
        call jevech('PINSTR', 'L', itemps)
!
!       RECUPERATION DES PRESSIONS AUX NOEUDS PARENTS
        nompar(1) = 'X'
        nompar(2) = 'Y'
        if (ndim .eq. 3) nompar(3) = 'Z'
        if (ndim .eq. 3) nompar(4) = 'INST'
        if (ndim .eq. 2) nompar(3) = 'INST'
        do i = 1, nnop
            do j = 1, ndim
                valpar(j) = zr(igeom+ndim*(i-1)+j-1)
            end do
            valpar(ndim+1) = zr(itemps)
            call fointe('FM', zk8(ipref), ndim+1, nompar, valpar, &
                        presn(i), icode)
        end do
    end if
!
!   SI LA VALEUR DE LA PRESSION EST NULLE SUR L'ÉLÉMENT, ON SORT
    compt = 0
    do i = 1, nnop
        if (option .eq. 'CALC_K_G_XFEM') pres = abs(zr(ipres-1+i))
        if (option .eq. 'CALC_K_G_XFEM_F') pres = abs(presn(i))
        if (pres .lt. r8prem()) compt = compt+1
    end do
    if (compt .eq. nnop) goto 999
!
98  continue
!
!   PARAMETRES PROPRES A X-FEM
    call jevech('PPINTER', 'L', jptint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PLONGCO', 'L', jlongc)
    call jevech('PBASECO', 'L', jbasec)
!
!   RÉCUPÉRATIONS DES DONNÉES SUR LA TOPOLOGIE DES FACETTES
    ninter = zi(jlongc-1+1)
    nface = zi(jlongc-1+2)
    nptf = zi(jlongc-1+3)
    if (ninter .lt. ndim) goto 999
!
    do i = 1, 30
        do j = 1, 6
            cface(i, j) = 0
        end do
    end do
!
    do i = 1, nface
        do j = 1, nptf
            cface(i, j) = zi(jcface-1+ndim*(i-1)+j)
        end do
    end do
!
!   RECUPERATION DES DONNEES MATERIAU AU 1ER POINT DE GAUSS !!
!   LE MATÉRIAU DOIT ETRE HOMOGENE DANS TOUT L'ELEMENT
    call rcvad2('XFEM', 1, 1, '+', zi(imate), &
                'ELAS', 3, nomres, valres, devres, &
                icodre)
    if ((icodre(1) .ne. 0) .or. (icodre(2) .ne. 0)) then
        call utmess('F', 'RUPTURE1_25')
    end if
    if (icodre(3) .ne. 0) then
        valres(3) = 0.d0
        devres(3) = 0.d0
    end if
!
!   BOUCLE SUR LES FACETTES
    do ifa = 1, nface
        call xsifle(ndim, ifa, jptint, cface, igeom, &
                    nfh, jheavn, singu, nfe, ddlc, &
                    ddlm, jlsn, jlst, jstno, ipres, &
                    ipref, itemps, idepl, nnop, valres, &
                    zr(jbaslo), ithet, nompar, option, igthet, &
                    jbasec, contac)
    end do
!
!
999 continue
end subroutine
