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
subroutine te0288(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/iselli.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/ltequa.h"
#include "asterfort/rcvad2.h"
#include "asterfort/teattr.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/xcgfvo.h"
#include "asterfort/xgelem.h"
#include "asterfort/xsifle.h"
#include "asterfort/xteini.h"
    character(len=16) :: option, nomte
!
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE POST-TRAITEMENT
!                          EN MÉCANIQUE DE LA RUPTURE
!                          POUR LES ÉLÉMENTS X-FEM (OPTION CALC_G)
!
!             ATTENTION : PAS D'ETAT INITIAL
!                         PAS DE GRANDES DEF, GRANDES ROT
!                         PAS DE COMP_INCR
!
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
!
    integer(kind=8) :: ndim, nno, nnop, npg, nptf, jbasec, nfiss
    integer(kind=8) :: nfh, nfe, ddlc, nse, ise, in, ino
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jbaslo, igeom, idepl
    integer(kind=8) :: ipres, ipref, itemps, jptint, jcface, jlongc, imate, jheavn
    integer(kind=8) :: ithet, i, j, compt, igthet, ibid, jlsn, jlst, icode
    integer(kind=8) :: ninter, nface, cface(30, 6), ifa, singu, jpmilt, irese, ddlm
    real(kind=8) :: thet, valres(3), devres(3), presn(27), valpar(4)
    real(kind=8) :: pres, fno(81), coorse(81)
    integer(kind=8) :: icodre(3), contac, iadzi, iazk24, jstno
    character(len=8) :: elrefp, elrese(6), fami(6), fami_se, nompar(4), enr
    character(len=16) :: nomres(3)
    character(len=16), pointer :: compor(:) => null()
    aster_logical :: incr
!
!
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'RIGI', 'XINT', 'BID', 'RIGI', 'XINT'/
    data nomres/'E', 'NU', 'ALPHA'/
!
!
    call elref1(elrefp)
    call jevech('PTHETAR', 'L', ithet)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
    call tecael(iadzi, iazk24, noms=0)
!
!     SI LA VALEUR DE THETA EST NULLE SUR L'ÉLÉMENT, ON SORT
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
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                ibid, ibid, ibid, ddlm, nfiss, &
                contac)
!
!
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PMATERC', 'L', imate)
    call jevech('PGTHETA', 'E', igthet)
!
    incr = compor(4) .eq. 'COMP_INCR'

!
!
!     ------------------------------------------------------------------
!              CALCUL DE G SUR L'ELEMENT MASSIF
!     ------------------------------------------------------------------
!
!     PARAMÈTRES PROPRES À X-FEM
!
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PBASLOR', 'L', jbaslo)
    if (nfe .gt. 0) call jevech('PSTANO', 'L', jstno)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
    call teattr('S', 'XFEM', enr, ibid)
    if (enr(1:2) .eq. 'XH') call jevech('PHEA_NO', 'L', jheavn)
!
!     PROPRES AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    if (ibid .eq. 0 .and. ltequa(elrefp, enr)) call jevech('PPMILTO', 'L', jpmilt)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO, NPG ET IVF
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    fami_se = fami(ndim+irese)
    if (nfe .gt. 0) then
        if (ndim .eq. 3 .and. count(zi((jstno-1+1):(jstno-1+nnop)) .eq. -2) .eq. 0) fami_se = 'XGEO'
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami_se, nno=nno, npg=npg)
!
!     CALCUL DES FORCES NODALES CORRESPONDANT AUX CHARGES VOLUMIQUES
    call xcgfvo(option, ndim, nnop, fno)
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = zi(jlonch-1+1)
!
!       BOUCLE D'INTEGRATION SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       BOUCLE SUR LES SOMMETS DU SOUS-TRIA (DU SOUS-SEG)
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
        call xgelem(elrefp, ndim, coorse, igeom, jheavt, &
                    ise, nfh, ddlc, ddlm, nfe, &
                    zr(jbaslo), nnop, idepl, zr(jlsn), zr(jlst), &
                    igthet, fno, nfiss, jheavn, jstno, &
                    incr)
!
!
    end do
!
!     ------------------------------------------------------------------
!              CALCUL DE G SUR LES LEVRES
!     ------------------------------------------------------------------
!
    if (option .eq. 'CALC_G_XFEM') then
!       SI LA PRESSION N'EST CONNUE SUR AUCUN NOEUD, ON LA PREND=0.
        call jevecd('PPRESSR', ipres, 0.d0)
    else if (option .eq. 'CALC_G_XFEM_F') then
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
            call fointe('FM', zk8(ipref), 4, nompar, valpar, &
                        presn(i), icode)
        end do
    end if
!
!     SI LA VALEUR DE LA PRESSION EST NULLE SUR L'ÉLÉMENT, ON SORT
    compt = 0
    do i = 1, nnop
        if (option .eq. 'CALC_G_XFEM') pres = abs(zr(ipres-1+i))
        if (option .eq. 'CALC_G_XFEM_F') pres = abs(presn(i))
        if (pres .lt. r8prem()) compt = compt+1
    end do
    if (compt .eq. nnop) goto 999
!
!     PARAMETRES PROPRES A X-FEM
    call jevech('PPINTER', 'L', jptint)
    call jevech('PCFACE', 'L', jcface)
    call jevech('PLONGCO', 'L', jlongc)
    call jevech('PBASECO', 'L', jbasec)
!
!     RÉCUPÉRATIONS DES DONNÉES SUR LA TOPOLOGIE DES FACETTES
    ninter = zi(jlongc-1+1)
    nface = zi(jlongc-1+2)
    nptf = zi(jlongc-1+3)
    if (ninter .lt. ndim) goto 999
!
    do i = 1, nface
        do j = 1, nptf
            cface(i, j) = zi(jcface-1+ndim*(i-1)+j)
        end do
    end do
!
!     RECUPERATION DES DONNEES MATERIAU AU 1ER POINT DE GAUSS !!
!     LE MATÉRIAU DOIT ETRE HOMOGENE DANS TOUT L'ELEMENT
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
!     BOUCLE SUR LES FACETTES
!
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
!
!
end subroutine
