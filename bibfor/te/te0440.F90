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
subroutine te0440(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/teattr.h"
#include "asterfort/tefrep.h"
#include "asterfort/tecach.h"
#include "asterfort/xfovol.h"
#include "asterfort/xteddl.h"
#include "asterfort/xteini.h"
!
    character(len=16), intent(in) :: option, nomte
!......................................................................
!
!    - FONCTION REALISEE:  CALCUL DES VECTEURS ELEMENTAIRES
!                          OPTIONS  CHAR_MECA_FR3D3D ET CHAR_MECA_FF3D3D
!                                   CHAR_MECA_FR2D2D ET CHAR_MECA_FF2D2D
!                          POUR LES ÉLÉMENTS X-FEM
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!
!
    integer(kind=8) :: j, ndim, nno, nnop, nnops, nnos, nnom, nddl, npg, singu
    integer(kind=8) :: nfh, ddls, nfe, ddlc, nse, ise, in, ino, ibid, ddlm
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, jlsn, jlst, jstno
    integer(kind=8) :: ivectu, iforc, itemps, igeom, jpmilt, irese
    integer(kind=8) :: nfiss, jfisno, jheavn
    integer(kind=8) :: imate, jbaslo
    integer(kind=8) :: ncompn, heavn(27, 5), iret, jtab(7), ig
    real(kind=8) :: he, coorse(81)
    character(len=8) :: elrefp, elrese(6), fami(6), enr, lag
    aster_logical :: fonc, lbid
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'RIGI', 'XINT', 'BID', 'RIGI', 'XINT'/
!
!-----------------------------------------------------------------------
!
!
!     ELEMENT DE REFERENCE PARENT
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop, nnos=nnops)
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO, NPG
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, nnos=nnos, &
                     npg=npg)
!
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                nnom, ddls, nddl, ddlm, nfiss, &
                ibid)
!
!     PARAMETRE DU VECTEUR ELEMENTAIRE
!-------------------------------------
    call jevech('PVECTUR', 'E', ivectu)
!
!     PARAMÈTRES PROPRES À X-FEM
!-------------------------------
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PLSN', 'L', jlsn)
    call jevech('PLST', 'L', jlst)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PSTANO', 'L', jstno)
    if (nfe .gt. 0) then
        call jevech('PMATERC', 'L', imate)
        call jevech('PBASLOR', 'L', jbaslo)
    end if
    call teattr('S', 'XFEM', enr, ibid)
!     RECUPERATION DE LA DEFINITION DES DDL HEAVISIDES
    if (nfh .gt. 0) then
        call jevech('PHEA_NO', 'L', jheavn)
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
        ASSERT(ncompn .eq. 5)
        do ino = 1, nnop
            do ig = 1, ncompn
                heavn(ino, ig) = zi(jheavn-1+ncompn*(ino-1)+ig)
            end do
        end do
    end if
!     PROPRE AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    if ((ibid .eq. 0) .and. (enr .eq. 'XH' .or. enr .eq. 'XHT' &
                             .or. enr .eq. 'XT' .or. enr .eq. 'XHC') &
        .and. .not. iselli(elrefp)) call jevech('PPMILTO', 'L', jpmilt)
    if (nfiss .gt. 1) call jevech('PFISNO', 'L', jfisno)
!
!     PARAMÈTRES DES FORCES VOLUMIQUES
!-------------------------------------
!
    if ((option .eq. 'CHAR_MECA_FR3D3D') .or. (option .eq. 'CHAR_MECA_FR2D2D')) then
!
        fonc = .false.
        if (ndim .eq. 3) call tefrep(option, 'PFR3D3D', iforc)
        if (ndim .eq. 2) call tefrep(option, 'PFR2D2D', iforc)
!
    elseif ((option .eq. 'CHAR_MECA_FF3D3D') .or. (option .eq. 'CHAR_MECA_FF2D2D')) then
!
        fonc = .true.
        call jevech('PINSTR', 'L', itemps)
        if (ndim .eq. 3) call jevech('PFF3D3D', 'L', iforc)
        if (ndim .eq. 2) call jevech('PFF2D2D', 'L', iforc)
!
    end if
!
!     RÉCUPÉRATION DE LA SUBDIVISION DE L'ÉLÉMENT EN NSE SOUS ELEMENT
    nse = zi(jlonch-1+1)
!
!       BOUCLE SUR LES NSE SOUS-ELEMENTS
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
!       FONCTION HEAVYSIDE CSTE SUR LE SS-ÉLT
        he = zi(jheavt-1+ise)
!
        call xfovol(elrefp, ndim, coorse, igeom, he, &
                    nfh*ndim, ddlc, singu, nnop, jlsn, &
                    jlst, heavn, iforc, itemps, ivectu, fonc, &
                    .true._1, zi(imate), jbaslo, jstno)
!
    end do
!
!     SUPPRESSION DES DDLS SUPERFLUS
    call teattr('C', 'XLAG', lag, ibid)
    if (ibid .eq. 0 .and. lag .eq. 'ARETE') then
        nnop = nnos
    end if
    call xteddl(ndim, nfh, nfe, ddls, nddl, &
                nnop, nnops, zi(jstno), .false._1, lbid, &
                option, nomte, ddlm, nfiss, jfisno, &
                vect=zr(ivectu))
!
end subroutine
