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
subroutine te0046(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/ltequa.h"
#include "asterfort/reeref.h"
#include "asterfort/teattr.h"
#include "asterfort/xteini.h"
    character(len=16) :: option, nomte
!
! person_in_charge: samuel.geniaut at edf.fr
!
!.......................................................................
!
!     BUT: CALCUL DES COORDONNEES DES POINTS DE GAUSS
!          DE LA FAMILLE X-FEM (POINTS DE GAUSS DES SOUS-ÉLÉMENTS)
!          DANS L'ESPACE DE L'ELEMENT PARENT DE REFERENCE
!
!          OPTIONS : 'XFEM_XPG'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!.......................................................................
!
!
!      CHARACTER*8   ELREF,FPG,ELC,NOMPAR(4)
!      INTEGER NDIM,NNO,NNOS,NPG,IPOIDS,IVF,IDFDE,JGANO
!      INTEGER NFH,NFE,SINGU,DDLC,DDLS,NDDL
!      INTEGER IGEOM,IPRES,ITEMPS,IFORC,IRET,IRES
!      INTEGER JLST,JPTINT,JAINT,JCFACE,JLONCH
!      INTEGER I,J,NINTER,NFACE,CFACE(5,3),IFA,NLI,IN(3),IG
!      INTEGER AR(12,3),NBAR,FAC(6,4),NBF,IBID2(12,3),IBID,INO,ILEV
!      INTEGER NNOF,NPGF,IPOIDF,IVFF,IDFDEF,IPGF,POS
!      REAL*8  MULT,PRES,CISA, FORREP(3,2),FF(27),G(3),JAC,ND(3),HE(2)
!      REAL*8  RR(2),LST,XG(4)
!      DATA    HE / -1.D0 , 1.D0/
!
    character(len=8) :: elrefp, elrese(6), fami(6), enr
    real(kind=8) :: xg(3), xe(3), ff(27), coorse(81)
    integer(kind=8) :: ibid, ndim, nnop, nno, npg, ivf
    integer(kind=8) :: nfh, nfe, singu, ddlc, jpmilt, irese
    integer(kind=8) :: jpintt, jcnset, jheavt, jlonch, igeom, jout
    integer(kind=8) :: i, j, nse, ise, in, ino, ipg, kpg
    aster_logical :: axi
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!
!-----------------------------------------------------------------------
!     INITIALISATIONS
!-----------------------------------------------------------------------
!
!     ELEMENT DE REFERENCE PARENT : RECUP DE NDIM ET NNOP
    call elref1(elrefp)
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nnop)
!
    axi = lteatt('AXIS', 'OUI')
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO, NPG ET IVF
    if (.not. iselli(elrefp)) then
        irese = 3
    else
        irese = 0
    end if
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), nno=nno, npg=npg, jvf=ivf)
!
!     INITIALISATION DES DIMENSIONS DES DDLS X-FEM
    call xteini(nomte, nfh, nfe, singu, ddlc, &
                ibid, ibid, ibid, ibid, ibid, &
                ibid)
!
!-----------------------------------------------------------------------
!     RECUPERATION DES ENTREES / SORTIE
!-----------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PHEAVTO', 'L', jheavt)
    call jevech('PLONCHA', 'L', jlonch)
!     PROPRES AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    call teattr('S', 'XFEM', enr, ibid)
    if ((ibid .eq. 0) .and. ltequa(elrefp, enr)) call jevech('PPMILTO', 'L', jpmilt)
!
    call jevech('PXFGEOM', 'E', jout)
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
!
!-----------------------------------------------------------------------
!         BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
!
        do kpg = 1, npg
!
!         COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
            xg(:) = 0.d0
            do i = 1, ndim
                do in = 1, nno
                    xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+in)*coorse(ndim*(in-1)+i)
                end do
            end do
!
!         COORDONNEES DU PG DANS L'ELEMENT DE REF PARENT : XE
            call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                        xe, ff)
!
!         NUMERO DE CE POINT DE GAUSS DANS LA FAMILLE 'XFEM'
            ipg = (ise-1)*npg+kpg
!
            do j = 1, ndim
                zr(jout-1+ndim*(ipg-1)+j) = xe(j)
            end do
!
!
        end do
!
!-----------------------------------------------------------------------
!         FIN DE LA BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
!
!
    end do
!
!
end subroutine
