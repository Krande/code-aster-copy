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

subroutine te0481(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iselli.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/ltequa.h"
    character(len=16) :: option, nomte
!
! person_in_charge: patrick.massin at edf.fr
!
!.......................................................................
!
!      BUT : CALCULATION OF GAUSS POINTS COORDINATES FOR X-FEM ELEMENTS 1D 2D 3D
!            BUT NOT BOUNDARY ELEMENTS
!
!      OPTION : COOR_ELGA
!
!      INPUT : OPTION (CALCULATION OPTION)
!              NOMTE : ELEMENT TYPE NAME
!.......................................................................
!
    character(len=8) :: elrefp, elrese(6), fami(6), enr
    real(kind=8) :: xg(3), coorse(81), jac, r
    integer(kind=8) :: ibid, ndim, nnop, nno, npg, ivf, ipoids
    integer(kind=8) :: jpmilt, irese
    integer(kind=8) :: jpintt, jcnset, jlonch, igeom, jcopg
    integer(kind=8) :: i, j, nse, ise, in, ino, ipg, kpg, idfde
    aster_logical :: axi
!
    data elrese/'SE2', 'TR3', 'TE4', 'SE3', 'TR6', 'T10'/
    data fami/'BID', 'XINT', 'XINT', 'BID', 'XINT', 'XINT'/
!
!
!.......................................................................
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
    call elrefe_info(elrefe=elrese(ndim+irese), fami=fami(ndim+irese), &
                     nno=nno, npg=npg, jvf=ivf, jdfde=idfde, jpoids=ipoids)
!
!-----------------------------------------------------------------------
!     RECUPERATION DES ENTREES / SORTIE
!-----------------------------------------------------------------------
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PPINTTO', 'L', jpintt)
    call jevech('PCNSETO', 'L', jcnset)
    call jevech('PLONCHA', 'L', jlonch)
    call jevech('PCOORPG', 'E', jcopg)
!     PROPRES AUX ELEMENTS 1D ET 2D (QUADRATIQUES)
    call teattr('S', 'XFEM', enr, ibid)
    if ((ibid .eq. 0) .and. ltequa(elrefp, enr)) &
        call jevech('PPMILTO', 'L', jpmilt)
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

!           CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!           AVEC LES COORDONNEES DU SOUS-ELEMENT
            if (ndim .eq. 2) then
                call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                            jac)
            else if (ndim .eq. 3) then
                call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                            jac)
            end if
!
! -         CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE):
            if (axi) then
                r = xg(1)
!
                ASSERT(r .gt. 0d0)
!               MODIFICATION DU JACOBIEN
                jac = jac*r
            end if

!         NUMERO DE CE POINT DE GAUSS DANS LA FAMILLE 'XFEM'
            ipg = (ise-1)*npg+kpg
!
            do i = 1, ndim
                zr(jcopg+(ndim+1)*(ipg-1)-1+i) = xg(i)
            end do
            zr(jcopg+(ndim+1)*(ipg-1)-1+ndim+1) = jac
!
        end do
!
!-----------------------------------------------------------------------
!         FIN DE LA BOUCLE SUR LES POINTS DE GAUSS DU SOUS-ELT
!-----------------------------------------------------------------------
!
    end do
!
end subroutine
