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
!
subroutine te0354(option, nomte)
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "blas/daxpy.h"
#include "blas/ddot.h"
    character(len=16) :: option, nomte
!
! ----------------------------------------------------------------------
!                   SOURCE THERMIQUE NON LINEAIRE
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'CHAR_THER_SOURNL'
!                                   'RESI_THER_SOURNL'
!                                   'MTAN_THER_SOURNL'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    integer :: nnomax
    parameter (nnomax=27)
!
    aster_logical :: axi, resi
!      INTEGER NDIM,NNO,NPG,NNOS,G,I,OS,OSM,M,N,IW,IVF,IDFDE,IRET,JGANO
    integer :: ndim, nno, npg, nnos, g, os, osm, m, n, iw, ivf, idfde, iret
    integer :: jgano
    integer :: igeom, itemps, ivect, imatr, isour, ither
    real(kind=8) :: dfdx(nnomax), dfdy(nnomax), dfdz(nnomax), w, rg, theta, sour
    real(kind=8) :: tg
    real(kind=8) :: dsdt, coef, coefop
!      CHARACTER*8 ELREFE,FCT
    character(len=8) :: elrefe
!
!
!    LECTURE DES PARAMETRES COMMUNS
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PSOURNL', 'L', isour)
    if (zk8(isour)(1:7) .eq. '&FOZERO') goto 999
    theta = zr(itemps+2)
!
!    LECTURE DES PARAMETRES SPECIFIQUES A CHAQUE OPTION
    if (option(1:4) .eq. 'CHAR') then
        resi = .true.
        coefop = 1-theta
        call jevech('PTEMPER', 'L', ither)
        call jevech('PVECTTR', 'E', ivect)
    else if (option(1:4).eq.'RESI') then
        resi = .true.
        coefop = -theta
        call jevech('PTEMPEI', 'L', ither)
        call jevech('PRESIDU', 'E', ivect)
    else
        resi = .false.
        coefop = -theta
        call jevech('PTEMPEI', 'L', ither)
        call jevech('PMATTTR', 'E', imatr)
    endif
!
!    ACCES AUX CARACTERISTIQUES DE L'ELEMENT FINI
    call elref1(elrefe)
    call elrefe_info(elrefe=elrefe, fami='RIGI', ndim=ndim, nno=nno, nnos=nnos,&
                     npg=npg, jpoids=iw, jvf=ivf, jdfde=idfde, jgano=jgano)
    axi = lteatt('AXIS','OUI')
!
    do g = 1, npg
        os = (g-1)*nno
!
!      CALCUL DU POIDS DU POINT DE GAUSS
        if (ndim .eq. 2) then
            call dfdm2d(nno, g, iw, idfde, zr(igeom),&
                        w, dfdx, dfdy)
            if (axi) then
                rg = ddot(nno,zr(igeom),2,zr(ivf+os),1)
                w = w*rg
            endif
        else
            call dfdm3d(nno, g, iw, idfde, zr(igeom),&
                        w, dfdx, dfdy, dfdz)
        endif
!
!      CALCUL DE LA TEMPERATURE AU POINT DE GAUSS
        tg = ddot(nno,zr(ither),1,zr(ivf+os),1)
!
!      CALCUL DU RESIDU
        if (resi) then
!
!        CALCUL DE LA SOURCE
            call fointe('FM', zk8(isour), 1, ['TEMP'], [tg],&
                        sour, iret)
            coef = w*sour*coefop
!
!        CONTRIBUTION AU RESIDU
            call daxpy(nno, coef, zr(ivf+os), 1, zr(ivect),&
                       1)
!
!      CALCUL DE LA MATRICE TANGENTE (STOCKAGE SYMETRIQUE)
        else
!
!        CALCUL DE LA DERIVEE DE LA SOURCE PAR RAPPORT A LA TEMPERATURE
            call foderi(zk8(isour), tg, sour, dsdt)
            coef = w*dsdt*coefop
!
!        CONTRIBUTION A LA MATRICE
            osm = 0
            do n = 0, nno-1
                do m = 0, n
                    zr(imatr+osm)=zr(imatr+osm)+coef*zr(ivf+os+n)*zr(&
                    ivf+os+m)
                    osm = osm + 1
                end do
            end do
        endif
    end do
!
999 continue
end subroutine
