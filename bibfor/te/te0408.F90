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
subroutine te0408(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dxtpif.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "asterfort/tecach.h"
!
    character(len=16) :: option, nomte
!
!     CHARACTER*32 JEXNUM,JEXNOM,JEXR8,JEXATR
    integer :: itabp(8), itempp, ino, nbcou, npgh, itemps, ibid
    integer :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano, isp
    integer :: ier, igauh, icou, jnbspi, iret, itempf, jresu, icacoq
    real(kind=8) :: tpinf, tpmoy, tpsup, cp1, cp2, cp3, tpc, zic, zmin
    real(kind=8) :: inst, valpu(2), hic, h
    aster_logical :: tempno
    character(len=8) :: nompu(2), alias8
!
    call teattr('S', 'ALIAS8', alias8, ibid)
    ASSERT(.not.lteatt('GRILLE','OUI'))
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    call jevech('PTEMPCR', 'E', jresu)
    call jevech('PNBSP_I', 'L', jnbspi)
!     NBCOU : NOMBRE DE COUCHES
    nbcou=zi(jnbspi-1+1)
    call jevech('PCACOQU', 'L', icacoq)
    h=zr(icacoq)
    hic=h/nbcou
!
!
!     1- SI LA TEMPERATURE EST AUX NOEUDS (TEMP/TEMP_SUP/TEMP_INF):
!        ------------------------------------------------------------
    call tecach('ONN', 'PTEMPER', 'L', iret, nval=8,&
                itab=itabp)
    if (iret .eq. 0 .or. iret .eq. 3) then
        tempno=.true.
        itempp=itabp(1)
!       -- CALCUL DES TEMPERATURES INF, SUP ET MOY
!          (MOYENNE DES NNO NOEUDS) ET DES COEF. DES POLY. DE DEGRE 2 :
!          ------------------------------------------------------------
        tpinf=0.d0
        tpmoy=0.d0
        tpsup=0.d0
        do ino = 1, nno
            call dxtpif(zr(itempp+3*(ino-1)), zl(itabp(8)+3*(ino-1)))
            tpmoy=tpmoy+zr(itempp-1+3*(ino-1)+1)/dble(nno)
            tpinf=tpinf+zr(itempp-1+3*(ino-1)+2)/dble(nno)
            tpsup=tpsup+zr(itempp-1+3*(ino-1)+3)/dble(nno)
        end do
!
        cp1=tpmoy
        cp2=(tpsup-tpinf)/h
        cp3=2.d0*(tpinf+tpsup-2.d0*tpmoy)/(h*h)
!
!
    else
!     2- SI LA TEMPERATURE EST UNE FONCTION DE 'INST' ET 'EPAIS'
!        -------------------------------------------------------
        call tecach('ONO', 'PTEMPEF', 'L', iret, iad=itempf)
        ASSERT(iret.eq.0)
        call jevech('PINST_R', 'L', itemps)
        inst=zr(itemps)
        tempno=.false.
        nompu(1)='INST'
        nompu(2)='EPAIS'
    endif
!
!
!     -- CALCUL DE LA TEMPERATURE SUR LES COUCHES :
!     ----------------------------------------------
!     NPGH  : NOMBRE DE POINTS PAR COUCHE
    npgh=3
!
!
    zmin=-h/2.d0
!
!
    do icou = 1, nbcou
        do igauh = 1, npgh
            isp=(icou-1)*npgh+igauh
!
            if (igauh .eq. 1) then
                zic=zmin+(icou-1)*hic
            else if (igauh.eq.2) then
                zic=zmin+hic/2.d0+(icou-1)*hic
            else
                zic=zmin+hic+(icou-1)*hic
            endif
!
            if (tempno) then
                tpc=cp3*zic*zic+cp2*zic+cp1
            else
                valpu(2)=zic
                valpu(1)=inst
                call fointe('FM', zk8(itempf), 2, nompu, valpu,&
                            tpc, ier)
            endif
!
            zr(jresu-1+isp)=tpc
        end do
    end do
!
end subroutine
