! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine te0050(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/pmfmats.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!    - fonction réalisée :  calcul des matrices élémentaires
!                          option :'RIGI_MECA_HYST'
!        pour tous les types d'éléments (sauf les éléments discrets)
!
!    - arguments:
!        données:      option       -->  option de calcul
!                      nomte        -->  nom du type élément
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nbres, nbpar
    parameter  ( nbres=2 )
    parameter  ( nbpar=3 )
!
    integer :: jgano, iret, nbval, idimge, npara
    integer :: i, k, mater, irigi
    integer :: iresu, imate
    integer :: idresu(5), idrigi(2), idgeo(5)
    integer :: ipoids, ivf, idfdx, igeom
    integer :: ndim, nno, nnos, npg1, ino
!
    real(kind=8) :: eta, valres(nbres), valpar(nbpar), vxyz
!
    integer :: icodre(nbres)
    character(len=8) :: nompar(nbpar), nomat
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'RIGI_MECA_HYST')
!
    call elrefe_info(fami='RIGI',ndim=ndim,nno=nno,nnos=nnos,&
                     npg=npg1,jpoids=ipoids,jvf=ivf,jdfde=idfdx,jgano=jgano)
!
!   récupération des champs paramètres et de leurs longueurs:
    call tecach('ONO', 'PMATUUC', 'E', iret, nval=5, itab=idresu)
!
    nbval= idresu(2)
!
    nompar(1)='X'
    nompar(2)='Y'
    nompar(3)='Z'
!
    call tecach('ONO', 'PGEOMER', 'L', iret, nval=5, itab=idgeo)
    igeom=idgeo(1)
    idimge=idgeo(2)/nno
!
    ASSERT(idimge.eq.2 .or. idimge.eq.3)
!
    npara=idimge
    do k = 1, idimge
        vxyz = 0.d0
        do ino = 1, nno
            vxyz = vxyz+zr(igeom + idimge*(ino-1) +k -1)
        enddo
        valpar(k) = vxyz/nno
    enddo
!
    call jevech('PMATERC', 'L', imate)
    mater=zi(imate)
    call rccoma(mater, 'ELAS', 0, phenom, icodre(1))
    if(.not.(phenom .eq. 'ELAS'       .or. phenom .eq. 'ELAS_COQMU' .or. &
             phenom .eq. 'ELAS_GLRC'  .or. phenom .eq. 'ELAS_DHRC'  .or. &
             phenom .eq. 'ELAS_ORTH')) then
        call utmess('F', 'PLATE1_1', nk=2, valk=[option, phenom])
    endif

!   si l'élément est multifibre, il faut prendre le materiau "section"
!   pour récupérer les coefficients de dilatation :
    call pmfmats(mater, nomat)
!
    call tecach('ONO', 'PRIGIEL', 'L', iret, nval=2, itab=idrigi)
    ASSERT(idrigi(2).eq.nbval)
!
!   récupération des coefficients fonctions de la géométrie :
    nomres(1)='AMOR_HYST'
    valres(1) = 0.d0
    call rcvalb('RIGI', 1, 1, '+', mater, nomat, phenom, npara, nompar, valpar, 1,&
                nomres, valres, icodre, 0, nan='NON')
!
!   calcul proprement dit
    iresu= idresu(1)
    irigi= idrigi(1)
    eta = valres(1)
    do i = 1, nbval
        zc(iresu-1+i)=dcmplx(zr(irigi-1+i),eta*zr(irigi-1+i))
    enddo
!
end subroutine
