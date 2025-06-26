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
subroutine te0408(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dxtpif.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/fonbpa.h"
#include "asterfort/get_value_mode_local.h"
#include "asterfort/jevech.h"
#include "asterfort/jeveuo.h"
#include "asterfort/provec.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!
! ------------------------------------------------------------------------------
!
    integer(kind=8) :: itabp(8), itempp, ino, nbcou, npgh, itemps
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, ivf, idfdx, jgano, isp
    integer(kind=8) :: ier, igauh, icou, jnbspi, iret, itempf, jresu, jgeom
    integer(kind=8) :: pta, ptb, ptc, ptd
!
    real(kind=8) :: tpinf, tpmoy, tpsup, cp1, cp2, cp3, tpc, zic, zmin
    real(kind=8) :: hic, epais, excent, norm2
!
    character(len=16) :: option, nomte
!
    aster_logical :: tempno, elem_dkt, CasIsOk
!
    integer(kind=8), parameter :: MaxPara = 4
    integer(kind=8) :: jprof, nbpara, casfct
    real(kind=8) :: valpu(MaxPara), xyz(3), cdg(3), vect(3), vectab(3), vectcd(3)
    character(len=8) :: nompu(MaxPara)
    character(len=16) :: nompara(MaxPara)
    character(len=19) :: nomfon
    character(len=24) :: chprol
!
    integer(kind=8) :: retpara(2)
    real(kind=8) :: valr(2)
    character(len=8) :: valp(2), k8bid
    blas_int :: b_incx, b_incy, b_n
!
! ------------------------------------------------------------------------------
!
    if (option .ne. 'PREP_VRC') then
        ASSERT(.false.)
    end if
    !
! Les éléments traités
!   dkt         :   MEDKQU4     MEDKTR3
!   dst         :   MEDSQU4     MEDSTR3
!   q4g         :   MEQ4QU4     MET3TR3
!   coque_3d    :   MEC3QU9H    MEC3TR7H
!   coque_axis  :   MECXSE3
    !
    elem_dkt = (nomte .eq. 'MEDKQU4') .or. (nomte .eq. 'MEDKTR3') .or. (nomte .eq. 'MEDSQU4') &
               .or. (nomte .eq. 'MEDSTR3') .or. (nomte .eq. 'MEQ4QU4') .or. &
               (nomte .eq. 'MET3TR3') .or. (nomte .eq. 'MEC3QU9H') .or. (nomte .eq. 'MEC3TR7H') &
               .or. (nomte .eq. 'MECXSE3')
    !
    ASSERT(elem_dkt)
    !
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    !
    call jevech('PTEMPCR', 'E', jresu)
! NBCOU : nombre de couches
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi)
    !
    epais = 0.0; excent = 0.0
    valp(1:2) = [character(8) :: 'EP', 'EXCENT']
    call get_value_mode_local('PCACOQU', valp, valr, iret, retpara_=retpara, &
                              nbpara_=2)
    if (retpara(1) .eq. 0) then
        epais = valr(1)
    else
        ASSERT(.FALSE.)
    end if
    if (retpara(2) .eq. 0) then
        excent = valr(2)
    else
        excent = 0.0
    end if
    !
! Si la température est aux noeuds
    call tecach('ONN', 'PTEMPER', 'L', iret, nval=8, &
                itab=itabp)
    if (iret .eq. 0 .or. iret .eq. 3) then
        tempno = .TRUE.
        itempp = itabp(1)
        !
        cdg(1:3) = 0.0
! Calcul des températures ELEM sur les couches INF, SUP, MOY
        tpinf = 0.d0; tpmoy = 0.d0; tpsup = 0.d0
        do ino = 1, nno
            call dxtpif(zr(itempp+3*(ino-1)), zl(itabp(8)+3*(ino-1)))
            tpmoy = tpmoy+zr(itempp-1+3*(ino-1)+1)/dble(nno)
            tpinf = tpinf+zr(itempp-1+3*(ino-1)+2)/dble(nno)
            tpsup = tpsup+zr(itempp-1+3*(ino-1)+3)/dble(nno)
        end do
! coefficients des polynomes de degré 2
        cp1 = tpmoy
        cp2 = (tpsup-tpinf)/epais
        cp3 = 2.d0*(tpinf+tpsup-2.d0*tpmoy)/(epais*epais)
! Dans ce cas pas de prise en compte de l'excentrement
        excent = 0.0
    else
! Si la température est issu d'un champ de fonctions
        call tecach('ONO', 'PTEMPEF', 'L', iret, iad=itempf)
        ASSERT(iret .eq. 0)
! Les paramètres de la fonction
        nomfon = zk8(itempf)
        chprol = nomfon//'.PROL'
        call jeveuo(chprol, 'L', jprof)
        call fonbpa(nomfon, zk24(jprof), k8bid, MaxPara, nbpara, &
                    nompara)
! C'est soit :
!   INST    EPAIS
!   INST    EXCEN
!   INST    X   Y   Z
        !
!   INST  EPAIS EXCEN  X   Y   Z   !!! Vive les entiers codés
!   1     2     4      8  16  32
        casfct = 0
        do ino = 1, nbpara
            select case (nompara(ino))
            case ('INST')
                casfct = casfct+1
            case ('EPAIS')
                casfct = casfct+2
            case ('EXCENT')
                casfct = casfct+4
            case ('X')
                casfct = casfct+8
            case ('Y')
                casfct = casfct+16
            case ('Z')
                casfct = casfct+32
            end select
        end do
        CasIsOk = (casfct .eq. (1+2)) .or. (casfct .eq. (1+4)) .or. (casfct .eq. (1+8+16+32))
        if (.not. CasIsOk) then
            call utmess('F', 'FONCT0_80', sk=nomfon)
        end if
        !
        tempno = .FALSE.
        call jevech('PINST_R', 'L', itemps)
        nompu(1) = 'INST'; valpu(1) = zr(itemps)
        !
! On va chercher la géométrie
        call jevech('PGEOMER', 'L', jgeom)
! Calcul du cdg et de la normale
        if (nnos .eq. 3) then
            pta = jgeom; ptb = jgeom+6; ptc = jgeom+3; ptd = jgeom
            cdg(1:3) = (zr(pta:pta+2)+zr(ptb:ptb+2)+zr(ptc:ptc+2))/3.0
        else if (nnos .eq. 4) then
            pta = jgeom; ptb = jgeom+6; ptc = jgeom+3; ptd = jgeom+9
            cdg(1:3) = (zr(pta:pta+2)+zr(ptb:ptb+2)+zr(ptc:ptc+2)+zr(ptd:ptd+2))/4.0
        else
            ASSERT(.FALSE.)
        end if
        vectab(1:3) = zr(ptb:ptb+2)-zr(pta:pta+2)
        vectcd(1:3) = zr(ptd:ptd+2)-zr(ptc:ptc+2)
        call provec(vectab, vectcd, vect)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norm2 = ddot(b_n, vect, b_incx, vect, b_incy)
        vect = vect(1:3)/sqrt(norm2)
    end if
    !
! COQUES :
!     En thermique les coques ne sont pas excentrées
!     En mécanique elles peuvent être excentrées
!     Si calcul thermique puis mécanique ( tempno )
!         Les températures INF, SUP, MOY ne doivent pas tenir compte de EXCENT
!     Si champ de fonction thermique puis calcul mécanique
!         EPAIS  c'est SANS la prise en compte de EXCENT
!         EXCENT c'est AVEC la prise en compte de EXCENT
    !
! hic    : épaisseur d'une couche
! npgh   : nombre de points par couche
    !
    hic = epais/nbcou
    npgh = 3
    !
! CALCUL DE LA TEMPERATURE SUR LES COUCHES
    zmin = -epais/2.d0
    !
    do icou = 1, nbcou
        do igauh = 1, npgh
            isp = (icou-1)*npgh+igauh
! zic : [-epais/2 ; epais/2 ]
            if (igauh .eq. 1) then
                zic = zmin+(icou-1)*hic
            else if (igauh .eq. 2) then
                zic = zmin+(icou-1)*hic+hic/2.0
            else
                zic = zmin+(icou-1)*hic+hic
            end if
            !
            if (tempno) then
                tpc = cp3*zic*zic+cp2*zic+cp1
            else
                if (casfct .eq. (1+2)) then
! Avec EPAIS donc SANS excentrement
                    nbpara = 2
                    nompu(2) = 'EPAIS'; valpu(2) = zic
                else if (casfct .eq. (1+4)) then
! Avec EXCEN donc AVEC excentrement
                    nbpara = 2
                    nompu(2) = 'EXCENT'; valpu(2) = excent+zic
                else if (casfct .eq. (1+8+16+32)) then
! Avec X Y Z
                    nbpara = 4
                    xyz(1:3) = cdg(1:3)+vect(1:3)*(excent+zic)
                    nompu(2) = 'X'; valpu(2) = xyz(1)
                    nompu(3) = 'Y'; valpu(3) = xyz(2)
                    nompu(4) = 'Z'; valpu(4) = xyz(3)
                end if
                call fointe(' ', zk8(itempf), nbpara, nompu, valpu, &
                            tpc, ier)
                if (ier .ne. 0) then
                    call utmess('F', 'FONCT0_80', sk=zk8(itempf))
                end if
            end if
            !
            zr(jresu-1+isp) = tpc
        end do
    end do
    !
end subroutine
