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
subroutine te0060(option, nomte)
!     BUT: CALCUL DES VECTEURS ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT AU TERME D'ECHANGE (FONCTION)
!          SUR DES FACES D'ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'CHAR_THER_TEXT_F/R'
!                   'CHAR_THER_RAYO_F/R'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    character(len=8) :: nompar(4)
    character(len=16) :: nomte, option
    real(kind=8) :: nx, ny, nz, sx(9, 9), sy(9, 9), sz(9, 9), jac, theta
    real(kind=8) :: valpar(4), xx, yy, zz, tem, echn, echnp1, texn, texnp1
    real(kind=8) :: sigm1, sigmn, eps1, epsn, tpf1, tpfn, tz0
    integer :: ipoids, ivf, idfdx, idfdy, igeom, itemps, i, ndim, nno, ipg, npg1
    integer :: ivectt, itext, iech, iray, ino, idec, jdec, kdec, ldec, itemp
    integer :: jno, j, ier, jgano, nnos
    aster_logical :: ltext
!
!
!====
! 1.1 PREALABLES: RECUPERATION ADRESSES FONCTIONS DE FORMES...
!====
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfdx, jgano=jgano)
    idfdy = idfdx + 1
!
    tz0 = r8t0()
!
! INITS.
    if (option(11:14) .eq. 'TEXT') then
        ltext = .true.
    else if (option(11:14).eq.'RAYO') then
        ltext = .false.
    else
!C OPTION DE CALCUL INVALIDE
        ASSERT(.false.)
    endif
!====
! 1.2 PREALABLES LIES AUX RECHERCHES DE DONNEES GENERALES
!====
    if (ltext) then
! CHAR_.._TEXT : 2 TYPES DE CALCUL
!
        call jevech('PTEMPER', 'L', itemp)
        call jevech('PCOEFHF', 'L', iech)
        call jevech('PT_EXTF', 'L', itext)
!
    else
! CHAR_..._RAYO: 4 TYPES DE CALCUL
!
! CHAMP DE RAYONNEMENT
        call jevech('PRAYONF', 'L', iray)
        call jevech('PTEMPER', 'L', itemp)
! FIN DU IF LTEXT
    endif
!
! TRONC COMMUN
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PTEMPSR', 'L', itemps)
    call jevech('PVECTTR', 'E', ivectt)
!
!====
! 1.3 PREALABLES LIES AUX CALCULS
!====
    theta = zr(itemps+2)
    nompar(1) = 'X'
    nompar(2) = 'Y'
    nompar(3) = 'Z'
    nompar(4) = 'INST'
    do i = 1, nno
        zr(ivectt+i-1) = 0.0d0
    end do
!
!    CALCUL DES PRODUITS VECTORIELS OMI X OMJ
!
    do ino = 1, nno
        i = igeom + 3*(ino-1) -1
        do jno = 1, nno
            j = igeom + 3*(jno-1) -1
            sx(ino,jno) = zr(i+2) * zr(j+3) - zr(i+3) * zr(j+2)
            sy(ino,jno) = zr(i+3) * zr(j+1) - zr(i+1) * zr(j+3)
            sz(ino,jno) = zr(i+1) * zr(j+2) - zr(i+2) * zr(j+1)
        end do
    end do
!
!====
! 2. CALCULS TERMES DE MASSE
!====
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do ipg = 1, npg1
        kdec = (ipg-1)*nno*ndim
        ldec = (ipg-1)*nno
        nx = 0.0d0
        ny = 0.0d0
        nz = 0.0d0
!
!    CALCUL DE LA NORMALE AU POINT DE GAUSS IPG
        do i = 1, nno
            idec = (i-1)*ndim
            do j = 1, nno
                jdec = (j-1)*ndim
                nx=nx+zr(idfdx+kdec+idec) * zr(idfdy+kdec+jdec) * sx(&
                i,j)
                ny=ny+zr(idfdx+kdec+idec) * zr(idfdy+kdec+jdec) * sy(&
                i,j)
                nz=nz+zr(idfdx+kdec+idec) * zr(idfdy+kdec+jdec) * sz(&
                i,j)
            end do
        end do
        jac = sqrt(nx*nx + ny*ny + nz*nz)
        tem = 0.d0
        xx = 0.d0
        yy = 0.d0
        zz = 0.d0
!
        if (itemp .ne. 0) then
            do i = 1, nno
! CALCUL DE T-
                tem = tem + zr(itemp+i-1) * zr(ivf+ldec+i-1)
            end do
        endif
!
        do i = 1, nno
! CALCUL DE LA POSITION GEOMETRIQUE DU PT DE GAUSS
            xx = xx + zr(igeom+3*i-3) * zr(ivf+ldec+i-1)
            yy = yy + zr(igeom+3*i-2) * zr(ivf+ldec+i-1)
            zz = zz + zr(igeom+3*i-1) * zr(ivf+ldec+i-1)
        end do
!
        valpar(1) = xx
        valpar(2) = yy
        valpar(3) = zz
        valpar(4) = zr(itemps)
!
!====
! 2.1 OPTION CHAR_THER_TEXT_F/R
!====
        if (ltext) then
!
            call fointe('FM', zk8(itext), 4, nompar, valpar,&
                        texnp1, ier)
!
            if (theta .ne. 1.0d0) then
                valpar(4) = zr(itemps)-zr(itemps+1)
                call fointe('FM', zk8(itext), 4, nompar, valpar,&
                            texn, ier)
!
            else
                texn = 0.d0
            endif
!
            valpar(4) = zr(itemps)
            call fointe('FM', zk8(iech), 4, nompar, valpar,&
                        echnp1, ier)
!
            if (theta .ne. 1.0d0) then
                valpar(4) = zr(itemps)-zr(itemps+1)
                call fointe('FM', zk8(iech), 4, nompar, valpar,&
                            echn, ier)
!
            else
                echn = 0.d0
            endif
!
            do i = 1, nno
                zr(ivectt+i-1) = zr(ivectt+i-1) + jac * zr(ipoids+ipg- 1) * zr(ivf+ldec+i-1) * ( &
                                 &theta*echnp1*texnp1+(1.0d0- theta)*echn*(texn-tem))
            end do
!
!====
! 2.2 OPTION CHAR_THER_RAYO_F/R
!====
        else
!
            call fointe('FM', zk8(iray), 4, nompar, valpar,&
                        sigm1, ier)
            if (theta .ne. 1.0d0) then
                valpar(4) = zr(itemps)-zr(itemps+1)
                call fointe('FM', zk8(iray), 4, nompar, valpar,&
                            sigmn, ier)
            else
                sigmn = 0.d0
            endif
!
            valpar(4) = zr(itemps)
            call fointe('FM', zk8(iray+1), 4, nompar, valpar,&
                        eps1, ier)
            if (theta .ne. 1.0d0) then
                valpar(4) = zr(itemps)-zr(itemps+1)
                call fointe('FM', zk8(iray+1), 4, nompar, valpar,&
                            epsn, ier)
            else
                epsn = 0.d0
            endif
!
            valpar(4) = zr(itemps)
            call fointe('FM', zk8(iray+2), 4, nompar, valpar,&
                        tpf1, ier)
            if (theta .ne. 1.0d0) then
                valpar(4) = zr(itemps)-zr(itemps+1)
                call fointe('FM', zk8(iray+2), 4, nompar, valpar,&
                            tpfn, ier)
            else
                tpfn = 0.d0
            endif
!
            do i = 1, nno
                zr(ivectt+i-1) = zr(ivectt+i-1) + jac * zr(ipoids+ipg- 1) * zr(ivf+ldec+i-1) * ( &
                                 &theta *sigm1*eps1* (tpf1+ tz0)**4+ (1.0d0-theta)*sigmn*epsn*((t&
                                 &pfn+tz0)**4-(tem+ tz0)**4))
            end do
! FIN DU IF LTEXT
        endif
! FIN BOUCLE SUR LES PTS DE GAUSS
    end do
end subroutine
