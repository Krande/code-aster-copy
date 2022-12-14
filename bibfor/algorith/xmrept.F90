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
subroutine xmrept(jcesd, jcesv, jcesl, izone, ndim,&
                  ds_contact, geom, statue, mmait, amait,&
                  nmait)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8gaem.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/cfdisi.h"
#include "asterfort/cfmmvd.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/xxmmvd.h"
!
!
    type(NL_DS_Contact), intent(in) :: ds_contact
    real(kind=8) :: geom(3)
    integer :: ndim, mmait, amait, nmait, statue, izone
    integer :: jcesd(10), jcesv(10), jcesl(10)
!
! ----------------------------------------------------------------------
!
! ROUTINE XFEM (CONTACT - GRANDS GLISSEMENTS)
!
! ON CHERCHE LE POINT D'INTERSECTION MAITRE LE PLUS PROCHE DE
! POINT DE CONTACT
!
! TRAVAIL EFFECTUE EN COLLABORATION AVEC I.F.P.
!
! ----------------------------------------------------------------------
!
!  JCES*(1)  : POINTEURS DE LA SD SIMPLE NB DE FACETTES ET DE PT D'INTER
!  JCES*(2)  : POINTEURS DE LA SD SIMPLE DES INFOS SUR ARETES COUPEES
!  JCES*(6)  : POINTEURS DE LA SD SIMPLE DES COOR DES PT D'INTER MAITRE
! IN  IZONE  : NUM??RO DE ZONE DE LA MAILLE ESCLAVE
! In  ds_contact       : datastructure for contact management
! IN  NDIM   : DIMENSION DU MODELE
! IN  STATUE : STATUT DE LA MAILLE ESCLAVE
! IN  GEOM   : COORDONNEES DU POINT DE CONTACT ESCLAVE
! OUT MMAIT  : LE NUM??RO GLOBAL DE LA MAILLE MA??TRE CONTENANT LE POINT
!              D'INTERSECTION DE PLUS PROCHE
! OUT AMAIT  : LE NUM??RO LOCAL DE L'ARRETE CONTENANT LE POINT
!              D'INTERSECTION DE PLUS PROCHE
! OUT NMAIT  : LE NUM??RO LOCAL DU NOEUD CONTENANT LE POINT
!              D'INTERSECTION DE PLUS PROCHE
!
! ----------------------------------------------------------------------
!
    integer :: nummai, ntmae, ima, iad, nbpt
    integer :: zmesx
    integer :: ini, j, ifiss
    real(kind=8) :: coord(3), dmin, dist
    character(len=24) :: maescx
    integer :: jmaesx
    integer :: zxain
!
!   tolerances --- absolue et relative --- pour determiner si deux distances sont egales
    real(kind=8), parameter :: atol=1.e-12
    real(kind=8), parameter :: rtol=1.e-12
    aster_logical :: near
!
! --------------------------------------------------------------------
!
    call jemarq()
!
! --- SI LA MAILLE ESCLAVE EST CRACK-TIP, ON SORT
!
    if (statue .eq. 2 .or. statue .lt. 0) goto 999
!
! --- INITIALISATIONS
!
    dmin = r8gaem()
    coord(1:3)=0.d0
    ntmae = cfdisi(ds_contact%sdcont_defi,'NTMAE')
!
! --- RECUPERATION DE QUELQUES DONNEES
!
    maescx = ds_contact%sdcont_defi(1:16)//'.MAESCX'
    call jeveuo(maescx, 'L', jmaesx)
    zmesx = cfmmvd('ZMESX')
    zxain = xxmmvd('ZXAIN')
!
! --- BOUCLE SUR LES MAILLES FISSUR??ES
!
    do ima = 1, ntmae
!
! --- SI CE N'EST PAS LA BONNE ZONE, ON SORT
!
        if (zi(jmaesx+zmesx*(ima-1)+2-1) .ne. izone) goto 100
!
! --- FAUX PT D'INTEGRATION
!
        if (zi(jmaesx+zmesx*(ima-1)+4-1) .lt. 0) goto 100
!
        nummai = zi(jmaesx+zmesx*(ima-1)+1-1)
        ifiss = zi(jmaesx+zmesx*(ima-1)+5-1)
!
! ---  RECUPERATION DU NOMBRE DE POINTS D'INTERSECTION DE LA MAILLE
        call cesexi('C', jcesd(1), jcesl(1), nummai, 1,&
                    ifiss, 3, iad)
        ASSERT(iad.gt.0)
        nbpt = zi(jcesv(1)-1+iad)
! ----- BOUCLE SUR LES POINTS D'INTERSECTION
!
        do ini = 1, nbpt
! ------- COORDONNEES GEOMETRIQUES DU POINT D'INTERSECTION
!
            do j = 1, ndim
                call cesexi('S', jcesd(6), jcesl(6), nummai, 1,&
                            ifiss, ndim*( ini-1)+j, iad)
                ASSERT(iad.gt.0)
                coord(j)=zr(jcesv(6)-1+iad)
            end do
!
! ------- CALCUL DE LA DISTANCE
            dist = sqrt( ( coord(1)-geom(1))**2+ (coord(2)-geom(2))**2+ (coord(3)-geom(3) )**2 )
            call cesexi('S', jcesd(2), jcesl(2), nummai, 1,&
                        ifiss, zxain*( ini-1)+2, iad)
            if (nint(zr(jcesv(2)-1+iad)) .eq. 0 .and. ini .eq. 3 .and. ndim .eq. 2) cycle
!
!           dist est-elle egale a dmin ?
            near = abs(dist-dmin) .le. (atol + dmin*rtol)
!
            if (dist .lt. dmin .and. .not. near) then
                call cesexi('S', jcesd(2), jcesl(2), nummai, 1,&
                            ifiss, zxain*(ini-1)+1, iad)
                ASSERT(iad.gt.0)
                if (nint(zr(jcesv(2)-1+iad)) .gt. 0) then
                    amait = nint(zr(jcesv(2)-1+iad))
                    nmait = 0
                    dmin = dist
                    mmait = nummai
                else
                    call cesexi('S', jcesd(2), jcesl(2), nummai, 1,&
                                ifiss, zxain*(ini-1)+2, iad)
                    ASSERT(iad.gt.0)
                    if (nint(zr(jcesv(2)-1+iad)) .gt. 0) then
                        amait = 0
                        nmait = nint(zr(jcesv(2)-1+iad))
                        dmin = dist
                        mmait = nummai
                    endif
                endif
            endif
        end do
100     continue
    end do
!
999 continue
!
    call jedema()
end subroutine
