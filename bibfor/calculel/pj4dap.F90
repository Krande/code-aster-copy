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
subroutine pj4dap(ino2, geom2, geom1, tria3, cobary, &
                  itr3, nbtrou, btdi, btvr, btnb, &
                  btlc, btco, l_dmax, dmax, dala, &
                  loin, dmin)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/pj3dgb.h"
#include "asterfort/pj4da2.h"
    real(kind=8) :: cobary(3), geom1(*), geom2(*), btvr(*)
    integer(kind=8) :: itr3, nbtrou, btdi(*), btnb(*), btlc(*), btco(*), tria3(*)
!  but :
!    trouver le tria3 qui servira a interpoler le noeud ino2
!    ainsi que les coordonnees barycentriques de ino2 dans ce tria3
!
!  in   ino2       i  : numero du noeud de m2 cherche
!  in   geom2(*)   r  : coordonnees des noeuds du maillage m2
!  in   geom1(*)   r  : coordonnees des noeuds du maillage m1
!  in   tria3(*)   i  : objet '&&pjxxco.tria3'
!  in   btdi(*)    i  : objet .bt3ddi de la sd boite_3d
!  in   btvr(*)    r  : objet .bt3dvr de la sd boite_3d
!  in   btnb(*)    i  : objet .bt3dnb de la sd boite_3d
!  in   btlc(*)    i  : objet .bt3dlc de la sd boite_3d
!  in   btco(*)    i  : objet .bt3dco de la sd boite_3d
!  in   l_dmax     l  : .true. : il faut prendre dmax en compte
!  in   dmax       r  : distance au dela de laquelle le noeud ino2
!                       ne sera pas projete.
!  in   dala       r  : distance au dela de laquelle le noeud ino2
!                       sera considere comme lointain
!  out  nbtrou     i  : 1 -> on a trouve 1 tria3 solution
!                     : 0 -> on n'a pas trouve de tria3 solution
!  out  itr3       i  : numero du tria3 solution
!  out  cobary(4)  r  : coordonnees barycentriques de ino2 dans itr3
!  out  dmin       r  : distance de ino2 au bord de itr3 si ino2 est
!                       exterieur a itr3.
!  out  loin       l  : .true. si dmin > 10% diametre(itr3) ou si dmin < dala
!
!  remarque :
!    si nbtrou=0, ino2 ne sera pas projete car il est au dela de dmax
!    alors : dmin=0, loin=.false.
! ----------------------------------------------------------------------
!
!
    real(kind=8) :: cobar2(3), dmin, d2, surf, rtr3
    integer(kind=8) :: p, q, r, p1, q1, p2, q2, r1, r2, ino2, i, k, iposi, nx, ny, ntrbt
!
    aster_logical :: l_dmax, loin
    real(kind=8) :: dmax, dala
! DEB ------------------------------------------------------------------
    nbtrou = 0
    loin = .false.
    dmin = 0.d0
!
    nx = btdi(1)
    ny = btdi(2)
!
!
!   --  on cherche le tria3 itr3 le plus proche de ino2 :
!   ------------------------------------------------------
    if (l_dmax) then
        dmin = dmax
    else
        dmin = r8maem()
    end if
!
!   -- on recherche la grosse boite candidate :
    call pj3dgb(ino2, geom2, geom1, tria3, 4, &
                btdi, btvr, btnb, btlc, btco, &
                p1, q1, r1, p2, q2, &
                r2)
    do p = p1, p2
        do q = q1, q2
            do r = r1, r2
                ntrbt = btnb((r-1)*nx*ny+(q-1)*nx+p)
                iposi = btlc((r-1)*nx*ny+(q-1)*nx+p)
                do k = 1, ntrbt
                    i = btco(iposi+k)
                    call pj4da2(ino2, geom2, i, geom1, tria3, &
                                cobar2, d2, surf)
                    if (sqrt(d2) .lt. dmin) then
                        rtr3 = surf
                        itr3 = i
                        dmin = sqrt(d2)
                        nbtrou = 1
                        cobary(1) = cobar2(1)
                        cobary(2) = cobar2(2)
                        cobary(3) = cobar2(3)
                    end if
                end do
            end do
        end do
    end do
!
!
!   -- calcul de loin :
    if (nbtrou .eq. 1) then
        if (dala .ge. 0.d0) then
            if (dmin .lt. dala) loin = .false.
        else
            if (rtr3 .eq. 0) then
                loin = .true.
            else
                rtr3 = rtr3**(1.d0/2.d0)
                if (dmin/rtr3 .gt. 1.d-1) loin = .true.
            end if
        end if
    else
        dmin = 0.d0
    end if
!
!
!
end subroutine
