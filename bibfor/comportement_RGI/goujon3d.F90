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

subroutine goujon3d(endo, nbrenf, numr, numf, vecr, &
                    deq, rhor, wpl3, vwpl33, vwpl33t, &
                    rc, sigrf, sigrfissp)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================

    implicit none
#include "rgi_module.h"
#include "asterfort/matmat3d.h"

!Calcul du vecteur force pour un renfort traversant une fissure
!incluant l effet goujon : sortie=sigrfissp(numr,numf,direction)
!en base principale de fissuration

    integer(kind=8)      :: i, k

! Variable externes
    aster_logical :: endo
    integer(kind=8) :: nbrenf, numr, numf
    real(kind=8) :: vecr(nbrenf, 3), deq(nbrenf), rhor(nbrenf)
!     real(kind=8) :: wplt6(6)
    real(kind=8) :: wpl3(3), vwpl33(3, 3), vwpl33t(3, 3)
    real(kind=8), intent(in) :: rc
    real(kind=8) :: sigrf
    real(kind=8) :: sigrfissp(nbrenf, 3, 3)

! Variable locale
    real(kind=8) :: cosa, ws, sigc, theta, st, ct, rmin, gmax
    real(kind=8) :: g3(3), t3(3), d3(3), sigrf3(3), sigrf3p(3)

!Artefact de calcul
    real(kind=8) :: cmin, sigs, thetax, xg

!Initialisation du vecteur, nul par defaut
    do k = 1, 3
        sigrfissp(numr, numf, k) = 0.d0
    end do

!Initialisation des variables locales
!      ws=0.d0
!      sigc=0.d0
!      theta=0.d0
!      st=0.d0
!      ct=0.d0
!      rmin=0.d0
!      gmax=0.d0
!      g3(:)=0.d0
!      d3(:)=0.d0
!      sigrf3(:)=0.d0
!      sigrf3p(:)=0.d0

!Si la matrice est endommage,la contribution de l acier s initie
    if (endo) then
!Cas de la fissure numf et de l armature numr
        i = numf
!Cos de l angle entre la normale de numf et direction de numr
        cosa = 0.d0
        do k = 1, 3
            cosa = cosa+vecr(numr, k)*vwpl33(k, i)
        end do
!La valeur ne peut etre negative
        if (cosa .lt. 0.) then
!             print*,'Pb orientation R/D dans goujon3d fiss',i,'arm',numr
!             print*,'Reorientation du vecteur armature car cosinus negatif'
            do k = 1, 3
!              print*,vwpl33(k,i),vecr(numr,k)
                vecr(numr, k) = -vecr(numr, k)
            end do
            cosa = -cosa

        end if
!Si la valeur est quasi nulle
        if (cosa .ne. 0.) then
            if (cosa .le. 1.d0-EPSIL) then
!            l angle entre fissure et armature est < a Pi/2 ok
!            ouverture dans l axe de l armature (deplacement axial)
                ws = wpl3(i)*cosa
!            vecteur deplacement transversal au renfort (base fixe) et sa norme
                xg = 0.d0
                do k = 1, 3
                    g3(k) = wpl3(i)*vwpl33(k, i)-ws*vecr(numr, k)
                    xg = xg+g3(k)**2
                end do
!            norme du deplacement radial
                xg = sqrt(xg)
                if (cosa .lt. 0.) then
                    print *, 'goujon3d fiss', i, 'arm', numr
                    print *, 'xg', xg, 'wpl3(i)', wpl3(i)
                end if
!            Contrainte absolue dans l acier
                sigs = abs(sigrf)
!Si la contrainte de l acier est plus grande que zero
                if (sigs .gt. 0.) then
!                rayon de courbure minimal pour devier l armature
                    rmin = 0.785398d0*deq(numr)*sigs/rc
!                glissement maxi pour pouvoir calculer l'angle
                    gmax = 2.d0*rmin
!                vecteur unitaire de deplacement  transversal au renfort (base fixe)
                    if (xg .gt. 0.d0) then
                        do k = 1, 3
                            t3(k) = g3(k)/xg
                        end do
                    else
!                    pas de deplacement transversal
                        do k = 1, 3
                            t3(k) = 0.d0
                        end do
                    end if
!                contrainte normale dans la matrice dans la direction
!                transversale a l armature : rc
                    sigc = rc
!                les conditions d apparition d une bielle transverse
!               sont reunies, on calcule l inclinaison du renfort dans
!                la fissure
                    if (xg .le. gmax) then
!                  la deplacement transversal est suffisant
                        theta = 2.d0*asin(sqrt(xg/gmax))
                    else
                        theta = 2.d0*asin(sqrt(0.5d0))
                    end if
!                   verif de la condition de non depassement de |VI.Vr|
                    cmin = 0.d0
                    do k = 1, 3
                        cmin = cmin+vecr(numr, k)*vwpl33(k, i)
                    end do
                    cmin = min(abs(cmin), 1.d0)
                    thetax = acos(cmin)
                    theta = max(thetax, theta)
!                direction de la tension du renfort dans la fissure (base principale)
                    st = sin(theta)
                    ct = cos(theta)
                    do k = 1, 3
!                    base fixe
                        d3(k) = st*t3(k)+ct*vecr(numr, k)
!                    cosa prend en compte le nombre de renfort par unite de surface de fissure
!                    base fixe
                        sigrf3(k) = sigrf*rhor(numr)*cosa*d3(k)
!                     print*,'sigrf3',sigrf3(k),'sigrf',sigrf,'rho',rhor(numr),'d3',d3(k)

                    end do
!                passage du vecteur dans la base principale des endommagements
!                C(NL,NC2)=A(NL,NC1)*B(NC1,NC2)
!                MATMAT3d(A,B,NL,NC1,NC2,C)
                    call matmat3d(vwpl33t, sigrf3, 3, 3, 1, sigrf3p)
!       stockage en base principale des endommagements
!       vecteur contrainte armature suivant normale'

                    do k = 1, 3
                        sigrfissp(numr, numf, k) = sigrf3p(k)
                    end do
                else
!Si la contrainte de l acier est équivalente à zero
!                la contrainte dans l acier est nulle
!                print*,'********renfort', numr, 'fissure', i
!                 print*,'contrainte nulle dans l acier pour goujon3d'
                    do k = 1, 3
                        sigrfissp(numr, numf, k) = 0.d0
                    end do
                end if
            else
!                 print*,'********renfort', numr, 'fissure', i
!                 print*,'Calculsansevo vecteur=',vecr(numr,3),'numr',numr
                do k = 1, 3
                    sigrfissp(numr, numf, k) = sigrf*rhor(numr)*vecr(numr, k)
                end do
            end if
        else
!L armature est parallele a la fissure
!           print*,'********renfort', numr, 'fissure', i
!           print*,'armature et fissure // ds goujon3d'
            do k = 1, 3
                sigrfissp(numr, numf, k) = 0.d0
            end do
        end if
    else
!Si pas d endo
!        print*,'********renfort', numr, 'fissure', i
!        print*,'pas d endommagement ds goujon3d'
        do k = 1, 3
            sigrfissp(numr, numf, k) = 0.d0
        end do
    end if

end subroutine
