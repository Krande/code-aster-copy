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
subroutine projax(jvecpg, nbvec, nbordr, proaxe, iflag, &
                  rmima, jraxe)
! person_in_charge: van-xuan.tran at edf.fr
    implicit none
#include "jeveux.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jerazo.h"
#include "asterfort/proax0.h"
#include "asterfort/raxini.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbvec, nbordr, iflag(nbvec), jvecpg, jraxe
    real(kind=8) :: rmima(4*nbvec)
    character(len=16) :: proaxe
! ----------------------------------------------------------------------
! BUT: PROJETER SUR UN OU DEUX AXES LES POINTS REPRESANTANT LE
!      CISAILLEMENT TAU DANS LE PLAN u, v.
! ----------------------------------------------------------------------
! ARGUMENTS:
! VECPG     IN   R  : VECTEUR DE TRAVAIL CONTENANT LES
!                     COMPOSANTES u ET v DU VECTEUR TAU (CISAILLEMENT),
!                     POUR TOUS LES VECTEURS NORMAUX (n) ET TOUS LES
!                     NUMEROS D'ORDRE.
! NBVEC     IN   I  : NOMBRE DE VECTEURS NORMAUX.
! NBORDR    IN   I  : NOMBRE DE NUMERO D'ORDRE STOCKE DANS LA
!                     STRUCTURE DE DONNEES RESULTAT.
! PROAXE    IN   K16: TYPE DE PROJECTION (UN OU DEUX AXES).
! IFLAG     IN   I  : VECTEUR DE DRAPEAUX INDIQUANT :
!                      - IFLAG(i) = 0 --> CAS GENERAL
!                      - IFLAG(i) = 1 --> CAS OU LES POINTS DANS LE
!                                         PLAN DE CISAILLEMENT SONT
!                                         ALIGNES VERTICALEMENT.
!                      - IFLAG(i) = 2 --> CAS OU LES POINTS DANS LE
!                                         PLAN DE CISAILLEMENT SONT
!                                         ALIGNES HORIZONTALEMENT.
!                      - IFLAG(i) = 3 --> CAS OU LES POINTS DANS LE
!                                         PLAN DE CISAILLEMENT SONT
!                                         CONTENUS DANS UN CADRE DE
!                                         COTES INFERIEURS A EPSILO.
! RMIMA     IN   R  : VECTEUR CONTENANT LES COORDONNEES DES POINTS
!                     EXTREMES DU CADRE (UMIN, UMAX, VMIN, VMAX).
!                     POUR TOUS LES VECTEURS NORMAUX.
! RAXE      OUT  R  : VECTEUR CONTENANT L'AMPLITUDE DES POINTS
!                     PROJETES.
!
!-----------------------------------------------------------------------
!     ------------------------------------------------------------------
    integer(kind=8) :: i, n1, ivect, iordr, jsec1, jsec2, jsec3, jsec4
    integer(kind=8) :: indsec, nbpts1, nbpts2, nbpts3, nbpts4, nptsec(4)
!
    real(kind=8) :: umin, umax, vmin, vmax, u0, v0, diamin
    real(kind=8) :: cstai, cstbi, cstas, cstbs, uaxe2, vaxe2
    real(kind=8) :: ui, vi, dists(4), ai, bi, as, bs, rpaxi, rpaxs
    real(kind=8) :: up, vp, val
!
    character(len=4) :: axeini
!
!-----------------------------------------------------------------------
!234567                                                              012
!
    call jemarq()
!
    call wkvect('&&PROJAX.SECT1', 'V V R', nbordr*2, jsec1)
    call wkvect('&&PROJAX.SECT2', 'V V R', nbordr*2, jsec2)
    call wkvect('&&PROJAX.SECT3', 'V V R', nbordr*2, jsec3)
    call wkvect('&&PROJAX.SECT4', 'V V R', nbordr*2, jsec4)
!
    n1 = 0
!
    do ivect = 1, nbvec
        call jerazo('&&PROJAX.SECT1', nbordr*2, 1)
        call jerazo('&&PROJAX.SECT2', nbordr*2, 1)
        call jerazo('&&PROJAX.SECT3', nbordr*2, 1)
        call jerazo('&&PROJAX.SECT4', nbordr*2, 1)
!
        umin = rmima(4*(ivect-1)+1)
        umax = rmima(4*(ivect-1)+2)
        vmin = rmima(4*(ivect-1)+3)
        vmax = rmima(4*(ivect-1)+4)
        u0 = umin+(umax-umin)/2.0d0
        v0 = vmin+(vmax-vmin)/2.0d0
!
        if ((umax-umin) .ge. (vmax-vmin)) then
            diamin = umax-umin
        else
            diamin = vmax-vmin
        end if
!
! 1. DETERMINATION DES POINTS SITUES DANS LES 4 SECTEURS
!
        nbpts1 = 0
        nbpts2 = 0
        nbpts3 = 0
        nbpts4 = 0
!
        nptsec(1) = 0
        nptsec(2) = 0
        nptsec(3) = 0
        nptsec(4) = 0
!
        do iordr = 1, nbordr
            n1 = n1+1
            ui = zr(jvecpg+2*n1-1)
            vi = zr(jvecpg+2*n1)
!
            dists(1) = sqrt((umin-ui)**2+(vmax-vi)**2)
            dists(2) = sqrt((umax-ui)**2+(vmax-vi)**2)
            dists(3) = sqrt((umax-ui)**2+(vmin-vi)**2)
            dists(4) = sqrt((umin-ui)**2+(vmin-vi)**2)
!
            indsec = 0
            do i = 1, 4
                if ((dists(i) .ge. diamin) .and. (indsec .eq. 0)) then
                    if (ui .ge. u0) then
                        if (vi .ge. v0) then
                            zr(jsec2+nbpts2*2) = ui
                            zr(jsec2+nbpts2*2+1) = vi
                            nbpts2 = nbpts2+1
                            indsec = 1
                        else
                            zr(jsec3+nbpts3*2) = ui
                            zr(jsec3+nbpts3*2+1) = vi
                            nbpts3 = nbpts3+1
                            indsec = 1
                        end if
                    else
                        if (vi .ge. v0) then
                            zr(jsec1+nbpts1*2) = ui
                            zr(jsec1+nbpts1*2+1) = vi
                            nbpts1 = nbpts1+1
                            indsec = 1
                        else
                            zr(jsec4+nbpts4*2) = ui
                            zr(jsec4+nbpts4*2+1) = vi
                            nbpts4 = nbpts4+1
                            indsec = 1
                        end if
                    end if
                end if
            end do
        end do
!
! 2. CHOIX DE L'AXE INITIAL (SI UN 2EME AXE EST DEMANDE,
!                            IL SERA DEDUIT DE L'AXE INITIAL.)
!
!      A ----------------------------- B
!        | Sect.1             Sect.2 |    PAR DEFINITION L'AXE 1 EST LE
!        |                           |    SEGMENT DB.
!        |                           |
!        |                           |    PAR DEFINITION L'AXE 2 EST LE
!        |                           |    SEGMENT AC.
!        | Sect.4             Sect.3 |
!      D ----------------------------- C
!
!
        if (iflag(ivect) .eq. 0) then
            n1 = n1-nbordr
!
            if (( &
                (nbpts1 .eq. 0) .and. (nbpts3 .eq. 0) .and. (nbpts2 .gt. 0) .and. &
                (nbpts4 .gt. 0) &
                ) &
                .or. &
                ( &
                (nbpts1 .eq. 0) .and. (nbpts3 .gt. 0) .and. (nbpts2 .gt. 0) .and. &
                (nbpts4 .gt. 0) &
                ) &
                .or. &
                ( &
                (nbpts1 .gt. 0) .and. (nbpts3 .eq. 0) .and. (nbpts2 .gt. 0) .and. &
                (nbpts4 .gt. 0) &
                )) then
                axeini = 'AXE1'
            elseif (((nbpts1 .gt. 0) .and. (nbpts3 .gt. 0) .and. &
                     (nbpts2 .eq. 0) .and. (nbpts4 .eq. 0)) .or. ((nbpts1 .gt. &
                                            0) .and. (nbpts3 .gt. 0) .and. (nbpts2 .eq. 0) .and. ( &
                                          nbpts4 .gt. 0)) .or. ((nbpts1 .gt. 0) .and. (nbpts3 .gt. &
                                               0) .and. (nbpts2 .gt. 0) .and. (nbpts4 .eq. 0))) then
                axeini = 'AXE2'
            elseif (((nbpts1 .eq. 0) .and. (nbpts3 .gt. 0) .and. &
                     (nbpts2 .eq. 0) .and. (nbpts4 .gt. 0)) .or. ((nbpts1 .eq. &
                                            0) .and. (nbpts3 .gt. 0) .and. (nbpts2 .gt. 0) .and. ( &
                                          nbpts4 .eq. 0)) .or. ((nbpts1 .gt. 0) .and. (nbpts3 .eq. &
                                              0) .and. (nbpts2 .eq. 0) .and. (nbpts4 .gt. 0)) .or. &
                    ((nbpts1 .gt. 0) .and. (nbpts3 .eq. 0) .and. (nbpts2 .gt. &
                                           0) .and. (nbpts4 .eq. 0)) .or. ((nbpts1 .gt. 0) .and. ( &
                                        nbpts3 .gt. 0) .and. (nbpts2 .gt. 0) .and. (nbpts4 .gt. 0) &
                                                                                             )) then
                nptsec(1) = nbpts1
                nptsec(2) = nbpts2
                nptsec(3) = nbpts3
                nptsec(4) = nbpts4
                call raxini(zr(jsec1), zr(jsec2), zr(jsec3), zr(jsec4), nptsec, &
                            nbordr, umin, umax, vmin, vmax, &
                            axeini)
            else
                call utmess('F', 'PREPOST4_58')
            end if
!
! 3. CALCUL DES CONSTANTES NECESSAIRES A LA PROJECTION SUR UN OU DEUX
!    AXES
!
            if (axeini .eq. 'AXE1') then
                cstai = umax-umin
                cstbi = vmax-vmin
                ai = (vmax-vmin)/(umax-umin)
                bi = (umax*vmin-umin*vmax)/(umax-umin)
            else if (axeini .eq. 'AXE2') then
                cstai = -(umax-umin)
                cstbi = vmax-vmin
                ai = -(vmax-vmin)/(umax-umin)
                bi = (umax*vmax-umin*vmin)/(umax-umin)
            end if
!
            if (proaxe .eq. 'DEUX_AXES') then
                uaxe2 = u0+1.0d0
                vaxe2 = v0-(cstai/cstbi)*(uaxe2-u0)
                cstas = uaxe2-u0
                cstbs = vaxe2-v0
                as = (vaxe2-v0)/(uaxe2-u0)
                bs = (uaxe2*v0-u0*vaxe2)/(uaxe2-u0)
            end if
!
! 4. PROJECTION SUR UN AXE OU DEUX AXES
!
            do iordr = 1, nbordr
                n1 = n1+1
!
                ui = zr(jvecpg+2*n1-1)
                vi = zr(jvecpg+2*n1)
!
! 4.1 PROJECTION SUR L'AXE INITIAL
!
                call proax0(ui, vi, cstai, cstbi, ai, &
                            bi, u0, v0, rpaxi)
!
! 4.2 PROJECTION SUR LE SECOND AXE SI CELA EST DEMANDE
!
                if (proaxe .eq. 'DEUX_AXES') then
                    call proax0(ui, vi, cstas, cstbs, as, &
                                bs, u0, v0, rpaxs)
                else
                    rpaxs = 0.0d0
                end if
!
! 4.3 CALCUL DU MODULE ET ATTRIBUTION DU SIGNE
!
                if (rpaxi .lt. 0.0d0) then
                    zr(jraxe+n1) = -sqrt(rpaxi**2+rpaxs**2)
                else
                    zr(jraxe+n1) = sqrt(rpaxi**2+rpaxs**2)
                end if
            end do
!
! LES POINTS SONT ALIGNES VERTICALEMENT (PAS DE PROJECTION)
        else if (iflag(ivect) .eq. 1) then
            n1 = n1-nbordr
!
            do iordr = 1, nbordr
                n1 = n1+1
!
                ui = zr(jvecpg+2*n1-1)
                vi = zr(jvecpg+2*n1)
                up = ui
                vp = vi
                val = sqrt((up-u0)**2+(vp-v0)**2)
                if (vp .lt. v0) then
                    val = -val
                end if
                zr(jraxe+n1) = val
            end do
!
! LES POINTS SONT ALIGNES HORIZONTALEMENT (PAS DE PROJECTION)
        else if (iflag(ivect) .eq. 2) then
            n1 = n1-nbordr
!
            do iordr = 1, nbordr
                n1 = n1+1
!
                ui = zr(jvecpg+2*n1-1)
                vi = zr(jvecpg+2*n1)
                up = ui
                vp = vi
                val = sqrt((up-u0)**2+(vp-v0)**2)
                if (up .lt. u0) then
                    val = -val
                end if
                zr(jraxe+n1) = val
            end do
!
! LES POINTS SONT DANS UN CADRE DONT LES COTES SONT INFERIEURS A EPSILO
! (PAS DE PROJECTION)
        else if (iflag(ivect) .eq. 3) then
            n1 = n1-nbordr
!
            do iordr = 1, nbordr
                n1 = n1+1
!
                ui = zr(jvecpg+2*n1-1)
                vi = zr(jvecpg+2*n1)
                up = ui
                vp = vi
                val = sqrt((up-u0)**2+(vp-v0)**2)
                if (up .lt. u0) then
                    val = -val
                end if
                zr(jraxe+n1) = val
            end do
!
        end if
!
    end do
!
    call jedetr('&&PROJAX.SECT1')
    call jedetr('&&PROJAX.SECT2')
    call jedetr('&&PROJAX.SECT3')
    call jedetr('&&PROJAX.SECT4')
!
    call jedema()
end subroutine
