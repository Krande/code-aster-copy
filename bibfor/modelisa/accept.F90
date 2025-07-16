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

subroutine accept(f, nbm, method, imode, jmode, &
                  jc, dir, uc, uct, &
                  l, lt, val_spec)
    implicit none
!     OPERATEUR PROJ_SPEC_BASE
!     PROJECTION D UN OU PLUSIEURS SPECTRES DE TURBULENCE SUR UNE BASE
!     MODALE PERTURBEE PAR PRISE EN COMPTE DU COUPLAGE FLUIDE STRUCTURE
!-----------------------------------------------------------------------
!  IN   : IEX = 0 => ON DEMANDE L EXECUTION DE LA COMMANDE
!         IEX = 1 => ON NE FAIT QUE VERIFIER LES PRECONDITIONS
!  OUT  : IER = 0 => TOUT S EST BIEN PASSE
!         IER > 0 => NOMBRE D ERREURS RENCONTREES
!-----------------------------------------------------------------------
!
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/coegen.h"
#include "asterfort/corcos.h"
#include "asterfort/coyang.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
    integer(kind=8) :: nbm, i, ind, jnd, ipg, jpg, iad1, iad2, itab, imode, jmode
    integer(kind=8) :: ispe, jspe, iorig, jorig, ntail, ntail1, ntail2
    integer(kind=8) :: jpgini, jpgfin
    real(kind=8) :: f, pi, deuxpi, omega, uc, kl, kt, dir(3, 3)
    real(kind=8) :: d1, d2, d3, mes(3), coeh, jc, l, local(3)
    real(kind=8) :: uct, lt, rayon, dist, dteta, jc1
    character(len=8) :: method
    character(len=24) :: val_spec
!
    data(local(i), i=1, 3)/3*0.d0/
!
!
!-----------------------------------------------------------------------
    call jemarq()
!
! QUELQUES CONSTANTES
    pi = r8pi()
    deuxpi = 2*pi
    omega = deuxpi*f
!
! Condition fictive pour que les bonnes longueurs de corrélation soient
! utilisees si spec_corr_conv1 est choisi avec corcos
    if (method(1:7) .ne. 'AU_YANG') then
        if (val_spec .eq. 'SPEC_CORR_CONV_1') then
            kl = 0.1d0*omega/uc
            kt = 0.55d0*omega/uc
            l = 1/kl
            lt = 1/kt
        end if
    end if
!
! BOUCLE SUR LES ELEMENTS DU MODELE
    jc = 0.d0
    jc1 = 0.d0

    call jeveuo('&&GROTAB.TAB', 'L', itab)
    call jelira('&&GROTAB.TAB', 'LONUTI', ntail)

    if (method(1:7) .eq. 'AU_YANG') rayon = zr(itab+ntail-1)
    ntail = ntail-1
    ntail1 = ntail/(6*nbm)
    ntail2 = ntail/nbm
    iorig = (imode-1)*ntail2
    jorig = (jmode-1)*ntail2
    jpgfin = (jmode-1)*ntail2+ntail1-1
    do ipg = (imode-1)*ntail2, (imode-1)*ntail2+ntail1-1
        if (jmode .eq. imode) then
            jpgini = ipg
            !jpgini = (imode-1)*ntail2
        else
            jpgini = (jmode-1)*ntail
        end if
        ispe = ipg-(imode-1)*ntail2
        iad1 = itab+iorig+6*ispe
        do jpg = jpgini, jpgfin
            jspe = jpg-(jmode-1)*ntail2
            iad2 = itab+jorig+6*jspe
! CALCUL DISTANCES INTER POINTS DE GAUSS
! ABSCISSES
            mes(1) = zr(iad1+1)-zr(iad2+1)
! ORDONNEES
            mes(2) = zr(iad1+2)-zr(iad2+2)
!
! COTE
            mes(3) = zr(iad1+3)-zr(iad2+3)
!
! COHERENCE CORCOS
            if (method(1:6) .eq. 'CORCOS') then
                do ind = 1, 3
                    local(ind) = 0.d0
                    do jnd = 1, 3
                        local(ind) = local(ind)+dir(ind, jnd)*mes(jnd)
                    end do
                end do
!
                d1 = abs(local(1))
                d2 = abs(local(2))
                d3 = abs(local(3))
!
                coeh = corcos(d1, d2, local(1), local(2), uc, uct, l, lt, omega)
!
! COHERENCE GENERALE
            else if (method .eq. 'GENERALE') then
                d1 = abs(mes(1))
                d2 = abs(mes(2))
                d3 = abs(mes(3))
                coeh = coegen(d1, d2, d3, l, omega, uc)
            else if (method(1:7) .eq. 'AU_YANG') then
                dist = zr(iad1+4)-zr(iad2+4)
                dteta = zr(iad1+5)-zr(iad2+5)
! on enleve les abs devant les distances
                coeh = coyang(dist, dteta, rayon, omega, uc, uct, l, lt)
            end if
            ! if (jmode .eq. imode .and. jpg .gt. ipg) then
            if (jmode .eq. imode .and. jpg .gt. ipg) then
                jc1 = jc1+coeh*zr(iad1)*zr(iad2)
            else
                jc = jc+coeh*zr(iad1)*zr(iad2)

            end if
        end do
    end do
    if (imode .eq. jmode) jc = jc+2*jc1
    !
    call jedema()
end subroutine
