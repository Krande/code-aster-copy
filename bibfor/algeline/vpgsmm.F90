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
subroutine vpgsmm(nbeq, nfreq, nconv, vect, alpha, &
                  lmatb, typeps, vaux, ddlexc, delta, &
                  dsor, omecor)
!---------------------------------------------------------------------
!
!
!     SUBROUTINE CHAPEAU DE VPGSKP ORTHONORMALISANT UNIQUEMENT AU SEIN
!     DES PAQUETS DE MODES MULTIPLES.
!     DOIT ETRE APPELES UNIQUEMENT APRES UN VPORDO (EN MODULE OU EN
!     VALEUR ABSOLUE).
!     CETTE ROUTINE N'EST VRAIMENT UTILE EN GAIN DE TEMPS QUE SI LE
!     PD SCALAIRE EST MATRICIEL, NBEQ GRANDS ET SURTOUT NCONV > 10.
!
!     RQ1: MEME PARAMETRES D'APPEL QUE VPGSKP +
!          DSOR: COMPOSANTES DES VECTEURS PROPRES (CF.VPSORN)
!          OMECOR: PARAMETRE POUR CRITERES MODES MULTIPLES (CF.VPSORN)
!     RQ2: CRITERES MODES MULTIPLES BASES SUR SEUILR/P.
!
!     NIVEAU DEVELOPPEUR: EN MODE SURCHARGE
!     RQ3: TESTER LA M-ORTHOGONALITE VIA LE LOGICAL LCHECK
!     RQ4: EN POSANT ALPHA=-ALPHA, M-ORTHOGONALISATION COMPLETE
!-----------------------------------------------------------------------
! CORPS DU PROGRAMME
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mrmult.h"
#include "asterfort/vpgskp.h"
#include "asterfort/wkvect.h"
!
!-----------------------------------------------------------------------
! DECLARATION PARAMETRES D'APPELS
    integer(kind=8) :: nbeq, nfreq, nconv, lmatb, ddlexc(nbeq), typeps
    real(kind=8) :: vect(nbeq, nconv), alpha, vaux(nbeq), delta(nconv), dsor(nfreq+1, 2), omecor
!
!-----------------------------------------------------------------------
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: ireor1, ireor2, i, j, k, compt1, compt2, compt3
    integer(kind=8) :: iauxi, iauxj, nconvl, ifm, niv
    real(kind=8) :: seuilr, seuilp, auxri, auxri1
    aster_logical :: ltest, lcheck
!
!-----------------------------------------------------------------------
    call jemarq()
    call infniv(ifm, niv)
!
! flag pour tester la m-orthogonalite avant/apres (niveau developpeur)
!      lcheck=.true.
    lcheck = .false.
! flag pour lancer la m-orthogonalite totale (niveau developpeur)
!      alpha=-alpha
    if (lcheck) then
        write (ifm, *) 'm-orthogonalite std avant vpgsmm'
        do i = 1, nconv
            call mrmult('ZERO', lmatb, vect(1, i), vaux, 1, &
                        .false._1)
            auxri = 0.d0
            do j = 1, i
                auxri1 = 0.d0
                do k = 1, nbeq
                    auxri1 = auxri1+vect(k, j)*vaux(k)*ddlexc(k)
                end do
                auxri1 = abs(auxri1)
                if (j .ne. i) auxri = auxri+auxri1
            end do
            write (ifm, *) 'sigma(j=1,i-1) abs<vj,mvi>=', i, auxri
            write (ifm, *) '               abs<vi,mvi>=     ', auxri1
        end do
    end if
!
    if (alpha .lt. 0.d0) then
        alpha = -alpha
        write (ifm, *) '<vpgsmm> m-reorthogonalisation totale'
        call vpgskp(nbeq, nconv, vect, alpha, lmatb, &
                    typeps, vaux, ddlexc, delta)
        goto 72
    end if
!
    call wkvect('&&VPGSMM.REORTHO.PART1', 'V V I', nconv, ireor1)
    call wkvect('&&VPGSMM.REORTHO.PART2', 'V V I', 2*nconv, ireor2)
    do i = 1, nconv
        zi(ireor1+i-1) = -999
        zi(ireor2+i-1) = -999
        zi(ireor2+nconv+i-1) = -999
!         write(6,*)i,dsor(i,1)
    end do
!
! Critères de test dependant de CALC_* / SEUIL_* (*=FREQ ou CHAR_CRIT):
! * en dessous de seuilr: modes rigides
! * si deux modes ne sont pas considérés comme des modes rigides :
!     modes multiples si leur écart relatif < seuilp
! rq: seuils volontairement lâches car il vaut mieux trop réorthogonaliser que l'inverse
    seuilr = 100.d0*omecor
    seuilp = omecor
    compt1 = 0
    do i = 1, nconv-1
        auxri = dsor(i, 1)
        auxri1 = dsor(i+1, 1)
        ltest = .false.
        if (abs(auxri) .lt. seuilr) then
            ltest = (abs(auxri1) .lt. seuilr)
        else
            if (abs(auxri1) .ge. seuilr) then
                ltest = 2.d0*abs(auxri-auxri1) .lt. abs(auxri+auxri1)*seuilp
            else
                ltest = .false.
            end if
        end if
        if (ltest) then
            if (zi(ireor1+i-1) .eq. -999) then
                compt1 = compt1+1
                zi(ireor1+i-1) = compt1
                zi(ireor1+i) = compt1
            else
                zi(ireor1+i) = zi(ireor1+i-1)
            end if
        end if
    end do
    i = 1
    compt1 = 0
    compt2 = 1
68  continue
    iauxi = zi(ireor1+i-1)
    if (iauxi .ne. -999) then
        compt1 = compt1+1
        do j = i+1, nconv
            iauxj = zi(ireor1+j-1)
            if (iauxj .ne. iauxi) then
                compt2 = 1
                goto 70
            else
                compt2 = compt2+1
                zi(ireor2+compt1-1) = i
                zi(ireor2+nconv+compt1-1) = compt2
            end if
        end do
70      continue
        i = j
        if (i .lt. nconv) goto 68
    else
        i = i+1
        if (i .lt. nconv) goto 68
    end if
!
    compt3 = nconv
    iauxj = 1
    if ((lcheck) .or. (niv .ge. 2)) then
        write (ifm, *) '<vpgsmm> m-reorthogonalisation partielle sur ', compt1, ' paquets'
    end if
    do i = 1, compt1
        iauxi = zi(ireor2+i-1)
        nconvl = zi(ireor2+nconv+i-1)
! gardes-fous
        if ((iauxi .lt. iauxj) .or. (iauxi .gt. (nconv-1))) then
            ASSERT(.false.)
        end if
        if ((nconvl .lt. 2) .or. (nconvl .gt. compt3)) then
            ASSERT(.false.)
        end if
!
        if ((lcheck) .or. (niv .ge. 2)) then
            write (ifm, *) '<vpgsmm> paquet modes multiples n', i, iauxi, nconvl
        end if
        call vpgskp(nbeq, nconvl, vect(1, iauxi), alpha, lmatb, &
                    typeps, vaux, ddlexc, delta)
!
! mise a jour des gardes-fous
        compt3 = compt3-nconvl
        if (compt3 .lt. 0) then
            ASSERT(.false.)
        end if
        iauxj = iauxi
    end do
    call jedetr('&&VPGSMM.REORTHO.PART1')
    call jedetr('&&VPGSMM.REORTHO.PART2')
!
72  continue
    if (lcheck) then
        write (ifm, *) 'm-orthogonalite std apres vpgsmm'
        do i = 1, nconv
            call mrmult('ZERO', lmatb, vect(1, i), vaux, 1, &
                        .false._1)
            auxri = 0.d0
            do j = 1, i
                auxri1 = 0.d0
                do k = 1, nbeq
                    auxri1 = auxri1+vect(k, j)*vaux(k)*ddlexc(k)
                end do
                auxri1 = abs(auxri1)
                if (j .ne. i) auxri = auxri+auxri1
            end do
            write (ifm, *) 'sigma(j=1,i-1) abs<vj,mvi>=', i, auxri
            write (ifm, *) '               abs<vi,mvi>=     ', auxri1
        end do
    end if
!
!     FIN DE VPGSMM
!
    call jedema()
end subroutine
