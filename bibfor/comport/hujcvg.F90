! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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

subroutine hujcvg(nmat, mater, nvi, vind, vinf, &
                  vins, nr, yd, dy, r, &
                  indi, iter, itmax, intg, toler, &
                  bnews, mtrac, ye, lreli, iret)
! person_in_charge: alexandre.foucault at edf.fr
! aslint: disable=W1306
    implicit none
!     ----------------------------------------------------------------
!     CRITERE DE CONVERGENCE SUR SOLUTION PROPOSEE PAR NEWTON
!     ----------------------------------------------------------------
!     IN:  MATER  : DONNEES MATERIAU
!          NMAT   : DIMENSION TABLEAU DONNEES MATERIAU
!          NVI    : NOMBRE DE VARIABLES INTERNES
!          VIND   : VARIABLES INTERNES A T
!          VINF   : VARIABLES INTERNES A T+DT
!          VINS   : VARIABLES INTERNES A T (NON MODIFIEES)
!          NR     : DIMENSION DU SYSTEME NL
!          YD     : VECTEUR SOLUTION A T
!          DY     : INCREMENT DES GRANDEURS INCONNUES A T+DT
!          R      : RESIDU DU SYSTEME NL
!          INDI   : INDICE DES MECANISMES POTENTIEL ACTIFS
!          ITER   : NUMERO ITERATION NEWTON LOCAL
!          ITMAX  : NB MAXI TERATIONS LOCALES
!          INTG   : NUMERO D'INTEGRATION COURANTE
!          TOLER  : TOLERANCE A CONVERGENCE
!          BNEWS  : INDICATEUR LIE A LA TRACTION
!          MTRAC  : INDEICATEUR LIE A LA TRACTION (BIS)
!          YE     : VALEURS DES INCONNUES APRES LCINIT
!          LRELI  : TYPE DE SCHEMA D'INTEGRATION
!     OUT: IRET   : = 0 CONVERGENCE
!                   = 1 ITERATION SUIVANTE
!                   = 2 RE-INTEGRATION
!                   = 3 REDECOUPAGE DU PAS DE TEMPS
!          BNEWS  : INDICATEUR LIE A LA TRACTION MAJ
!     ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterc/r8prem.h"
#include "asterfort/hujprj.h"
#include "asterfort/hujpxd.h"

    integer :: nvi, nr, nmat
    real(kind=8) :: mater(nmat, 2), vind(nvi), vinf(nvi), vins(nvi)
    real(kind=8) :: yd(nr), dy(nr), r(nr), toler, ye(nr)
    integer :: indi(7), iter, iret, itmax, intg
    aster_logical :: bnews(3), mtrac, lreli
!
    integer :: nbmeca, nbmect, ndt, k, i, j, msup(4), resi, imin, ndi
    integer :: kk, jj
    real(kind=8) :: err, yf(nr), ydt(nr), yft(nr), dev(3), pf, qf
    real(kind=8) :: e0, pref, rtrac, ptrac, zero, maxi, ratio
    real(kind=8) :: un, deux, yet(nr), lamin, cinq, prob(3), matert(22, 2)
    aster_logical :: tracti, noconv, negtra, prox(4), proxc(4), modif
    aster_logical :: neglam(3), cycl, euler, ltry, impose, probt
!
    parameter(ndi=3)
    parameter(ndt=6)
    parameter(zero=0.d0)
    parameter(un=1.d0)
    parameter(deux=2.d0)
    parameter(cinq=5.d0)
!     ----------------------------------------------------------------
! --- PARAMETRES MATERIAU
    e0 = mater(1, 1)
    pref = mater(8, 2)
    ptrac = mater(21, 2)
    rtrac = abs(pref*1.d-6)
!
    tracti = .false.
    neglam(1) = .false.
    neglam(2) = .false.
    neglam(3) = .false.
!
! --- DETERMINATION DU NOMBRE DE MECANISMES POTENTIELLEMENT ACTIFS
    nbmeca = 0
    nbmect = 0
    do k = 1, 7
        if (indi(k) .gt. 0) then
            nbmect = nbmect+1
            if (indi(k) .le. 8) nbmeca = nbmeca+1
        end if
    end do
!
! --- MISE A L'ECHELLE DES CONTRAINTES ET VARIABLES CONTENUES DANS YF
! --- SIGMA * E0, R * E0/PREF
    yf(1:nr) = yd(1:nr)+dy(1:nr)
!
! --- COPIE A PARTIR DU TRAITEMENT DE HUJMID
    ydt(1:nr) = yd(1:nr)
    yft(1:nr) = yf(1:nr)
    yet(1:nr) = ye(1:nr)
!
    do i = 1, ndt
        ydt(i) = yd(i)*e0
        yft(i) = yf(i)*e0
        yet(i) = ye(i)*e0
    end do
    do i = 1, nbmeca
        ydt(ndt+1+i) = yd(ndt+1+i)*e0/abs(pref)
        yft(ndt+1+i) = yf(ndt+1+i)*e0/abs(pref)
        yet(ndt+1+i) = ye(ndt+1+i)*e0/abs(pref)
    end do
! --- CONTROLE DU NIVEAU DE CONVERGENCE ATTEINT
    err = norm2(r)
!
! --- INCREMENT DU NOMBRE D'INTEGRATION COURANTE TENTEE
    if (iter .eq. 1) intg = intg+1
!
! --- INCREMENT DU NOMBRE D'ITERATIONS LOCALES
    vinf(35) = vinf(35)+1
!
! --- SI ITERATION MAXIMALE ATTEINTE
! --- POSSIBILITE DE RELANCER SUIVANT ALGO DE HUJMID
!
    if (iter .le. itmax) then
! ----------------------------------
! --- CONTROLE DU NOMBRE DE BOUCLES
! ----------------------------------
        if (intg .gt. 5) then
            iret = 3
            goto 999
        end if
! -------------------------
! ----   CONVERVENCE   ----
! -------------------------
        if (err .lt. toler) then
            goto 60
        else
! --------------------------------------------
! --- NON CONVERGENCE : ITERATION SUIVANTE ---
! --------------------------------------------
            iret = 1
            if ((nbmeca .ne. nbmect) .and. (nbmeca .eq. 0)) then
                if (err .gt. 1d5) then
                    iret = 3
                    goto 999
                end if
            else
                do i = 1, 3
                    call hujprj(i, yft, dev, pf, qf)
                    if (((pf+rtrac-ptrac)/abs(pref)) .ge. -r8prem()) then
                        do j = 1, nbmeca
                            if ((indi(j) .eq. i) .or. (indi(j) .eq. (i+4))) then
                                tracti = .true.
                                goto 999
                            end if
                        end do
                    end if
                    if (abs(pf) .gt. e0*1.d1) then
                        iret = 3
                        goto 999
                    end if
                end do
            end if
            goto 1000
        end if
    else
! --- SI NEWTON_RELI ACTIF - ON LE RETIRE
        if (lreli) then
            iret = 2
            lreli = .false.
            goto 1000
        else
! --- ITERATIONS MAXI ATTEINT
            iret = 3
            goto 999
        end if
    end if
!
60  continue
! ---------------------------------------
! --- CONTROLE DE LA SOLUTION OBTENUE ---
! ---------------------------------------
! -------------------------------------------------
! ---- VERIFICATION DES MULTIPLICATEURS PLASTIQUES
! -------------------------------------------------
    maxi = toler
    do k = 1, nbmect
        if (yft(ndt+1+nbmeca+k) .gt. maxi) maxi = yft(ndt+1+nbmeca+k)
    end do
!
    negtra = .false.
!
    do k = 1, nbmect
        ratio = yft(ndt+1+nbmeca+k)/maxi
        if (ratio .lt. (-toler)) then
            if (indi(k) .gt. 8) then
! ----------------------------------------------
! ---> MECANISME DE TRACTION
! LAMBDA < 0 --> DESACTIVATION DU MECANISME POUR
! LA PROCHAINE TENTATIVE D'INTEGRATION
! ----------------------------------------------
                bnews(indi(k)-8) = .true.
                negtra = .true.
            end if
        end if
    end do
!
! -------------------------------------------------------
! ---> MECANISME DE TRACTION
! LAMBDA < 0 --> TENTATIVE SUPPL D'INTEGRATION SI COMPT<5
! -------------------------------------------------------
    if (negtra) then
        if (intg .gt. 5) then
            iret = 3
            goto 1000
        else if (nbmeca .ne. 0) then
            vinf(1:nvi) = vind(1:nvi)
        end if
        iret = 2
    end if
!
! -------------------------------------------------------
! --- ON CONTROLE QUE LES PRESSIONS ISOTROPES PAR PLAN
!     NE PRESENTENT PAS DE TRACTION
! -------------------------------------------------------
    noconv = .false.
    do i = 1, ndi
        call hujprj(i, yft, dev, pf, qf)
        if (((pf+rtrac-ptrac)/abs(pref)) .gt. zero) then
            noconv = .true.
            bnews(i) = .false.
            iret = 2
        end if
    end do
!
! -------------------------------------------------------
! --- SI TRACTION DETECTEE ET NON CONVERGENCE, ON IMPOSE
! --- ETAT DE CONTRAINTES ISOTROPE
! -------------------------------------------------------
    if ((noconv) .and. (nbmect .ne. nbmeca)) then
        noconv = .false.
        do i = 1, 3
            vind(23+i) = zero
            vind(27+i) = zero
            vind(5+4*i) = zero
            vind(6+4*i) = zero
            vind(7+4*i) = zero
            vind(8+4*i) = zero
        end do
        do i = 1, ndi
            dy(i) = -yd(i)-deux*rtrac/e0
            dy(ndi+i) = -yd(i)
        end do
        do i = 1, nbmeca+nbmect+1
            dy(ndt+i) = zero
        end do
        vinf(1:nvi) = vind(1:nvi)
        iret = 0
    end if
!
    goto 1000
!
999 continue
! -----------------------------------------------
! --- NOMBRE DE TENTATIVES DE RELANCE DEPASSE ---
! -----------------------------------------------
    if (intg .gt. 5) then
! --- ON REGARDE SI L'ETAT INITIAL MATERIAU AVAIT SOLLICITE
! --- UN MECANISME DE TRACTION : ETAT INIT = SIGD, VINS
        iret = 3
        noconv = .true.
        do i = 1, ndi
            call hujprj(i, ydt, dev, pf, qf)
            if (((pf+deux*rtrac-ptrac)/abs(pref)) .gt. -r8prem()) then
                noconv = .false.
                iret = 0
            end if
        end do
        if (.not. noconv) then
! --- EN POSANT NOCONV = .TRUE., ON CONDUIT L'ALGORITHME PRESENT
! --- A IMPOSER UN ETAT DE CONTRAINTES ISOTROPE COMMUN
! --- AUX 3 SEUILS PLASTIQUES DE TRACTION
            do i = 1, nr
                dy(i) = zero
            end do
            vind(1:nvi) = vins(1:nvi)
            vinf(1:nvi) = vins(1:nvi)
            iret = 0
        else
!
! -------------------------------------------------------
! --- SI PROCHE TRACTION ET NON CONVERGENCE, ON IMPOSE
! --- ETAT DE CONTRAINTES ISOTROPE
! -------------------------------------------------------
            impose = .false.
            do i = 1, ndi
                call hujprj(i, yft, dev, pf, qf)
                if ((abs(pf-ptrac)/abs(pref)) .lt. cinq*rtrac/abs(pref)) then
                    impose = .true.
                end if
            end do
            if (impose) then
                noconv = .false.
                do i = 1, 3
                    vind(23+i) = zero
                    vind(27+i) = zero
                    vind(5+4*i) = zero
                    vind(6+4*i) = zero
                    vind(7+4*i) = zero
                    vind(8+4*i) = zero
                end do
                do i = 1, ndi
                    dy(i) = -yd(i)-deux*rtrac/e0
                    dy(ndi+i) = -yd(ndi+i)
                end do
                do i = 1, nbmeca+nbmect+1
                    dy(ndt+i) = zero
                end do
                vinf(1:nvi) = vind(1:nvi)
                iret = 0
            end if
        end if
!
        goto 1000
    end if
!
! -------------------------------------------
! --- GESTION DES NON-CONVERGENCE LOCALES ---
! -------------------------------------------
! --- Y AVAIT IL UN MECANISME CYCLIQUE DEJA DESACTIVE
!     DURANT CETTE TENTATIVE?
    msup(1) = 0
    msup(2) = 0
    msup(3) = 0
    msup(4) = 0
    jj = 0
    do i = 5, 8
        if ((vind(23+i) .ne. vins(23+i)) .and. (vind(23+i) .eq. zero)) then
            jj = jj+1
            msup(jj) = i
        end if
    end do
!
! --- EXISTE-T-IL UN PB DE TANGENCE ENTRE MECANISMES
    do k = 1, 4
        prox(k) = .false.
        proxc(k) = .false.
    end do
!
    do i = 1, 22
        matert(i, 1) = mater(i, 1)
        matert(i, 2) = mater(i, 2)
    end do
!
    do k = 1, nbmeca
        if ((indi(k) .gt. 4) .and. (indi(k) .lt. 8)) then
            kk = indi(k)-4
            call hujpxd(indi(k), matert, yft, vind, prox(kk), &
                        proxc(kk))
        end if
    end do
!
    probt = .false.
    do i = 1, 3
        prob(i) = zero
        if (prox(i)) then
            prob(i) = un
            probt = .true.
        else if (proxc(i)) then
            prob(i) = deux
            probt = .true.
        end if
    end do
!
    if (probt) then
        vind(1:nvi) = vins(1:nvi)
        do i = 1, 3
            if (prob(i) .eq. un) then
                vind(i+4) = mater(18, 2)
                vind(23+i) = un
                vind(27+i) = zero
                vind(4*i+5) = zero
                vind(4*i+6) = zero
                vind(4*i+7) = zero
                vind(4*i+8) = zero
                vind(5*i+31) = zero
                vind(5*i+32) = zero
                vind(5*i+33) = zero
                vind(5*i+34) = zero
                vind(5*i+35) = mater(18, 2)
            else if (prob(i) .eq. deux) then
                vind(27+i) = zero
            end if
        end do
        iret = 0
        probt = .false.
!
! --- MECANISME CYCLIQUE A DESACTIVE
! --- ET DEJA DESACTIVE ANTERIEUREMENT
        if (jj .ne. 0) then
            do i = 1, jj
                vind(23+msup(i)) = zero
            end do
        end if
!
        vinf(1:nvi) = vind(1:nvi)
        iret = 2
        goto 1000
    end if
!
    if (tracti) then
        vind(1:nvi) = vins(1:nvi)
        modif = .false.
        do i = 1, nbmect
            if (yet(ndt+1+nbmeca+i) .eq. zero) then
                modif = .true.
                if (indi(i) .le. 8) then
                    if (indi(i) .lt. 5) then
                        if ((abs(vind(4*indi(i)+5)) .gt. r8prem()) .or. &
                            (abs(vind(4*indi(i)+6)) .gt. r8prem())) then
                            vind(23+indi(i)) = -un
                        else
                            vind(23+indi(i)) = zero
                        end if
                    else
                        vind(23+indi(i)) = zero
                    end if
                else
                    bnews(indi(i)-8) = .true.
                    neglam(indi(i)-8) = .true.
                end if
            end if
        end do
!
        do i = 1, nbmect
            if (indi(i) .eq. 8) then
                vind(23+indi(i)) = zero
                modif = .true.
            end if
        end do
!
        mtrac = .false.
        do i = 1, 3
! --- ON NE DOIT PAS REACTIVE UN MECANISME DE TRACTION QUI DONNE
!     COMME PREDICTEUR UN MULTIPLICATEUR PLASTIQUE NEGATIF
            if (.not. neglam(i)) then
                call hujprj(i, yft, dev, pf, qf)
! ----------------------------------------------------
! ---> ACTIVATION MECANISMES DE TRACTION NECESSAIRES
! ----------------------------------------------------
                if (((pf+deux*rtrac-ptrac)/abs(pref)) .gt. -r8prem()) then
                    bnews(i) = .false.
                    if (.not. modif) mtrac = .true.
                end if
            end if
        end do
        vinf(1:nvi) = vind(1:nvi)
        iret = 2
        goto 1000
    end if
!
!-----------------------------------------------------------
! --- ESSAIS HEURISTIQUES POUR RELANCER LA RESOLUTION LOCALE
!-----------------------------------------------------------
    maxi = zero
    resi = 0
    do i = 1, nr
        if (abs(r(i)) .gt. maxi) then
            maxi = abs(r(i))
            resi = i
        end if
    end do
    cycl = .false.
    do i = 1, nbmeca
        if ((indi(i) .gt. 4) .and. (indi(i) .lt. 8) .and. (vind(indi(i)) .eq. mater(18, 2))) then
            cycl = .true.
        end if
    end do
!
! ---------------------------------------------------------------
! --- SI RESIDU LOCAL MAXI PORTE PAR RDEV_CYC => MECANISME RETIRE
! ---------------------------------------------------------------
!
    if ((resi .gt. 7) .and. (resi .le. 7+nbmeca)) then
        resi = resi-7
        if ((indi(resi) .gt. 4) .and. (indi(resi) .lt. 8)) then
!
            vind(1:nvi) = vins(1:nvi)
            vind(23+indi(resi)) = zero
            if (jj .ne. 0) then
                do i = 1, jj
                    vind(23+msup(i)) = zero
                end do
            end if
!
! --- EXISTE-T-IL UN MECANISME DEVIATOIRE AYANT LE MEME COMPORTEMENT
!     QUE CELUI IDENTIFIE PRECEDEMMENT COMME POSANT PROBLEME ?
            do i = 1, nbmeca
                if ((indi(i) .gt. 4) .and. (indi(i) .lt. 8) .and. &
                    (((maxi-abs(r(7+i)))/toler) .lt. toler) .and. (i .ne. resi)) then
                    vind(23+indi(i)) = zero
                end if
            end do
!
            iret = 2
            vinf(1:nvi) = vind(1:nvi)
            goto 1000
        else
            iret = 3
        end if
    end if
!
! ---------------------------------------------------------------
! --- SI MECA CYCLIQUE ALORS ILS SONT RETIRES
! ---------------------------------------------------------------
!
    if (cycl) then
        vind(1:nvi) = vins(1:nvi)
        do i = 1, nbmeca
           if ((indi(i) .gt. 4) .and. (indi(i) .lt. 8) .and. (vind(indi(i)) .eq. mater(18, 2))) then
                vind(23+indi(i)) = zero
            end if
        end do
        iret = 2
        vinf(1:nvi) = vind(1:nvi)
        goto 1000
    end if
!
! ---------------------------------------------------------------
! --- SI MECANISME TRACTION ACTIF => RETIRE DE MPOT
! ---------------------------------------------------------------
!
    if (nbmect .ne. nbmeca) then
        vind(1:nvi) = vins(1:nvi)
        iret = 2
        do i = nbmeca+1, nbmect
            if (yet(ndt+1+nbmeca+i) .eq. zero) then
                bnews(indi(i)-8) = .true.
            end if
        end do
        vinf(1:nvi) = vind(1:nvi)
        goto 1000
    end if
!
! ---------------------------------------------------------------
! --- CONTROLE DU PREDICTEUR ELASTIQUE: YE(LAMBDA)
! ---------------------------------------------------------------
!
    vind(1:nvi) = vins(1:nvi)
    euler = .true.
    lamin = 1.d2
    imin = 0
    do i = 1, nbmeca
        if (yet(ndt+1+nbmeca+i) .eq. zero) then
            if ((indi(i) .gt. 4) .and. (indi(i) .lt. 9)) then
                vind(indi(i)+23) = 0
                euler = .false.
            else if (indi(i) .lt. 5) then
                if ((abs(vind(4*indi(i)+5)) .gt. r8prem()) .or. &
                    (abs(vind(4*indi(i)+6)) .gt. r8prem())) then
                    vind(23+indi(i)) = -un
                else
                    vind(23+indi(i)) = zero
                end if
                euler = .false.
            end if
        else if (yet(ndt+1+nbmeca+i) .lt. lamin) then
            lamin = yet(ndt+1+nbmeca+i)
            imin = i
        end if
    end do
!
    if (.not. euler) then
! --- MECANISME CYCLIQUE A DESACTIVE
! --- ET DEJA DESACTIVE ANTERIEUREMENT
        if (jj .ne. 0) then
            do i = 1, jj
                vind(23+msup(i)) = zero
            end do
        end if
!
        vinf(1:nvi) = vind(1:nvi)
        iret = 2
        goto 1000
    else if (imin .gt. 0) then
        if (indi(imin) .lt. 5) then
            vind(23+indi(imin)) = -un
        else
            vind(23+indi(imin)) = zero
        end if
        vinf(1:nvi) = vind(1:nvi)
        iret = 2
        goto 1000
    end if
!
! ---------------------------------------------------------------
! --- DERNIER ESSAI: VALEUR DES CONTRAINTES PRE, DURANT ET POST
! ---------------------------------------------------------------
    ltry = .false.
    do i = 1, ndi
        call hujprj(i, yd, dev, pf, qf)
        if (((pf+deux*rtrac-ptrac)/abs(pref)) .gt. -r8prem()) then
            bnews(i) = .false.
            ltry = .true.
        end if
        call hujprj(i, yet, dev, pf, qf)
        if (((pf+deux*rtrac-ptrac)/abs(pref)) .gt. -r8prem()) then
            bnews(i) = .false.
            ltry = .true.
        end if
        call hujprj(i, yft, dev, pf, qf)
        if (((pf+deux*rtrac-ptrac)/abs(pref)) .gt. -r8prem()) then
            bnews(i) = .false.
            ltry = .true.
        end if
    end do
!
    if (ltry) then
        vinf(1:nvi) = vind(1:nvi)
        iret = 2
        goto 1000
    else
        iret = 3
    end if
!
1000 continue
end subroutine
