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
subroutine mpinvr(nbmesu, nbmode, nbabs, phi, rmesu, &
                  coef, xabs, lfonct, reta)
!     PROJ_MESU_MODAL : RESOLUTION DU SYSTEME PAR SVD OU PAR LU
!                       (DONNEES FREQUENTIELLES REELLES)
!
!     IN  : NBMESU : NOMBRE DE MESURE
!     IN  : NBMODE : NOMBRE DE MODES
!     IN  : NBABS  : NOMBRE D ABSCISSES
!     IN  : PHI    : BASE DE PROJECTION
!     IN  : RMESU  : VALEURS MESUREES
!     IN  : COEF   : COEFFICIENTS DE PONDERATION
!     IN  : XABS   : LISTE REELLE D ABSCISSES
!     IN  : LFONCT : .TRUE. SI COEFFICIENTS DE PONDERATION
!                    DEPENDENT DE LA FREQUENCE
!     OUT : RETA   : DEPLACEMENT GENERALISE   ( MATRICE )
!
    implicit none
!
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/mtcrog.h"
#include "asterfort/rslsvd.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nbmesu, nbmode, nbabs
    integer(kind=8) :: vali
    real(kind=8) :: phi(nbmesu, nbmode), xabs(nbabs), coef(*)
    real(kind=8) :: valr
    real(kind=8) :: rmesu(nbmesu, nbabs), reta(nbmode, nbabs)
    aster_logical :: lfonct
!
    integer(kind=8) :: imod, jmod, imes, iabs, ierr, ibid, jmes
    integer(kind=8) :: lsecmb, lwks, lphiph, lphitp, lmatsy, lwork, leta, lvals, lu, lv
    real(kind=8) :: alpha, eps
    aster_logical :: nul
    character(len=3) :: method
    character(len=8) :: regul
    character(len=16) :: nomcha
    character(len=16) :: secmb, wks, phipht, phitph, matsys, work, eta, vals, u
    character(len=16) :: v
!
! ----------------------------------------------------------------------
!
    call jemarq()
!
! CREATION DES VECTEURS DE TRAVAIL
    secmb = '&SECMB'
    wks = '&WKS'
    phipht = '&PHIPHT'
    phitph = '&PHITPH'
    matsys = '&MATSYS'
    work = '&WORK'
    eta = '&ETA'
    vals = '&VALS'
    u = '&U'
    v = '&V'
!
    call wkvect(secmb, 'V V R', nbmode, lsecmb)
    call wkvect(wks, 'V V R', nbmode, lwks)
    call wkvect(phipht, 'V V R', nbmesu*nbmesu, lphiph)
    call wkvect(phitph, 'V V R', nbmode*nbmode, lphitp)
    call wkvect(matsys, 'V V R', nbmode*nbmode, lmatsy)
    call wkvect(work, 'V V R', nbmode, lwork)
    call wkvect(eta, 'V V R', nbmode, leta)
    call wkvect(vals, 'V V R', nbmode, lvals)
    call wkvect(u, 'V V R', nbmode*nbmode, lu)
    call wkvect(v, 'V V R', nbmode*nbmode, lv)
!
! METHODE DE RESOLUTION : LU / SVD
    call getvtx('RESOLUTION', 'METHODE', iocc=1, scal=method, nbret=ibid)
    if (ibid .eq. 0) method = 'LU'
!
    if (method .eq. 'SVD') then
        call getvr8('RESOLUTION', 'EPS', iocc=1, scal=eps, nbret=ibid)
        if (ibid .eq. 0) eps = 0.d0
    end if
!
! REGULARISATION : NON / NORM_MIN / TIK_RELA
    call getvtx('RESOLUTION', 'REGUL', iocc=1, scal=regul, nbret=ibid)
    if (ibid .eq. 0) regul = 'NON'
!
    call getvtx('MODELE_MESURE', 'NOM_CHAM', iocc=1, scal=nomcha, nbret=ibid)
!
! ===============================
! CALCUL DE PHI_TRANSPOSEE * PHI
! ===============================
    do imod = 1, nbmode
        do jmod = 1, nbmode
            zr(lphitp-1+imod+nbmode*(jmod-1)) = 0.d0
            do imes = 1, nbmesu
                zr(lphitp-1+imod+nbmode*(jmod-1)) = zr( &
                                                    lphitp-1+imod+nbmode*(jmod-1))+phi(imes, &
                                                                              imod)*phi(imes, jmod &
                                                                                                 )
            end do
        end do
    end do
!
    if (nbmesu .lt. nbmode) then
! ===============================
! CALCUL DE PHI * PHI_TRANSPOSEE
! ===============================
        do imes = 1, nbmesu
            do jmes = 1, nbmesu
                zr(lphiph-1+imes+nbmesu*(jmes-1)) = 0.d0
                do imod = 1, nbmode
                    zr(lphiph-1+imes+nbmesu*(jmes-1)) = zr( &
                                                        lphiph-1+imes+nbmesu*(jmes-1) &
                                                        )+phi(imes, &
                                                              imod)*phi(jmes, imod &
                                                                        )
                end do
            end do
        end do
    end if
!
! =======================================
! CALCUL DE LA REPONSE GENERALISEE : RETA
! =======================================
!
! DEBUT DE LA BOUCLE SUR LES ABSCISSES (FREQUENCE)
! *****************************
    do iabs = 1, nbabs
!
        nul = .true.
!
! DEBUT DE LA BOUCLE SUR LES MODES
! ********************************
        do imod = 1, nbmode
!
! RECHERCHE DU COEFFICIENT DE PONDERATION
! ***************************************
            if (lfonct) then
!             -> ALPHA DEPENDANT DES ABSCISSES
                alpha = coef(nbmode*(iabs-1)+imod)
            else
!             -> ALPHA INDEPENDANT DES ABSCISSES
                alpha = coef(imod)
            end if
!
!         -> ON VERIFIE QUE ALPHA > 0, SINON ARRET
            if (alpha .lt. 0.d0) then
                vali = iabs
                call utmess('F', 'ALGORITH15_24', si=vali)
            else if (alpha .gt. r8prem()) then
                nul = .false.
            end if
!
! DETERMINATION DE LA MATRICE A INVERSER :
! MATSYS(IABS) = PHI_T*PHI + ALPHA(IABS)
! ****************************************
            do jmod = 1, nbmode
                zr(lmatsy-1+imod+nbmode*(jmod-1)) = zr(lphitp-1+imod+nbmode*(jmod-1))
            end do
!
            zr(lmatsy-1+imod+nbmode*(imod-1)) = zr(lmatsy-1+imod+nbmode*(imod-1))+alpha
!
! DETERMINATION DU SECOND MEMBRE :
! SCDMB(IABS) = PHI_T*Q + ALPHA(IABS)*RETA(IABS-1)
! RQ : A IABS=1, RETA(0)=0 (LA SOLUTION A PRIORI EST NULLE)
! ********************************************************
            zr(lsecmb-1+imod) = 0.d0
!
            do imes = 1, nbmesu
                zr(lsecmb-1+imod) = zr(lsecmb-1+imod)+phi(imes, imod)*rmesu(imes, iabs)
            end do
!
            if ((regul .eq. 'TIK_RELA') .and. (iabs .gt. 1)) then
                ibid = iabs-1
                zr(lsecmb-1+imod) = zr(lsecmb-1+imod)+alpha*reta(imod, ibid)
            end if
!
! FIN DE LA BOUCLE SUR LES MODES
! ******************************
        end do
!
!
! RESOLUTION DU SYSTEME :
! MATSYS(IABS) * ETA(IABS) = SCDMB(IABS)
! **************************************
!       -> ALARME SI ALPHA NUL ET NBMESU<NBMODE : MOINDRE NORME
        if ((nbmesu .lt. nbmode) .and. (nul)) then
            call utmess('A', 'ALGORITH15_25')
!
            if (regul .eq. 'NON') then
! CALCUL MOINDRE NORME
                call utmess('A', 'ALGORITH15_26')
                do imes = 1, nbmesu
                    zr(lsecmb-1+imes) = rmesu(imes, iabs)
                    do jmes = 1, nbmesu
                        zr(lmatsy-1+imes+nbmode*(jmes-1)) = zr(lphiph-1+imes+nbmesu*(jmes-1))
                    end do
                end do
!
! CHOIX POUR LA METHODE D INVERSION
                if (method .eq. 'SVD') then
! METHODE SVD
! CREATION DU VECTEUR SECOND MEMBRE
                    do jmes = 1, nbmesu
                        zr(leta-1+jmes) = zr(lsecmb-1+jmes)
                    end do
!
                    call rslsvd(nbmode, nbmesu, nbmesu, zr(lmatsy), zr(lvals), &
                                zr(lu), zr(lv), 1, zr(leta), eps, &
                                ierr, zr(lwork))
                    if (ierr .ne. 0) then
                        vali = iabs
                        valr = xabs(iabs)
                        call utmess('F', 'ALGORITH15_27', si=vali, sr=valr)
                    end if
!
                else
! METHODE DE CROUT
                    call mtcrog(zr(lmatsy), zr(lsecmb), nbmode, nbmesu, 1, &
                                zr(leta), zr(lwks), ierr)
                    if (ierr .ne. 0) then
                        vali = iabs
                        valr = xabs(iabs)
                        call utmess('F', 'ALGORITH15_28', si=vali, sr=valr)
                    end if
                end if
!
! COPIE DES RESULTATS DANS RETA
                do jmod = 1, nbmode
                    do jmes = 1, nbmesu
                        reta(jmod, iabs) = phi(jmes, jmod)*zr(leta-1+jmes)
                    end do
                end do
!
                goto 100
            end if
        end if
! FIN CALCUL MOINDRE NORME
!
!
! RESOLUTION EN FONCTION DE LA METHODE DE RESOLUTION
        if (method .eq. 'SVD') then
! METHODE SVD
! CREATION DU VECTEUR SECOND MEMBRE
            do jmod = 1, nbmode
                zr(leta-1+jmod) = zr(lsecmb-1+jmod)
            end do
!
            call rslsvd(nbmode, nbmode, nbmode, zr(lmatsy), zr(lvals), &
                        zr(lu), zr(lv), 1, zr(leta), eps, &
                        ierr, zr(lwork))
            if (ierr .ne. 0) then
                vali = iabs
                valr = xabs(iabs)
                call utmess('F', 'ALGORITH15_27', si=vali, sr=valr)
            end if
!
        else
! METHODE DE CROUT
            call mtcrog(zr(lmatsy), zr(lsecmb), nbmode, nbmode, 1, &
                        zr(leta), zr(lwks), ierr)
            if (ierr .ne. 0) then
                vali = iabs
                valr = xabs(iabs)
                call utmess('F', 'ALGORITH15_28', si=vali, sr=valr)
            end if
        end if
!
! RECUPERATION DES RESULTATS
        do jmod = 1, nbmode
            reta(jmod, iabs) = zr(leta-1+jmod)
        end do
!
! FIN DE LA BOUCLE SUR LES ABSCISSES (FREQUENCE)
! ***************************
100     continue
    end do
!
!
! DESTRUCTION DES VECTEURS DE TRAVAIL
!
    call jedetr(secmb)
    call jedetr(wks)
    call jedetr(phipht)
    call jedetr(phitph)
    call jedetr(matsys)
    call jedetr(work)
    call jedetr(eta)
    call jedetr(vals)
    call jedetr(u)
    call jedetr(v)
!
    call jedema()
end subroutine
