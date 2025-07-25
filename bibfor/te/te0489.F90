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
subroutine te0489(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/dchapg.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/norsig.h"
#include "asterfort/radipg.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!     BUT:
!       CALCUL DES INDICATEURS LOCAUX DE DECHARGE
!                 ET DE PERTE DE RADIALITE POUR LES NORMES :
!              VMIS      : SECOND INVARIANT DU TENSEUR DES CONTRAINTES
!                          DEVIATORIQUES
!              TOTAL     : SECOND INVARIANT DU TENSEUR DES CONTRAINTES
!              VMIS_CINE : SECOND INVARIANT DU DEVIATEUR DU TENSEUR
!                          SIGMA - X
!                          OU SIGMA EST LE TENSEUR DES CONTRAINTES
!                          ET X     EST LE TENSEUR DE RAPPEL
!              TOTAL_CINE: SECOND INVARIANT DU TENSEUR SIGMA - X
!    ON NOTE SIGMA2 = SIGMA(M,T+DT)
!            SIGMA1 = SIGMA(M,T)
!            DSIGMA = SIGMA2 - SIGMA1
!    A)LES INDICATEURS LOCAUX DE DECHARGE :
!      I = (NORME(SIGMA2) - NORME(SIGMA1)/NORME(SIGMA2)
!               SONT CALCULES :
!    B)LES INDICATEURS LOCAUX DE PERTE DE RADIALITE :
!      I = 1 - ABS(DSIGMA : SIGMA1)/(NORME(DSIGMA)*NORME(SIGMA1))
!               SONT CALCULES :
!          ELEMENTS ISOPARAMETRIQUES 2D ET 3D
!          OPTIONS : 'DERA_ELGA'
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
! ......................................................................
!
!
!
    integer(kind=8) :: mxcmel
    parameter(mxcmel=162)
    integer(kind=8) :: nbpgmx
    parameter(nbpgmx=27)
!
    integer(kind=8) :: jtab(7)
    integer(kind=8) :: iret, npg1, ipoids, ivf, idfde, jgano, icompo
    integer(kind=8) :: idera1, idera2
    integer(kind=8) :: nbcmp, imate, idecal, isig
    integer(kind=8) :: nno, nbsig, nnos, isigtm, isigtp, idvar1
    integer(kind=8) :: idvar2, nbvarint, npg, i, k, ndim, igau, icodre(3)
    integer(kind=8) :: nbvari, jprolp, jvalep, nbvalp, ibid
    real(kind=8) :: sigma1(mxcmel), sigma2(mxcmel)
    real(kind=8) :: sigt1(mxcmel), sigt2(mxcmel)
    real(kind=8) :: dchav(nbpgmx), dchat(nbpgmx)
    real(kind=8) :: dchax(nbpgmx), dchay(nbpgmx)
    real(kind=8) :: radiv(nbpgmx), radit(nbpgmx)
    real(kind=8) :: cosang(mxcmel)
    real(kind=8) :: xrapel(nbpgmx*6), sigmx(6)
    real(kind=8) :: pm, pp, dp, rp, rp0
    real(kind=8) :: trx1, trx2
    real(kind=8) :: x1(mxcmel), x2(mxcmel)
    real(kind=8) :: trsig1, trsig2
    real(kind=8) :: e, dsde, sigy, alfafa, coco, unsurn
    real(kind=8) :: zero, untier, cst1, tp, para_vale
    real(kind=8) :: zernor, dchaxm
    real(kind=8) :: valres(3)
    real(kind=8) :: trace, sqrt
    real(kind=8) :: sigeqn
    character(len=4) :: fami
    character(len=8) :: para_type
    character(len=16) :: nomres(3)
    character(len=16) :: compor
! ----------------------------------------------------------------------
!
! ---- INITIALISATIONS :
!      ---------------
    zero = 0.0d0
    untier = 1.0d0/3.0d0
    zernor = 1000.0d0*r8prem()
    fami = 'RIGI'
    nbvalp = 0
    ibid = 0
!
    do i = 1, mxcmel
        sigt1(i) = zero
        sigt2(i) = zero
        sigma1(i) = zero
        sigma2(i) = zero
        x1(i) = zero
        x2(i) = zero
        xrapel(i) = zero
    end do
!
    do i = 1, nbpgmx
        dchav(i) = zero
        dchat(i) = zero
        dchax(i) = zero
        dchay(i) = zero
        radiv(i) = zero
        radit(i) = zero
        cosang(i) = zero
    end do
!
! ----     DIMENSION DE L'ELEMENT :
! ----     NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT :
    nbsig = nbsigm()
!
! ---- RECUPERATION DU NOMBRE DE COMPOSANTES
    call tecach('OOO', 'PDERAMG', 'L', iret, nval=7, &
                itab=jtab)
    idera1 = jtab(1)
    nbcmp = jtab(2)/jtab(3)
!
! ---- RECUPERATION DE LA CARTE DE COMPORTEMENT :
!      -----------------------------------------------------
    call jevech('PCOMPOR', 'L', icompo)
    compor = zk16(icompo)
!
! ---- RECUPERATION DES CONTRAINTES A L'INSTANT T :
!
    call tecach('OOO', 'PCONTMR', 'L', iret, nval=3, &
                itab=jtab)
    npg = jtab(3)
    isigtm = jtab(1)
!
! ---- RECUPERATION DES CONTRAINTES A L'INSTANT T+DT :
!
    call jevech('PCONTPR', 'L', isigtp)
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    ASSERT(npg .eq. npg1)
!
! ---- AFFECTATION DES VECTEURS DE TRAVAIL SIGMA1 ET SIGMA2
! ---- REPRESENTANT LES TENSEURS DE CONTRAINTES RESPECTIVEMENT
! ---- AUX INSTANTS T ET T+DT :
!      ----------------------
    k = 0
!
    do igau = 1, npg
        do i = 1, nbsig
            k = k+1
            sigt1(i+(igau-1)*nbsig) = zr(isigtm+k-1)
            sigt2(i+(igau-1)*nbsig) = zr(isigtp+k-1)
            sigma1(i+(igau-1)*nbsig) = zr(isigtm+k-1)
            sigma2(i+(igau-1)*nbsig) = zr(isigtp+k-1)
        end do
    end do
!
    if (compor(1:14) .ne. 'VMIS_CINE_LINE') then
!
! ---- CAS DE LA NORME VMIS : CALCUL DES DEVIATEURS DES CONTRAINTES
! ---- SIGMA1 ET SIGMA2 :
!      ----------------
!
        do igau = 1, npg
!
            trsig1 = sigma1( &
                     1+(igau-1)*nbsig)+sigma1(2+(igau-1)*nbsig)+sigma1(3+(igau-1)*nbsig)
!
            trsig2 = sigma2( &
                     1+(igau-1)*nbsig)+sigma2(2+(igau-1)*nbsig)+sigma2(3+(igau-1)*nbsig)
            sigma1(1+(igau-1)*nbsig) = sigma1(1+(igau-1)*nbsig)-untier*trsig1
            sigma1(2+(igau-1)*nbsig) = sigma1(2+(igau-1)*nbsig)-untier*trsig1
            sigma1(3+(igau-1)*nbsig) = sigma1(3+(igau-1)*nbsig)-untier*trsig1
!
            sigma2(1+(igau-1)*nbsig) = sigma2(1+(igau-1)*nbsig)-untier*trsig2
            sigma2(2+(igau-1)*nbsig) = sigma2(2+(igau-1)*nbsig)-untier*trsig2
            sigma2(3+(igau-1)*nbsig) = sigma2(3+(igau-1)*nbsig)-untier*trsig2
        end do
!
! ---- DANS LE CAS D'UN ECROUISSAGE CINEMATIQUE
! ---- RECUPERATION DES COMPOSANTES DU TENSEUR DE RAPPEL
! ---- AUX INSTANTS T ET T+DT :
!      ----------------------
    else if (compor(1:14) .eq. 'VMIS_CINE_LINE') then
!
        call jevech('PVARIMR', 'L', idvar1)
        call jevech('PVARIPR', 'L', idvar2)
!
! ----   AFFECTATION DES VECTEURS DE TRAVAIL X1 ET X2 REPRESENTANT
! ----   LES TENSEURS DE RAPPEL RESPECTIVEMENT AUX INSTANTS T
! ----   ET T+DT :
!        -------
        nbvarint = 7
!
        do igau = 1, npg
            do i = 1, nbsig
                x1(i+(igau-1)*nbsig) = zr(idvar1+i+(igau-1)*nbvarint-1)
                x2(i+(igau-1)*nbsig) = zr(idvar2+i+(igau-1)*nbvarint-1)
            end do
        end do
!
! ---- CAS DE LA NORME TOTAL_CINE : CALCUL DES CONTRAINTES
! ---- (SIGMA1-X1) ET (SIGMA2-X2) :
!      -------------------------
!
        k = 0
        do igau = 1, npg
            do i = 1, nbsig
                k = k+1
!
                sigt1(i+(igau-1)*nbsig) = zr(isigtm+k-1)
!     &                                -X1(I+ (IGAU-1)*NBSIG)
                sigt2(i+(igau-1)*nbsig) = zr(isigtp+k-1)
!     &                                -X2(I+ (IGAU-1)*NBSIG)
                sigma1(i+(igau-1)*nbsig) = zr(isigtm+k-1)-x1(i+(igau-1)*nbsig)
                sigma2(i+(igau-1)*nbsig) = zr(isigtp+k-1)-x2(i+(igau-1)*nbsig)
            end do
        end do
!
! ---- CAS DE LA NORME VMIS_CINE : CALCUL DES DEVIATEURS DES CONTRAINTES
! ---- (SIGMA1-X1) ET (SIGMA2-X2) :
!      -------------------------
!
        do igau = 1, npg
!
            trsig1 = sigma1( &
                     1+(igau-1)*nbsig)+sigma1(2+(igau-1)*nbsig)+sigma1(3+(igau-1)*nbsig)
!
            trsig2 = sigma2( &
                     1+(igau-1)*nbsig)+sigma2(2+(igau-1)*nbsig)+sigma2(3+(igau-1)*nbsig)
!
            trx1 = x1(1+(igau-1)*nbsig)+x1(2+(igau-1)*nbsig)+x1(3+(igau-1)*nbsig)
!
            trx2 = x2(1+(igau-1)*nbsig)+x2(2+(igau-1)*nbsig)+x2(3+(igau-1)*nbsig)
!
            sigma1(1+(igau-1)*nbsig) = sigma1(1+(igau-1)*nbsig)-untier*(trsig1-trx1)
            sigma1(2+(igau-1)*nbsig) = sigma1(2+(igau-1)*nbsig)-untier*(trsig1-trx1)
            sigma1(3+(igau-1)*nbsig) = sigma1(3+(igau-1)*nbsig)-untier*(trsig1-trx1)
            sigma1(4+(igau-1)*nbsig) = sigma1(4+(igau-1)*nbsig)
!
            sigma2(1+(igau-1)*nbsig) = sigma2(1+(igau-1)*nbsig)-untier*(trsig2-trx2)
            sigma2(2+(igau-1)*nbsig) = sigma2(2+(igau-1)*nbsig)-untier*(trsig2-trx2)
            sigma2(3+(igau-1)*nbsig) = sigma2(3+(igau-1)*nbsig)-untier*(trsig2-trx2)
            sigma2(4+(igau-1)*nbsig) = sigma2(4+(igau-1)*nbsig)
!
            if (ndim .eq. 3) then
                sigma1(5+(igau-1)*nbsig) = sigma1(5+(igau-1)*nbsig)
                sigma1(6+(igau-1)*nbsig) = sigma1(6+(igau-1)*nbsig)
!
                sigma2(5+(igau-1)*nbsig) = sigma2(5+(igau-1)*nbsig)
                sigma2(6+(igau-1)*nbsig) = sigma2(6+(igau-1)*nbsig)
            end if
        end do
!
    end if
!
! ---- CALCUL DE L'INDICATEUR LOCAL DE DECHARGE DCHA :
! ----  I = (NORME(SIGMA2) - NORME(SIGMA1))/NORME(SIGMA2) :
!      --------------------------------------------------
    call dchapg(sigt1, sigt2, npg, nbsig, dchat)
    call dchapg(sigma1, sigma2, npg, nbsig, dchav)
!
!
! ---- CALCUL DE L'INDICATEUR LOCAL DE PERTE DE RADIALITE ERR_RADI
    read (zk16(icompo-1+2), '(I16)') nbvari
    call jevech('PVARIMR', 'L', idvar1)
    call jevech('PVARIPR', 'L', idvar2)
    call jevech('PMATERC', 'L', imate)
    call radipg(sigt1, sigt2, npg, nbsig, radit, &
                cosang, 1, compor, zi(imate), nbvari, &
                zr(idvar1), zr(idvar2))
!
    if (compor(1:9) .eq. 'VMIS_ISOT') then
!
! - ON RECUPERE
!      > LES INDICATEURS A L'INSTANT T-DT
!      > LES VARIABLES INTERNES AUX INSTANTS T et T+DT
!      > LES PARAMETERS MATERIAUX A INSTANT T
!
!
! ---- CALCUL DE L'INDICATEUR LOCAL DE PERTE DE RADIALITE RADI:
! ----  I = 1- ABS(SIGMA1:DSIGMA)/(NORME(SIGMA1)*NORME(DSIGMA) :
!      -------------------------------------------------------
        call radipg(sigma1, sigma2, npg, nbsig, radiv, &
                    cosang, 0, compor, imate, nbvari, &
                    zr(idvar1), zr(idvar2))
!
! ---- CALCUL DE L'INDICATEUR LOCAL IND_DCHA et VAL_DCHA
! ---- CALCUL DU TENSEUR DE RAPPEL X
!      -------------------------------------------------------
!
! --- ON NE TRAITE QUE LES LOIS SUIVANTES
!        - VMIS_ISOT_LINE
!        - VMIS_ISOT_TRAC
!        - VMIS_ISOT_PUIS
!
!
! --- BOUCLE SUR LES POINTS DE GAUSS
!
        do igau = 1, npg
            idecal = (igau-1)*6
! --- DEFORMATION PLASTIQUE CUMULEE
            pm = zr(idvar1-1+(igau-1)*nbvari+1)
            pp = zr(idvar2-1+(igau-1)*nbvari+1)
            dp = pp-pm
            dchaxm = zr(idera1-1+(igau-1)*nbcmp+3)
!
            if (pm .le. zernor) then
                dchax(igau) = 1.d0
            else if (abs(dchaxm+2.d0) .gt. zernor) then
                if (dp .gt. zernor) then
                    if (cosang(igau) .gt. zernor) then
                        dchax(igau) = 2.d0
                    else
                        dchax(igau) = -2.d0
                    end if
                else if (dp .le. zernor) then
                    if (abs(dchaxm+1.d0) .gt. zernor) then
! --- RECUPERATION DES CARACTERISTIQUES DE LA LOI DE COMPORTEMENT
                        if (compor .eq. 'VMIS_ISOT_LINE') then
                            nomres(1) = 'D_SIGM_EPSI'
                            nomres(2) = 'SY'
                            call rcvalb(fami, igau, 1, '+', zi(imate), &
                                        ' ', 'ECRO_LINE', ibid, ' ', [0.d0], &
                                        2, nomres, valres, icodre, 2)
                            dsde = valres(1)
                            rp0 = valres(2)
!
                            nomres(1) = 'E'
                            call rcvalb(fami, igau, 1, '+', zi(imate), &
                                        ' ', 'ELAS', ibid, ' ', [0.d0], &
                                        1, nomres, valres, icodre, 2)
                            e = valres(1)
                            rp = dsde*e/(e-dsde)*pm+rp0
!
                        else if (compor(10:14) .eq. '_PUIS') then
                            nomres(1) = 'E'
                            call rcvalb(fami, igau, 1, '+', zi(imate), &
                                        ' ', 'ELAS', ibid, ' ', [0.d0], &
                                        1, nomres, valres, icodre, 2)
                            e = valres(1)
!
                            nomres(1) = 'SY'
                            nomres(2) = 'A_PUIS'
                            nomres(3) = 'N_PUIS'
                            call rcvalb(fami, igau, 1, '+', zi(imate), &
                                        ' ', 'ECRO_PUIS', ibid, ' ', [0.d0], &
                                        3, nomres, valres, icodre, 2)
                            sigy = valres(1)
                            alfafa = valres(2)
                            coco = e/alfafa/sigy
                            unsurn = 1.d0/valres(3)
                            rp = sigy*(coco*pp)**unsurn+sigy
!
                        else if (compor(10:14) .eq. '_TRAC') then
                            call rcvarc(' ', 'TEMP', '-', fami, igau, &
                                        1, tp, iret)
                            call rctype(zi(imate), 1, 'TEMP', [tp], para_vale, &
                                        para_type)
                            if (para_type .eq. 'TEMP') then
                                call utmess('F', 'COMPOR5_5', sk=para_type)
                            end if
                            call rctrac(zi(imate), 1, 'SIGM', para_vale, jprolp, &
                                        jvalep, nbvalp, e)
                            call rcfonc('S', 1, jprolp, jvalep, nbvalp, &
                                        sigy=rp0)
                            call rcfonc('V', 1, jprolp, jvalep, nbvalp, &
                                        p=pm, rp=rp)
                        else
                            call utmess('F', 'ELEMENTS_32')
                        end if
! --- CALCUL DE X
                        cst1 = (rp-rp0)/rp
                        do isig = 1, nbsig
                            xrapel(idecal+isig) = sigt1(isig+(igau-1)* &
                                                        nbsig)*cst1
                        end do
!
                        do isig = 1, nbsig
                            sigmx(isig) = sigt1(isig+(igau-1)*nbsig)-xrapel(idecal+isig)
                        end do
!
                        trace = (sigmx(1)+sigmx(2)+sigmx(3))/3.d0
                        sigmx(1) = sigmx(1)-trace
                        sigmx(2) = sigmx(2)-trace
                        sigmx(3) = sigmx(3)-trace
                        sigeqn = sqrt(1.5d0)*norsig(sigmx, nbsig)
!
                        if (abs(sigeqn-rp0) .le. zernor) then
                            dchax(igau) = -1.d0
                        else
                            dchax(igau) = -2.d0
                            dchay(igau) = norsig(sigmx, nbsig)/rp0
                        end if
                    else if (abs(dchaxm+1.d0) .le. zernor) then
                        do isig = 1, nbsig
                            sigmx(isig) = sigt1(isig)-zr(idera1-1+(igau-1)*nbcmp+2+isig)
                        end do
                        trace = (sigmx(1)+sigmx(2)+sigmx(3))/3.d0
                        sigmx(1) = sigmx(1)-trace
                        sigmx(2) = sigmx(2)-trace
                        sigmx(3) = sigmx(3)-trace
!
                        sigeqn = sqrt(1.5d0)*norsig(sigmx, nbsig)
                        if (compor .eq. 'VMIS_ISOT_LINE') then
                            nomres(1) = 'SY'
                            call rcvalb(fami, igau, 1, '+', zi(imate), &
                                        ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                                        1, nomres, valres, icodre, 2)
                            rp0 = valres(1)
                        else if (compor(10:14) .eq. '_PUIS') then
                            nomres(1) = 'SY'
                            call rcvalb(fami, igau, 1, '+', zi(imate), &
                                        ' ', 'ECRO_PUIS', 0, ' ', [0.d0], &
                                        1, nomres, valres, icodre, 2)
                            rp0 = valres(1)
                        else if (compor(10:14) .eq. '_TRAC') then
                            call rcvarc(' ', 'TEMP', '-', fami, igau, &
                                        1, tp, iret)
                            call rctype(zi(imate), 1, 'TEMP', [tp], para_vale, &
                                        para_type)
                            if (para_type .eq. 'TEMP') then
                                call utmess('F', 'COMPOR5_5', sk=para_type)
                            end if
                            call rctrac(zi(imate), 1, 'SIGM', para_vale, jprolp, &
                                        jvalep, nbvalp, e)
                            call rcfonc('S', 1, jprolp, jvalep, nbvalp, &
                                        sigy=rp0)
                            call rcfonc('V', 1, jprolp, jvalep, nbvalp, &
                                        p=pm, rp=rp)
                        else
                            call utmess('F', 'ELEMENTS_32')
                        end if
                        if (sigeqn .le. rp0) then
                            dchax(igau) = -1.d0
                        else
                            dchax(igau) = -2.d0
                            dchay(igau) = norsig(sigmx, nbsig)/rp0
                            do isig = 1, nbsig
                                xrapel((igau-1)*6+isig) = zr(idera1+( &
                                                             igau-1)*nbcmp+3+isig)
                            end do
                        end if
                    end if
                end if
            else
!  IND_DCHA = -2 AU PAS PRECEDENT T-DT
!        ON STOCKE DANS LE PAS COURANT T, LES VALEURS OBTENUES
!        AU PAS PRECEDENT (T-DT)
                dchax(igau) = -2.d0
                dchay(igau) = zr(idera1+(igau-1)*nbcmp-1+4)
                do isig = 1, nbsig
                    xrapel((igau-1)*6+isig) = zr(idera1+(igau-1)*nbcmp+ &
                                                 3+isig)
                end do
            end if
        end do
    end if
!
! ---- RECUPERATION ET AFFECTATION DU VECTEUR EN SORTIE
! ---- AVEC LE VECTEUR DES INDICATEURS LOCAUX :
!      --------------------------------------
    call jevech('PDERAPG', 'E', idera2)
    do igau = 1, npg
        zr(idera2+(igau-1)*nbcmp-1+1) = dchav(igau)
        zr(idera2+(igau-1)*nbcmp-1+2) = dchat(igau)
        zr(idera2+(igau-1)*nbcmp-1+3) = dchax(igau)
        zr(idera2+(igau-1)*nbcmp-1+4) = dchay(igau)
        do isig = 1, nbsig
            zr(idera2+(igau-1)*nbcmp+3+isig) = xrapel((igau-1)*6+isig)
        end do
        zr(idera2+(igau-1)*nbcmp-1+11) = radiv(igau)
        zr(idera2+(igau-1)*nbcmp-1+12) = radit(igau)
    end do
end subroutine
