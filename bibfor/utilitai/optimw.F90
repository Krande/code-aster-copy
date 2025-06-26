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
subroutine optimw(method, nrupt, x, y, prob, &
                  sigw, nt, nur, nbres, calm, &
                  cals, mk, sk, mkp, skp, &
                  impr, ifm, dept, indtp, nbtp)
    implicit none
#include "asterf_types.h"
#include "asterfort/ntweib.h"
#include "asterfort/utmess.h"
    integer(kind=8) :: nrupt, nt(*), nbres, nur(*), ir, indtp(*), nbtp, ifm
    real(kind=8) :: x(*), y(*), sigw(*), mk, sk(*), mkp, skp(*), prob(*)
    character(len=16) :: method
    aster_logical :: calm, cals, impr, dept
!     AUTEUR : M. BONNAMY
!     ----------------------------------------------------------------
!
!     BUT: CALCUL DE RECALAGE DES PARAMETRES DE WEIBULL
!
!     ----------------------------------------------------------------
!
!     METHOD       /IN/:METHODE DE CALAGE
!     NRUPT        /IN/:NOMBRE TOTAL DE CONTRAINTES DE WEIBULL
!     SIGW         /IN/:CONTRAINTES DE WEIBULL AUX INSTANTS DE RUPTURE
!     NT           /IN/:DIMENSION DE LA SOUS-BASE CORRESPONDANT A LA
!                       TEMPERATURE T
!     NUR          /IN/:NUMERO DE RESULTAT ASSOCIEE A
!                       LA CONTRAINTE SIGW(I)
!     NBRES        /IN/:NOMBRE DE BASES DE RESULTATS
!     MK           /IN/:PARAMETRE M(K)DE WEIBULL
!     SK           /IN/:PARAMETRE SIGMA-U(K) DE WEIBULL
!     CALM         /IN/:TRUE SI M EST FIXE
!     CALS         /IN/:TRUE SI SIGMA_U EST FIXE
!     IMPR         /IN/:IMPRESSION DETAILLEE
!     DEPT         /IN/:DEPENDANCE EN TEMPERATURE POUR SIGMA-U
!     INDTP        /IN/:INDICE DE TEMPERATURE POUR CHAQUE RESULTAT
!     NBTP         /IN/:NOMBRE DE TEMPERATURE DIFFERENTES
!
!     X,Y          /OUT/:VALEUR DES FONCTIONS LOG(SIGMAW)
!                        ET LOG(LOG(1/(1-PF)))
!     PROB         /OUT/:PROBABILITE THEORIQUE POUR REGRESSION LINEAIRE
!     MKP          /OUT/:PARAMETRE M(K+1)DE WEIBULL
!     SKP          /OUT/:PARAMETRE SIGMA-U(K+1) DE WEIBULL
!
!     ----------------------------------------------------------------
!
    real(kind=8) :: syi, sxi, sxixi, sxiyi, sxiyj, sxixj, unsurn, unsurm
    real(kind=8) :: swm, prec, mg, md, prov, snt, s1, s2
    integer(kind=8) :: i, j, k, itp, irg
!     ----------------------------------------------------------------
!
    if (method(1:9) .eq. 'REGR_LINE') then
!
!        METHODE DE REGRESSION LINEAIRE
!
!
        if (.not. dept) then
!
!       UN SEUL RESU : PAS DE DEPENDANCE EN TEMPERATURE
!
            syi = 0.d0
            sxi = 0.d0
            sxixi = 0.d0
            sxiyi = 0.d0
!
            do i = 1, nrupt
!
                prob(i) = i
                s1 = nrupt
                prob(i) = prob(i)/(s1+1.d0)
                y(i) = log(log(1.d0/(1.d0-prob(i))))
                x(i) = log(sigw(i))
                syi = syi+y(i)
                sxi = sxi+x(i)
                sxixi = sxixi+x(i)*x(i)
                sxiyi = sxiyi+x(i)*y(i)
!
            end do
!
            sxiyj = 0.d0
            sxixj = 0.d0
!
            do i = 1, nrupt
!
                do j = 1, nrupt
!
                    sxiyj = sxiyj+x(i)*y(j)
                    sxixj = sxixj+x(i)*x(j)
!
                end do
!
            end do
!
            unsurn = nrupt
            unsurn = 1.d0/unsurn
!
            if ((.not. calm) .and. (.not. cals)) then
                mkp = (unsurn*sxiyj-sxiyi)/(unsurn*sxixj-sxixi)
                skp(1) = exp(unsurn*(sxi-(1.d0/mkp)*syi))
            else if (calm) then
                mkp = mk
                skp(1) = exp(unsurn*(sxi-(1.d0/mkp)*syi))
            else if (cals) then
                skp(1) = sk(1)
                mkp = sxiyi/(sxixi-log(skp(1))*sxi)
            end if
            if (impr) write (ifm, *) 'M(K) =', mkp, 'SIGU(K) = ', skp(1)
!
            sxi = 0.d0
            do j = 1, nrupt
!
                prov = (1.d0-exp(-(sigw(j)/sk(1))**mk))
                if (prov .ne. 1.d0) prov = log(log(1.d0/(1.d0-prov)))
                sxi = sxi+(y(j)-prov)**2
!
            end do
            if (impr) write (ifm, *) 'ECART THEORIE-EXPERIENCE AU DEBUT DE L''ITERATION : ', sxi
!
        else
!
!       DEPENDANCE EN TEMPERATURE SIGMA_U(T)
!
            sxixi = 0.d0
            sxiyi = 0.d0
!
            do itp = 1, nbtp
!
                snt = 0.d0
!
                do ir = 1, nbres
!
                    if (indtp(ir) .eq. itp) snt = snt+nt(ir)
!
                end do
!
                do i = 1, nrupt
!
                    irg = 1
                    do k = 1, i-1
                        if (indtp(nur(k)) .eq. itp) then
                            irg = irg+1
                        end if
                    end do
!
                    if (indtp(nur(i)) .eq. itp) then
!
                        prob(i) = irg
                        prob(i) = prob(i)/(snt+1.d0)
                        y(i) = log(log(1.d0/(1.d0-prob(i))))
                        x(i) = log(sigw(i))
                        sxixi = sxixi+x(i)*x(i)
                        sxiyi = sxiyi+x(i)*y(i)
!
                    end if
!
                end do
!
            end do
!
            s1 = 0.d0
            s2 = 0.d0
!
            do itp = 1, nbtp
!
                sxiyj = 0.d0
                sxixj = 0.d0
                snt = 0.d0
!
                do ir = 1, nbres
!
                    if (indtp(ir) .eq. itp) snt = snt+nt(ir)
!
                end do
!
                do i = 1, nrupt
!
                    do j = 1, nrupt
!
                        if (indtp(nur(i)) .eq. itp .and. indtp(nur(j)) .eq. itp) then
                            sxiyj = sxiyj+x(i)*y(j)
                            sxixj = sxixj+x(i)*x(j)
                        end if
!
                    end do
!
                end do
                s1 = s1+sxiyj/snt
                s2 = s2+sxixj/snt
!
            end do
!
            if ((.not. calm)) then
                mkp = (s1-sxiyi)/(s2-sxixi)
            else if (calm) then
                mkp = mk
            end if
!
            if (impr) write (ifm, *) 'M(K) =', mkp
!
            if (((.not. calm) .and. (.not. cals)) .or. calm) then
!
!          (M ET SIGMA-U) OU (SIGMA-U) SONT A RECALER
!
                do itp = 1, nbtp
!
                    snt = 0.d0
!
                    do ir = 1, nbres
!
                        if (indtp(ir) .eq. itp) snt = snt+nt(ir)
!
                    end do
!
                    sxi = 0.d0
                    syi = 0.d0
!
                    do i = 1, nrupt
!
                        if (indtp(nur(i)) .eq. itp) then
                            syi = syi+y(i)
                            sxi = sxi+x(i)
                        end if
!
                    end do
!
                    skp(itp) = exp((sxi-(1.d0/mkp)*syi)/snt)
                    if (impr) write (ifm, *) 'S(K) (', itp, ')=', skp(itp)
!
                end do
!
            else if (cals) then
!
                do ir = 1, nbtp
!
                    skp(ir) = sk(ir)
                    if (impr) write (ifm, *) 'S(K) (', ir, ')=', skp(ir)
!
                end do
!
            end if
!
        end if
!
    else if (method(1:9) .eq. 'MAXI_VRAI') then
!
!        METHODE DU MAXIMUM DE VRAISSEMBLANCE
!
        if ((.not. calm) .and. (.not. cals)) then
!
!        M ET SIGMA-U SONT A RECALER
!
            prec = 1.d-8
            mg = 1.d0
            md = mk
            swm = 0.d0
            unsurn = nrupt
            unsurn = 1.d0/unsurn
            if (impr) write (ifm, *) 'RESOLUTION F(M)=0 PAR NEWTON'
!
!          RESOLUTION DE L'EQUATION NON LINEAIRE F(M)=0
!
            call ntweib(nrupt, cals, sk, sigw, nur, &
                        nt, nbres, mg, md, prec, &
                        mkp, impr, ifm, indtp, nbtp)
!
            unsurm = 1.d0/mkp
!
!          CALCUL DU SIGMA-U
!
            do itp = 1, nbtp
!
                snt = 0.d0
                do ir = 1, nbres
!
                    if (indtp(ir) .eq. itp) snt = snt+nt(ir)
!
                end do
!
                swm = 0.d0
                do i = 1, nrupt
!
                    if (indtp(nur(i)) .eq. itp) then
                        swm = swm+sigw(i)**mkp
                    end if
!
                end do
!
                skp(itp) = (swm/snt)**(unsurm)
!
            end do
!
        else if (calm) then
!
!        M EST CALE
!
            mkp = mk
            unsurm = 1.d0/mkp
!
            do itp = 1, nbtp
!
                snt = 0.d0
                do ir = 1, nbres
                    if (indtp(ir) .eq. itp) snt = snt+nt(ir)
                end do
!
                swm = 0.d0
                do i = 1, nrupt
!
                    if (indtp(nur(i)) .eq. itp) then
                        swm = swm+sigw(i)**mkp
                    end if
!
                end do
!
                skp(itp) = (swm/snt)**(unsurm)
!
            end do
!
        else if (cals) then
!
!        SIGMA-U EST CALE
!
            do ir = 1, nbtp
                skp(ir) = sk(ir)
            end do
            prec = 1.d-8
            mg = 1.d0
            md = mk
!
!          RESOLUTION DE L'EQUATION NON LINEAIRE F(M)=0
!
            if (impr) write (ifm, *) 'RESOLUTION F(M)=0 PAR NEWTON'
            call ntweib(nrupt, cals, sk, sigw, nur, &
                        nt, nbres, mg, md, prec, &
                        mkp, impr, ifm, indtp, nbtp)
!
        end if
!
        if (impr) then
            write (ifm, *) 'M(K) =', mkp
            do ir = 1, nbtp
                write (ifm, *) 'S(K) (', ir, ')=', skp(ir)
            end do
        end if
!
    end if
!
    if (mkp .lt. 1.d0) then
        call utmess('F', 'UTILITAI3_36')
    end if
!
end subroutine
