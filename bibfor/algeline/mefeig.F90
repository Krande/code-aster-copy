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
subroutine mefeig(ndim, nbmod, matm, matr, mata, &
                  fre, ksi, mavr, alfr, alfi, &
                  mat1, mavi, w, z, ind)
    implicit none
!
#include "asterc/r8pi.h"
#include "asterfort/iunifi.h"
#include "asterfort/mtcrog.h"
#include "asterfort/utmess.h"
#include "asterfort/vphqrp.h"
    integer(kind=8) :: nbmod, ndim(14), ind(2*nbmod)
    real(kind=8) :: matm(nbmod, nbmod), matr(nbmod, nbmod)
    real(kind=8) :: mat1(2*nbmod, 2*nbmod), mavi(2*nbmod, 2*nbmod)
    real(kind=8) :: mata(nbmod, nbmod), fre(nbmod), ksi(nbmod)
    real(kind=8) :: mavr(2*nbmod, 2*nbmod), alfr(2*nbmod), alfi(2*nbmod)
    real(kind=8) :: w(4*nbmod), z(4*nbmod, 2*nbmod)
!     RESOLUTION DU PROBLEME SOUS ECOULEMENT - CALCUL DES VALEURS ET
!     VECTEURS PROPRES DY SYSTEME GENERALISE
!     OPERATEUR APPELANT : OP0144 , FLUST3, MEFIST
! ----------------------------------------------------------------------
!     OPTION DE CALCUL   : CALC_FLUI_STRU , CALCUL DES PARAMETRES DE
!     COUPLAGE FLUIDE-STRUCTURE POUR UNE CONFIGURATION DE TYPE "FAISCEAU
!     DE TUBES SOUS ECOULEMENT AXIAL"
! ----------------------------------------------------------------------
! IN  : NDIM   : TABLEAU DES DIMENSIONS
! IN  : NBMOD  : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
! IN  : MATM   : MATRICE DE MASSE AJOUTEE REPRESENTANT LA PROJECTION DES
!                EFFORTS FLUIDES INERTIELS DANS LA BASE DES DEFORMEES
!                MODALES DES CYLINDRES
! IN  : MATR   : MATRICE DE RAIDEUR  AJOUTEE REPRESENTANT LA PROJECTION
!                DES EFFORTS FLUIDES DE RAIDEUR DANS LA BASE DES
!                DEFORMEES MODALES DES CYLINDRES
! IN  : MATA   : MATRICE D AMORTISSEMENT AJOUTEE REPRESENTANT LA
!                PROJECTION DES EFFORTS FLUIDES D AMORTISSEMENT DANS LA
!                BASE DES DEFORMEES MODALES DES CYLINDRES
! OUT : FRE    : FREQUENCES SOUS ECOULEMENT
! OUT : KSI    : AMORTISSEMENTS SOUS ECOULEMENT
! OUT : MAVR   : TABLEAU DE TRAVAIL ET DEFORMEES MODALES SOUS ECOULEMENT
!                EN BASE MODALE - PARTIE REEL DES VECTEURS PROPRES
! --  : ALFR   : TABLEAU DE TRAVAIL: PARTIE REELE DE LA FREQUENCE
!                COMPLEXE
! --  : ALFI   : TABLEAU DE TRAVAIL: PARTIE IMAGINAIRE DE LA FREQUENCE
!                COMPLEXE
! --  : MAT1   : TABLEAU DE TRAVAIL: MATRICE A DU SYSTEME A.X = L.X
! --  : MAVI   : TABLEAU DE TRAVAIL POUR L INVERSION DE LA MATRICE DE
!                MASSE - PARTIE IMAGINAIRE DES VECTEURS PROPRES
! --  : W      : TABLEAU DE TRAVAIL: VALEURS PROPRES COMPLEXES APRES
!                RESOLUTION DU SYSTEME
! --  : Z      : TABLEAU DE TRAVAIL: VECTEURS PROPRES COMPLEXES APRES
!                RESOLUTION DU SYSTEME
! OUT : IND    : TABLEAU DE TRAVAIL: INDICES DES VALEURS PROPRES
!                RETENUES, CORRESPONDANT A UNE FREQUENCE POSITIVE, ET
!                ORDONNEES DE FACON CROISSANTE
! ----------------------------------------------------------------------
    integer(kind=8) :: i, j, k
    character(len=3) :: note
    real(kind=8) :: norema, noremi, noimma, noimmi, normr, normi, norm
    real(kind=8) :: temp
! ----------------------------------------------------------------------
!
! --- LECTURE DES DIMENSIONS
!-----------------------------------------------------------------------
    integer(kind=8) :: icode, ier, ifm, ihh, jhh, n
    integer(kind=8) :: nitqr
    real(kind=8) :: a, adiv, alpha, alphai, alphar, b, bdiv
    real(kind=8) :: beta, betai, betar, c, det, pi
    real(kind=8) :: snor, u, v, vnorma, vnormb
!-----------------------------------------------------------------------
    nbmod = ndim(2)
!
!
    pi = r8pi()
!
!
! --- INVERSION DE LA MATRICE DE MASSE MAT1 = MATM
!     MAVI EST L IDENTITE
!
    do i = 1, nbmod
        do j = 1, nbmod
            mat1(i, j) = matm(i, j)
            mavi(i, j) = 0.d0
        end do
        mavi(i, i) = 1.d0
    end do
!
    ier = 1
    call mtcrog(mat1, mavi, 2*nbmod, nbmod, nbmod, &
                mavr, alfr, ier)
    if (ier .ne. 0) then
        call utmess('F', 'ALGELINE_77')
    end if
!
!
! --- CALCUL DU PRODUIT DE L INVERSE DE LA MATRICE DE MASSE PAR LES
! --- MATRICES DE RAIDEUR ET D AMORTISSEMENT
! --- D AMORTISSEMENT
! --- CONSTRUCTION DE LA MATRICE DU PROBLEME MODALE GENERALISE MAR1 = A
    do i = 1, nbmod
        do j = 1, nbmod
            mat1(i, nbmod+j) = 0.d0
            mat1(nbmod+i, nbmod+j) = 0.d0
            mat1(i, j) = 0.d0
            mat1(nbmod+i, j) = 0.d0
            do k = 1, nbmod
                mat1(nbmod+i, j) = mat1(nbmod+i, j)-mavr(i, k)*matr(k, j)
                mat1(nbmod+i, nbmod+j) = mat1(nbmod+i, nbmod+j)-mavr(i, k)*mata(k, j)
            end do
        end do
        mat1(i, nbmod+i) = 1.d0
    end do
!
!
! --- RESOLUTION DU SYSTEME A.X = L.X
!
    ier = 1
    icode = 1
    call vphqrp(mat1, 2*nbmod, 2*nbmod, icode, w, &
                z, 2*nbmod, mavr, 30, ier, &
                nitqr)
!
    if (ier .eq. 1) then
        call utmess('F', 'ALGELINE_78')
    end if
!
! --- ALFR: PARTIES REELLES DES VALEURS PROPRES
! --- ALFI: PARTIES IMAGINAIRES DES VALEURS PROPRES
! --- MAVR: PARTIES REELLES DES VECTEURS PROPRES
! --- MAVI: PARTIES IMAGINAIRES DES VECTEURS PROPRES
    do ihh = 1, 2*nbmod
        alfr(ihh) = w(2*(ihh-1)+1)
        alfi(ihh) = w(2*(ihh-1)+2)
        do jhh = 1, 2*nbmod
            mavr(jhh, ihh) = z(2*(jhh-1)+1, ihh)
            mavi(jhh, ihh) = z(2*(jhh-1)+2, ihh)
!         MAVR(JHH,IHH) = Z(2*(JHH+(IHH-1)*NBMOD-1)+1)
!         MAVI(JHH,IHH) = Z(2*(JHH+(IHH-1)*NBMOD-1)+2)
        end do
    end do
!
!
! --- MINIMISATION DE LA PARTIE IMAGINAIRE DES VECTEURS PROPRES
! --- PAR MULTIPLICATION PAR UNE CONSTANTE COMPLEXE
! --- LES MISES A JOUR SONT FAITES SUR LA TOTALITE DES VECTEURS
! --- LES CALCULS DES NORMES SONT FAITS SUR LES NBMOD PREMIERS TERMES
!
! --  BOUCLE SUR LES VECTEURS PROPRES
    do i = 1, 2*nbmod
!
! --  ON EFFECTUE LA MINIMISATION QUE SUR LES VECTEURS PROPRES POUR
! --  LESQUELS LES PARTIES IMAGINAIRES DES VALEURS PROPRES
! --  CORRESPONDANTES SONT POSITIVES
        if (alfi(i) .ge. 0.d0) then
!
! --        CALCUL DE LA SOMME DES CARRES DES PARTIES REELLES, DES
!           PARTIES IMAGINAIRES ET DES DOUBLES PRODUITS CROISES
            a = 0.d0
            b = 0.d0
            c = 0.d0
            do j = 1, nbmod
                a = a+mavr(j, i)*mavr(j, i)
                b = b+mavi(j, i)*mavi(j, i)
                c = c+mavr(j, i)*mavi(j, i)*2
            end do
!
            if (a .ne. 0.d0 .and. b .ne. 0.d0 .and. c .ne. 0.d0) then
!
! --        CALCUL DES COMPLEXES MULTIPLICATEURS ALPHA ET BETA
!           ALPHA = ALPHAR + I . ALPHAI
!           BETA  = BETAR  + I . BETAI
                u = a*a-b*b
                v = a*c+b*c
                det = u*u+v*v
                alpha = (u+sqrt(det))/v
                adiv = sqrt(alpha*alpha+1.d0)
                alphar = 1.d0/adiv
                alphai = alpha/adiv
                beta = (u-sqrt(det))/v
                bdiv = sqrt(beta*beta+1.d0)
                betar = 1.d0/bdiv
                betai = beta/bdiv
!
! --        CALCUL DES PRODUITS DES VECTEURS PROPRES PAR ALPHA ET BETA
                do j = 1, 2*nbmod
                    z(j, 1) = (mavr(j, i)*alphar-mavi(j, i)*alphai)
                    z(j+2*nbmod, 1) = (mavr(j, i)*alphai+mavi(j, i)*alphar)
                    z(j, 2) = (mavr(j, i)*betar-mavi(j, i)*betai)
                    z(j+2*nbmod, 2) = (mavr(j, i)*betai+mavi(j, i)*betar)
                end do
!
! --        CALCUL DE LA SOMME DES CARRES DES PARTIES IMAGINAIRES
                vnorma = 0.d0
                vnormb = 0.d0
                do j = 1, nbmod
                    vnorma = vnorma+z(j+2*nbmod, 1)*z(j+2*nbmod, 1)
                    vnormb = vnormb+z(j+2*nbmod, 2)*z(j+2*nbmod, 2)
                end do
!
! --        MISE A JOUR DU VECTEUR PROPRE OBTENUS AVEC LE COEFFICIENT
!           ALPHA OU BETA MINIMISANT LA NORME DE LA PARTIE IMAGINAIRE
                if (vnorma .lt. vnormb) then
                    do j = 1, 2*nbmod
                        mavr(j, i) = z(j, 1)
                        mavi(j, i) = z(j+2*nbmod, 1)
                    end do
                else
                    do j = 1, 2*nbmod
                        mavr(j, i) = z(j, 2)
                        mavi(j, i) = z(j+2*nbmod, 2)
                    end do
                end if
!
! --        TRAITEMENT DU CAS OU A OU B OU C = 0
! --        SI B = 0 ON A UN VECTEUR PROPRE REEL: ON NE FAIT RIEN
! --        SI A = 0 ON A UN VECTEUR PROPRE IMAGINAIRE PUR, ON MULTIPLIE
! --                 CE DERNIER PAR -I
! --        SI C=0 ET B>A, ON MULTIPLIE LE VECTEUR PROPRE PAR -I
            else if (a .eq. 0.d0) then
                do j = 1, 2*nbmod
                    mavr(j, i) = mavi(j, i)
                    mavi(j, i) = 0.d0
                end do
            else if (b .gt. a .and. c .eq. 0.d0) then
                do j = 1, 2*nbmod
                    u = mavi(j, i)
                    mavi(j, i) = -mavr(j, i)
                    mavr(j, i) = u
                end do
            end if
        end if
!
! --  FIN DE BOUCLE SUR LES VECTEURS PROPRES
    end do
!
! --- NORMALISATION DES VECTEURS PROPRES COMPLEXES (SUR LES NBMOD
! --- PREMIERES COMPOSANTES )
!
    do i = 1, 2*nbmod
        snor = 0.d0
        do j = 1, nbmod
            snor = snor+mavr(j, i)**2+mavi(j, i)**2
        end do
        snor = snor**0.5d0
        if (snor .ne. 0.d0) then
            do j = 1, 2*nbmod
                mavr(j, i) = mavr(j, i)/snor
                mavi(j, i) = mavi(j, i)/snor
            end do
        end if
    end do
!
! --- VALEUR PROPRE (J) = (ALFR(J)+I*ALFI(J))
! --- FREQUENCE COMPLEXE : S = ALFR(J) + I*ALFI(J)
! --- PARTIE IMAGINAIRE ALFI(J) = 2*PI*FREQUENCE
! ---                           * RACINE(1-AMORTISSEMENT**2)
! --- PARTIE REELLE     ALFR(J) = -PULSATION*AMORTISSEMENT
!
!
! --- SUPPRESSION DES MODES MATHEMATIQUES DE FREQUENCE NEGATIVE
! --- RECHERCHE DES INDICES DES VALEURS PROPRES CORRESPONDANT A UNE
! --- FREQUENCE POSITIVE
!
    n = 0
    do i = 1, 2*nbmod
        if (alfi(i) .ge. 0.d0) then
            n = n+1
            ind(n) = i
        end if
    end do
!
    if (n .ne. nbmod) then
        write (note(1:3), '(I3.3)') n
        call utmess('F', 'ALGELINE_79', sk=note)
    end if
!
! --- FREQUENCE ET AMORTISSEMENT REELS
!
    do i = 1, nbmod
        fre(i) = sqrt(alfr(ind(i))*alfr(ind(i))+alfi(ind(i))*alfi(ind(i)))/2.d0/pi
        ksi(i) = -alfr(ind(i))/2.d0/pi/fre(i)
    end do
!
! ---- CLASSEMENT PAR ORDRE CROISSANT DE FREQUENCES
!
    do i = 1, nbmod-1
        do j = i+1, nbmod
            if (fre(j) .lt. fre(i)) then
                temp = fre(i)
                fre(i) = fre(j)
                fre(j) = temp
                temp = ksi(i)
                ksi(i) = ksi(j)
                ksi(j) = temp
                k = ind(i)
                ind(i) = ind(j)
                ind(j) = k
            end if
        end do
    end do
!
! --- CALCUL DES NORMES DES PARTIES REELLES ET IMAGINAIRES DES
! --- VECTEURS PROPRES ET IMPRESSION SUR LE FICHIER RESULTAT
! --- RECHERCHE DES MAXIMA ET MINIMA
!
    ifm = iunifi('MESSAGE')
    write (ifm, *)
    write (ifm, *) '==============================================='
    write (ifm, *)
    write (ifm, *) '     NORMES RELATIVES DES PARTIES REELLES '
    write (ifm, *)
    write (ifm, *) '     ET IMAGINAIRES  DES VECTEURS PROPRES'
    write (ifm, *)
    write (ifm, *) '==============================================='
    write (ifm, *)
    write (ifm, 6000)
    write (ifm, 6001)
    write (ifm, 6002)
    write (ifm, 6003)
    write (ifm, 6001)
    write (ifm, 6000)
    write (ifm, 6000)
!
    noremi = 1.0d0
    norema = 0.0d0
    noimmi = 1.0d0
    noimma = 0.0d0
!
    do i = 1, nbmod
        normr = 0.d0
        normi = 0.d0
        do j = 1, nbmod
            normr = normr+mavr(j, ind(i))*mavr(j, ind(i))
            normi = normi+mavi(j, ind(i))*mavi(j, ind(i))
        end do
        norm = normi+normr
        norm = sqrt(norm)
        normr = sqrt(normr)/norm
        normi = sqrt(normi)/norm
        if (noremi .gt. normr) noremi = normr
        if (norema .lt. normr) norema = normr
        if (noimmi .gt. normi) noimmi = normi
        if (noimma .lt. normi) noimma = normi
    end do
    write (ifm, 6005) norema, noimmi
    write (ifm, 6001)
    write (ifm, 6006) noremi, noimma
    write (ifm, 6001)
    write (ifm, 6000)
    write (ifm, *)
!
!
! --- FORMATS
!
6000 format(1x, '***************************************************',&
&       '*****************')
6001 format(1x, '*              *                         *          ',&
&       '               *')
6002 format(1x, '*    NUMERO    *   NORME RELATIVE DE LA  *   NORME '&
&       , 'RELATIVE DE LA  *')
6003 format(1x, '*  DE VECTEUR  *     PARTIE REELLE       *    PARTIE'&
&       , ' IMAGINAIRE    *')
6005 format(1p, 1x, '*   MAX :    ', '  * ', 5x, d13.6, 4x, '  *  '&
&       , 4x, d13.6, 5x, ' *')
6006 format(1p, 1x, '*   MIN :    ', '  * ', 5x, d13.6, 4x, '  *  '&
&       , 4x, d13.6, 5x, ' *')
!
!
!
end subroutine
