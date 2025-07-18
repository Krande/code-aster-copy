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
subroutine pbflkc(umoy, rhof, hmoy, rmoy, long, &
                  cf0, mcf0, icoq, imod, nbm, &
                  rkip, tcoef, s1, s2, ki, &
                  lambda, kcalcu, passag)
    implicit none
! COUPLAGE FLUIDELASTIQUE, CONFIGURATIONS DU TYPE "COQUE_COAX"
! RESOLUTION DU PROBLEME FLUIDE INSTATIONNAIRE : INITIALISATION DE
! KCALCU(3,4) DANS LE CAS OU UMOY <> 0
! APPELANT : PBFLUI
!-----------------------------------------------------------------------
!  IN : UMOY   : VITESSE DE L'ECOULEMENT MOYEN
!  IN : RHOF   : MASSE VOLUMIQUE DU FLUIDE
!  IN : HMOY   : JEU ANNULAIRE MOYEN
!  IN : RMOY   : RAYON MOYEN
!  IN : LONG   : LONGUEUR DU DOMAINE DE RECOUVREMENT DES DEUX COQUES
!  IN : CF0    : COEFFICIENT DE FROTTEMENT VISQUEUX
!  IN : MCF0   : EXPOSANT VIS-A-VIS DU NOMBRE DE REYNOLDS
!  IN : ICOQ   : INDICE CARACTERISANT LA COQUE SUR LAQUELLE ON TRAVAILLE
!                ICOQ=1 COQUE INTERNE  ICOQ=2 COQUE EXTERNE
!  IN : IMOD   : INDICE DU MODE CONSIDERE
!  IN : NBM    : NOMBRE DE MODES PRIS EN COMPTE POUR LE COUPLAGE
!  IN : RKIP   : ORDRE DE COQUE DU MODE CONSIDERE, PONDERE PAR LA VALEUR
!                MOYENNE DU PROFIL DE PRESSION
!  IN : TCOEF  : TABLEAU DES COEFFICIENTS DES DEFORMEES AXIALES
!  IN : S1     : PARTIE REELLE     DE LA FREQUENCE COMPLEXE
!  IN : S2     : PARTIE IMAGINAIRE DE LA FREQUENCE COMPLEXE
!  IN : KI     : TABLEAU DE TRAVAIL
!  IN : LAMBDA : VALEURS PROPRES DE L'OPERATEUR DIFFERENTIEL
! OUT : KCALCU : MATRICE RECTANGULAIRE A COEFFICIENTS CONSTANTS
!                PERMETTANT DE CALCULER UNE SOLUTION PARTICULIERE DU
!                PROBLEME FLUIDE INSTATIONNAIRE, LORSQUE UMOY <> 0
! OUT : PASSAG : MATRICE DONT LES COLONNES SONT LES VECTEURS PROPRES DE
!                L'OPERATEUR DIFFERENTIEL
!-----------------------------------------------------------------------
!
    real(kind=8) :: umoy, rhof, hmoy, rmoy, long, cf0, mcf0
    integer(kind=8) :: icoq, imod, nbm
    real(kind=8) :: rkip, tcoef(10, nbm), s1, s2
    complex(kind=8) :: ki(4, 3), lambda(3), kcalcu(3, 4), passag(3, 3)
!
    real(kind=8) :: ln
    complex(kind=8) :: c1, c2, c3, c4, d1, d2, d3, d4, e1, e2, e3, e4, f1, f2
    complex(kind=8) :: f3, f4
    complex(kind=8) :: g1, g2, g3, g4, q11, q12, q13, s, p1, p2, p3, p4, p5, t
    complex(kind=8) :: u, v, j
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: itab, m1, m2
    real(kind=8) :: a1, a2, a3, a4, b1, b2, b3
    real(kind=8) :: b4, poids, w
!-----------------------------------------------------------------------
    itab = 0
    poids = -1.d0
    if (icoq .eq. 2) then
        itab = 5
        poids = 1.d0
    end if
    ln = tcoef(1+itab, imod)
    a1 = tcoef(2+itab, imod)*poids
    a2 = tcoef(3+itab, imod)*poids
    a3 = tcoef(4+itab, imod)*poids
    a4 = tcoef(5+itab, imod)*poids
    b1 = tcoef(2+itab, imod)/2.d0
    b2 = tcoef(3+itab, imod)/2.d0
    b3 = tcoef(4+itab, imod)/2.d0
    b4 = tcoef(5+itab, imod)/2.d0
!
    s = dcmplx(s1, s2)
    j = dcmplx(0.d0, 1.d0)
!
    p1 = dcmplx(ln/(long*hmoy))
    p2 = s/(umoy*hmoy)
    p3 = dcmplx(ln/(long*rmoy))
    p4 = s/(umoy*rmoy)
    p5 = (dcmplx(cf0/hmoy)+s/umoy)/hmoy
!
    c1 = p1*a2+p2*a1+p3*b2+p4*b1
    c2 = -1.d0*p1*a1+p2*a2-p3*b1+p4*b2
    c3 = p1*a4+p2*a3+p3*b4+p4*b3
    c4 = p1*a3+p2*a4+p3*b3+p4*b4
!
    d1 = p1*a2+p5*a1+p3*b2+p4*b1
    d2 = -1.d0*p1*a1+p5*a2-p3*b1+p4*b2
    d3 = p1*a4+p5*a3+p3*b4+p4*b3
    d4 = p1*a3+p5*a4+p3*b3+p4*b4
!
    u = 0.5d0*((rkip/rmoy)**2)*(s+dcmplx((cf0/hmoy)*(mcf0+2.d0)*umoy))
    v = -0.5d0*((rkip/rmoy)**2)*dcmplx(umoy)
!
    e1 = (u*(j*c2-c1))+v*lambda(1)*(d1-j*d2)
    e2 = (u*(-1.d0*j*c2-c1))+v*lambda(1)*(d1+j*d2)
    e3 = (-1.d0*u*(c3+c4))+v*lambda(1)*(d3+d4)
    e4 = (u*(c4-c3))+v*lambda(1)*(d3-d4)
!
    f1 = (u*(j*c2-c1))+v*lambda(2)*(d1-j*d2)
    f2 = (u*(-1.d0*j*c2-c1))+v*lambda(2)*(d1+j*d2)
    f3 = (-1.d0*u*(c3+c4))+v*lambda(2)*(d3+d4)
    f4 = (u*(c4-c3))+v*lambda(2)*(d3-d4)
!
    g1 = (u*(j*c2-c1))+v*lambda(3)*(d1-j*d2)
    g2 = (u*(-1.d0*j*c2-c1))+v*lambda(3)*(d1+j*d2)
    g3 = (-1.d0*u*(c3+c4))+v*lambda(3)*(d3+d4)
    g4 = (u*(c4-c3))+v*lambda(3)*(d3-d4)
!
    t = -1.d0*(s/umoy+dcmplx(cf0/hmoy))
    u = 3.d0*((rkip/rmoy)**2)*(s/umoy+dcmplx((cf0/hmoy)*(mcf0+2.d0)))
    v = dcmplx(2.d0*((rkip/rmoy)**2))
!
    q11 = (u+lambda(1)*(t*lambda(1)+v))**(-1.d0)
    q12 = (u+lambda(2)*(t*lambda(2)+v))**(-1.d0)
    q13 = (u+lambda(3)*(t*lambda(3)+v))**(-1.d0)
!
    w = ln/long
!
    ki(1, 1) = q11*e1/(j*w-lambda(1))
    ki(1, 2) = q12*f1/(j*w-lambda(2))
    ki(1, 3) = q13*g1/(j*w-lambda(3))
    ki(2, 1) = -1.d0*q11*e2/(j*w+lambda(1))
    ki(2, 2) = -1.d0*q12*f2/(j*w+lambda(2))
    ki(2, 3) = -1.d0*q13*g2/(j*w+lambda(3))
    ki(3, 1) = q11*e3/(dcmplx(w)-lambda(1))
    ki(3, 2) = q12*f3/(dcmplx(w)-lambda(2))
    ki(3, 3) = q13*g3/(dcmplx(w)-lambda(3))
    ki(4, 1) = -1.d0*q11*e4/(dcmplx(w)+lambda(1))
    ki(4, 2) = -1.d0*q12*f4/(dcmplx(w)+lambda(2))
    ki(4, 3) = -1.d0*q13*g4/(dcmplx(w)+lambda(3))
!
    do m2 = 1, 3
        passag(1, m2) = dcmplx(1.d0, 0.d0)
        passag(2, m2) = -1.d0*lambda(m2)*rmoy
        passag(3, m2) = -1.d0*lambda(m2)*(lambda(m2)-t)*(rmoy*rmoy*rhof*umoy)
    end do
!
    do m2 = 1, 4
        do m1 = 1, 3
            kcalcu(m1, m2) = passag(m1, 1)*ki(m2, 1)+passag(m1, 2)*ki(m2, 2)
            kcalcu(m1, m2) = kcalcu(m1, m2)+passag(m1, 3)*ki(m2, 3)
        end do
    end do
!
end subroutine
