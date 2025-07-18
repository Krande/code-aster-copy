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
! aslint: disable=W0413
! => real zero (init by calcul.F90)
!
subroutine forngr(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/btsig.h"
#include "asterfort/jacbm1.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/jm1dn1.h"
#include "asterfort/jm1dn2.h"
#include "asterfort/matbmn.h"
#include "asterfort/matbmr.h"
#include "asterfort/matbsr.h"
#include "asterfort/matbsu.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "asterfort/terefe.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "asterfort/vectpe.h"
#include "asterfort/vectrn.h"
#include "blas/daxpy.h"
    character(len=16) :: option, nomte
!
!     FONCTION  :  FORC_NODA DES COQUE_3D
!                  GEOMETRIQUE AVEC GRANDES ROTATIONS
!
!     ARGUMENTS :
!     DONNEES   :      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
!     OPTIONS   :     FORC_NODA      : FORCE INTERNE
!
! ......................................................................
!
!
!---- DECLARATIONS BIDONS
!
    real(kind=8) :: bid33(3, 3)
!
!---- DECLARATIONS LOCALES
!
    integer(kind=8) :: i, j, in, ii, nval, k1, iret, itab(7)
!
!---- DECLARATIONS RIGIDITE GEOMETRIQUE
!
    real(kind=8) :: stild(5)
!
!---- DECLARATIONS STANDARDS
!
    integer(kind=8) :: igeom, icontm, ivectu
    integer(kind=8) :: lzi, lzr, jcara
    integer(kind=8) :: nb1, nb2
!
!---- DECLARATIONS PROPRES COQUE_3D NON LINEAIRE
!
    integer(kind=8) :: inte, intsr, intsn
    real(kind=8) :: eptot
    integer(kind=8) :: npge, npgsr, npgsn
    parameter(npge=3)
    real(kind=8) :: vecta(9, 2, 3)
    real(kind=8) :: vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vecnph(9, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: jm1(3, 3), detj
    real(kind=8) :: jdn1ri(9, 51), jdn1rc(9, 51)
    real(kind=8) :: jdn1ni(9, 51), jdn1nc(9, 51)
    real(kind=8) :: jdn2rc(9, 51)
    real(kind=8) :: jdn2nc(9, 51)
    real(kind=8) :: ksi3s2
!
!---- DECLARATIONS COUCHES
!
    integer(kind=8) :: nbcou, nbsp
    integer(kind=8) :: icou
    real(kind=8) :: zic, zmin, epais, coef
! --- CONTRAINTE DE REFERENCE POUR REFE_FORC_NODA
    real(kind=8) :: sigref
!
!---- DECLARATIONS COQUE NON LINEAIRE
!
    integer(kind=8) :: jvDisp
    real(kind=8) :: b1su(5, 51), b2su(5, 51)
    real(kind=8) :: b1src(2, 51, 4)
    real(kind=8) :: b2src(2, 51, 4)
    real(kind=8) :: b1mnc(3, 51), b1mni(3, 51)
    real(kind=8) :: b2mnc(3, 51), b2mni(3, 51)
    real(kind=8) :: b1mri(3, 51, 4)
    real(kind=8) :: b2mri(3, 51, 4)
    real(kind=8) :: dudxri(9), dudxni(9)
    real(kind=8) :: dudxrc(9), dudxnc(9)
    real(kind=8) :: vecu(8, 3), vecthe(9, 3)
    real(kind=8) :: vecpe(51)
!    POUR_RESI_REFE_RELA
    real(kind=8) :: sigtmp(5), ftemp(51), effint(51)
    character(len=16) :: kmess(2)
!
!---- DECLARATIONS ROTATION GLOBAL LOCAL AU NOEUDS
!
    integer(kind=8) :: jnbspi
!
    real(kind=8) :: blam(9, 3, 3)
    blas_int :: b_incx, b_incy, b_n
!
!
! DEB
!
!---- LE NOMBRE DE COUCHES
!
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
!
    if (nbcou .le. 0) then
        call utmess('F', 'ELEMENTS_12')
    end if
!
!______________________________________________________________________
!
!---- RECUPERATION DES POINTEURS ( L : LECTURE )
!______________________________________________________________________
!
!....... GEOMETRIE INITIALE ( COORDONNEES INITIALE DES NOEUDS )
!
    call jevech('PGEOMER', 'L', igeom)
!
!---- RECUPERATION DES OBJETS INITIALISES
!
!....... LES ENTIERS
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
!
!------- NOMBRE DE NOEUDS ( NB1 : SERENDIP , NB2 : LAGRANGE )
!
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
!
!------- NBRE POINTS INTEGRATIONS ( NPGSR : REDUITE , NPGSN : NORMALE )
!
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
!
!....... LES REELS ( FONCTIONS DE FORMES, DERIVEES ET POIDS )
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
!
!------- CONTRAINTES DE CAUCHY AUX POINTS DE GAUSS
!
    if (option .eq. 'FORC_NODA') then
!
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, &
                    itab=itab)
        icontm = itab(1)
        nbsp = itab(7)
        if (nbsp .ne. npge*nbcou) then
            call utmess('F', 'ELEMENTS_4')
        end if
!
    else if (option .eq. 'REFE_FORC_NODA') then
!
        call terefe('SIGM_REFE', 'MECA_COQUE3D', sigref)
!
    end if
!
!______________________________________________________________________
!
!---- RECUPERATION DES POINTEURS ( E : ECRITURE ) SELON OPTION
!______________________________________________________________________
!
!------- VECTEUR DES FORCES INTERNES
!
    call jevech('PVECTUR', 'E', ivectu)
!
!______________________________________________________________________
!
!---- CARACTERISTIQUES DE COQUE
!
    call jevech('PCACOQU', 'L', jcara)
!
!---- EPAISSEUR TOTALE
!
    eptot = zr(jcara)
!
!---- COORDONNEE MINIMALE SUIVANT L EPAISSEUR
!
    zmin = -eptot/2.d0
!
!---- EPAISSEUR D UNE COUCHE
!
    epais = eptot/nbcou
!
!______________________________________________________________________
!
!---- RECUPERATION DE L ADRESSE DES VARIABLES NODALES TOTALES
!
!---- A L INSTANT MOINS  ( PAS PRECEDENT )
!
    if (option .eq. "FORC_NODA") then
        call jevech('PDEPLAR', 'L', jvDisp)
    else
        call jevech('PDEPLMR', 'L', jvDisp)
    end if
!
!______________________________________________________________________
!
!---- REPERE LOCAUX AUX NOEUDS SUR LA CONFIGURATION INITIALE
!
    call vectan(nb1, nb2, zr(igeom), zr(lzr), vecta, &
                vectn, vectpt)
!
!---- DEPLACEMENT TOTAL AUX NOEUDS DE SERENDIP
!
    call r8inir(8*3, 0.d0, vecu, 1)
!
    do in = 1, nb1
        do ii = 1, 3
            vecu(in, ii) = zr(jvDisp-1+6*(in-1)+ii)
        end do
    end do
!
!---- ROTATION TOTALE AUX NOEUDS
!
    call r8inir(9*3, 0.d0, vecthe, 1)
!
!
!------- EN ACCORD AVEC LA MISE A JOUR DES GRANDES ROTATIONS AUFAURE
!
!------- NOEUDS DE SERENDIP
!
    do in = 1, nb1
        do ii = 1, 3
            vecthe(in, ii) = zr(jvDisp-1+6*(in-1)+ii+3)
        end do
    end do
!
!--------- SUPERNOEUD
!
    do ii = 1, 3
        vecthe(nb2, ii) = zr(jvDisp-1+6*nb1+ii)
    end do
!
!
!---- TRANSFORMEES NORMALES ET MATRICES DE ROTATION AUX NOEUDS
!
    call vectrn(nb2, vectpt, vectn, vecthe, vecnph, &
                blam)
!
!---- VECTEUR PE DES VARIABLES NODALES TOTALES GENERALISEES
!
    call vectpe(nb1, nb2, vecu, vectn, vecnph, &
                vecpe)
!
!______________________________________________________________________
!
!---- INITIALISATION DES OPERATEURS DE DEFORMATION A EXTRAPOLER
!
!---- MEMBRANE REDUIT INCOMPLET
!
    call r8inir(3*51*4, 0.d0, b1mri, 1)
!
    call r8inir(3*51*4, 0.d0, b2mri, 1)
!
!---- SHEAR    REDUIT   COMPLET
!
    call r8inir(2*51*4, 0.d0, b1src, 1)
!
    call r8inir(2*51*4, 0.d0, b2src, 1)
!
! POUR RESI_REFE_RELA
!
    if (option .eq. 'REFE_FORC_NODA') then
!
        call r8inir(51, 0.d0, ftemp, 1)
!
    end if
!
!
!==== BOUCLE SUR LES COUCHES
!
    do icou = 1, nbcou
!
!======= BOUCLE SUR LES POINTS D INTEGRATION SUR L EPAISSEUR
!
        do inte = 1, npge
!
!---------- POSITION SUR L EPAISSEUR ET POIDS D INTEGRATION
!
            if (inte .eq. 1) then
                zic = zmin+(icou-1)*epais
                coef = 1.d0/3.d0
            else if (inte .eq. 2) then
                zic = zmin+epais/2.d0+(icou-1)*epais
                coef = 4.d0/3.d0
            else
                zic = zmin+epais+(icou-1)*epais
                coef = 1.d0/3.d0
            end if
!
!---------- COORDONNEE ISOP.  SUR L EPAISSEUR  DIVISEE PAR DEUX
!
            ksi3s2 = zic/epais
!
!========== 1 ERE BOUCLE SUR POINTS INTEGRATION REDUITE SURFACE MOYENNE
!
            do intsr = 1, npgsr
!
                call vectgt(0, nb1, zr(igeom), ksi3s2, intsr, &
                            zr(lzr), epais, vectn, vectg, vectt)
!
                call jacbm1(epais, vectg, vectt, bid33, jm1, &
                            detj)
!
!------------- J1DN1RI ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 0 INCOMPLET
                call jm1dn1(0, 0, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsr, jm1, jdn1ri)
!
!------------- CALCUL DE    DUDXRI ( 9 ) REDUIT INCOMPLET
!
                call promat(jdn1ri, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxri)
!
!+++++++++++++ B1MRI ( 3 , 51 , 4 ) MEMBRANE REDUIT INCOMPLET
!              B2MRI ( 3 , 51 , 4 )
!
                call matbmr(nb1, vectt, dudxri, intsr, jdn1ri, &
                            b1mri, b2mri)
!
!------------- J1DN1RC ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 1 COMPLET
!
                call jm1dn1(0, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsr, jm1, jdn1rc)
!
!------------- CALCUL DE    DUDXRC ( 9 ) REDUIT COMPLET
!
                call promat(jdn1rc, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxrc)
!
!------------- J1DN2RC ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 1 COMPLET
!
                call jm1dn2(0, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsr, vecnph, jm1, &
                            jdn2rc)
!
!+++++++++++++ B1SRC ( 2 , 51 , 4 ) SHEAR REDUIT COMPLET
!              B2SRC ( 2 , 51 , 4 )
!
                call matbsr(nb1, vectt, dudxrc, intsr, jdn1rc, &
                            jdn2rc, b1src, b2src)
!
!========== FIN 1 ERE BOUCLE NPGSR
!
            end do
!
!========== BOUCLE SUR POINTS INTEGRATION NORMALE SURFACE MOYENNE
!
            do intsn = 1, npgsn
!
!
                call vectgt(1, nb1, zr(igeom), ksi3s2, intsn, &
                            zr(lzr), epais, vectn, vectg, vectt)
!
                call jacbm1(epais, vectg, vectt, bid33, jm1, &
                            detj)
!
!------------- J1DN1NC ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 1 COMPLET
!
                call jm1dn1(1, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsn, jm1, jdn1nc)
!
!------------- CALCUL DE     DUDXNC ( 9 ) NORMAL COMPLET
!
                call promat(jdn1nc, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxnc)
!
!------------- J1DN2NC ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 1 COMPLET
!
                call jm1dn2(1, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsn, vecnph, jm1, &
                            jdn2nc)
!
!+++++++++++++ B1MNC ( 3 , 51 ) MEMBRANE NORMAL COMPLET
!              B2MNC ( 3 , 51 )
!
                call matbmn(nb1, vectt, dudxnc, jdn1nc, jdn2nc, &
                            b1mnc, b2mnc)
!
!------------- J1DN1NI ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 0 INCOMPLET
!
                call jm1dn1(1, 0, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsn, jm1, jdn1ni)
!
!------------- CALCUL DE     DUDXNI ( 9 ) NORMAL INCOMPLET
!
                call promat(jdn1ni, 9, 9, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, dudxni)
!
!+++++++++++++ B1MNI ( 3 , 51 ) MEMBRANE NORMAL INCOMPLET
!              B2MNI ( 3 , 51 )
!
                call matbmn(nb1, vectt, dudxni, jdn1ni, jdn1ni, &
                            b1mni, b2mni)
!
!============= B1SU ( 5 , 51 ) SUBSTITUTION TOTAL
!              B2SU ( 5 , 51 ) SUBSTITUTION DIFFERENTIEL
!
                call matbsu(nb1, zr(lzr), npgsr, intsn, b1mnc, &
                            b2mnc, b1mni, b2mni, b1mri, b2mri, &
                            b1src, b2src, b1su, b2su)
!
                if (option .eq. 'FORC_NODA') then
!
!------- CONTRAINTES DE CAUCHY = PK2 AUX POINTS DE GAUSS
!
                    k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
                    stild(1) = zr(icontm-1+k1+1)
                    stild(2) = zr(icontm-1+k1+2)
                    stild(3) = zr(icontm-1+k1+4)
                    stild(4) = zr(icontm-1+k1+5)
                    stild(5) = zr(icontm-1+k1+6)
!
!------------- FINT ( 6 * NB1 + 3 )  =     INTEGRALE  DE
!              ( B2SU ( 5 , 6 * NB1 + 3 ) ) T * STILD ( 5 ) *
!              POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btsig(6*nb1+3, 5, zr(lzr-1+127+intsn-1)*detj*coef, b2su, stild, &
                               zr(ivectu))
!
!------------- VARIABLES INTERNES INACTIVES COMPORTEMENT NON PLASTIQUE
!
                else if (option .eq. 'REFE_FORC_NODA') then
!
!            CALCUL DES FORCES NODALES DE REFERENCE EN AFFECTANT
!            LA VALEUR SIGM_REFE A CHAQUE CMP SUCCESSIVEMENT
!            POUR CHAQUE POINT D'INTEGRATION
!
                    call r8inir(5, 0.d0, sigtmp, 1)
                    call r8inir(51, 0.d0, effint, 1)
!
                    do i = 1, 5
                        sigtmp(i) = sigref
                        call btsig(6*nb1+3, 5, zr(lzr-1+127+intsn-1)*detj*coef, b2su, sigtmp, &
                                   effint)
                        sigtmp(i) = 0.d0
                        do j = 1, 51
                            ftemp(j) = ftemp(j)+abs(effint(j))
                        end do
                    end do
                end if
            end do
        end do
    end do
!
!      ON PREND LA VALEUR MOYENNE DES FORCES NODALES DE REFERENCE
!
    if (option .eq. 'REFE_FORC_NODA') then
        nval = nbcou*npge*npgsn*5
        b_n = to_blas_int(51)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/nval, ftemp, b_incx, zr(ivectu), &
                   b_incy)
        do j = 1, 51
            if (zr(ivectu+j-1) .eq. 0.) then
                kmess(1) = 'COQUE3D'
                kmess(2) = 'SIGM_REFE'
                call utmess('F', 'MECANONLINE5_59', nk=2, valk=kmess)
            end if
        end do
    end if
!
end subroutine
