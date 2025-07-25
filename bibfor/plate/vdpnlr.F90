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
! aslint: disable=W1501
!
subroutine vdpnlr(option, nomte, codret)
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/antisy.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/btdbma.h"
#include "asterfort/btsig.h"
#include "asterfort/gdt.h"
#include "asterfort/hsaco.h"
#include "asterfort/jacbm1.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/jm1dn1.h"
#include "asterfort/jm1dn2.h"
#include "asterfort/jm1dn3.h"
#include "asterfort/matbmn.h"
#include "asterfort/matbmr.h"
#include "asterfort/matbsr.h"
#include "asterfort/matbsu.h"
#include "asterfort/nmcomp.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rogllo.h"
#include "asterfort/tecach.h"
#include "asterfort/tilbar.h"
#include "asterfort/transp.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "asterfort/vectpe.h"
#include "asterfort/vectrn.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
    integer(kind=8) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!     FONCTION  :  CALCUL DES OBJETS ELEMENTS FINIS EN NON LINEAIRE
!                  GEOMETRIQUE AVEC GRANDES ROTATIONS
!                  COQUE_3D
!
!     ARGUMENTS :
!     DONNEES   :      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
!     OPTIONS   :
!                  RIGI_MECA_TANG : MATRICE TANGENTE DE RIGDITE
!                                   PHASE DE PREDICTION AU DEBUT DE
!                                   CHAQUE PAS
!
!                  RAPH_MECA      : CONTRAINTES CAUCHY ET FORCE INTERNE
!                                   SANS MATRICE TANGENTE DE RIGIDITE
!                                   REAC_ITER : 0 DANS STAT_NON_LINE
!
!                  FULL_MECA      : RAPH_MECA + RIGI_MECA_TANG
!                                   ITERATION TYPIQUE DE NEWTON
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: bid33(3, 3)
    integer(kind=8) :: i, j
    integer(kind=8) :: in
    integer(kind=8) :: jd
    integer(kind=8) :: ii, jj
    integer(kind=8) :: k1, jnbspi
    real(kind=8) :: etild(5), stild(5)
    real(kind=8) :: etildm(5)
    real(kind=8) :: eps2d(4), deps2d(4)
    real(kind=8) :: sign(4), sigma(4), dsidep(6, 6)
    real(kind=8) :: detild(5)
    real(kind=8) :: gxz, gyz
    real(kind=8) :: stlis(5, 4)
    real(kind=8) :: bars(9, 9)
    real(kind=8) :: vecni(3), antni(3, 3)
    real(kind=8) :: veczn(27)
    real(kind=8) :: antzi(3, 3)
    real(kind=8) :: rignc(3, 3)
    integer(kind=8) :: igeom, icontp, imatun, ivectu, ivarip, cod
    integer(kind=8) :: icontm, ivarix
    integer(kind=8) :: lzi, lzr, jcara
    integer(kind=8) :: nb1, nb2, ndimv
    real(kind=8) :: matc(5, 5)
    real(kind=8) :: dtild(5, 5)
    integer(kind=8) :: inte, intsr, intsn
    integer(kind=8) :: kntsr
    real(kind=8) :: eptot, kappa, ctor
    integer(kind=8) :: npge, npgsr, npgsn, ksp
    parameter(npge=3)
    real(kind=8) :: vecta(9, 2, 3)
    real(kind=8) :: vectn(9, 3), vectpt(9, 2, 3)
    real(kind=8) :: vecnph(9, 3)
    real(kind=8) :: vecphm(9, 3)
    real(kind=8) :: vectg(2, 3), vectt(3, 3)
    real(kind=8) :: jm1(3, 3), detj
    real(kind=8) :: hsc(5, 9)
    real(kind=8) :: jdn1ri(9, 51), jdn1rc(9, 51)
    real(kind=8) :: jdn1ni(9, 51), jdn1nc(9, 51)
    real(kind=8) :: jdn2rc(9, 51)
    real(kind=8) :: jdn2nc(9, 51)
    real(kind=8) :: jd2rcm(9, 51)
    real(kind=8) :: jd2ncm(9, 51)
    real(kind=8) :: j1dn3(9, 27)
    real(kind=8) :: btild3(5, 27)
    real(kind=8) :: ksi3s2
    integer(kind=8) :: nbcou
    integer(kind=8) :: icou
    real(kind=8) :: zic, zmin, epais, coef
    real(kind=8) :: vrignc(2601), vrigni(2601)
    real(kind=8) :: vrigrc(2601), vrigri(2601)
    real(kind=8) :: knn
    integer(kind=8) :: iup, ium
    real(kind=8) :: b1su(5, 51), b2su(5, 51)
    real(kind=8) :: b1sum(5, 51), b2sum(5, 51)
    real(kind=8) :: b1src(2, 51, 4)
    real(kind=8) :: b2src(2, 51, 4)
    real(kind=8) :: b1srcm(2, 51, 4)
    real(kind=8) :: b2srcm(2, 51, 4)
    real(kind=8) :: b1mnc(3, 51), b1mni(3, 51)
    real(kind=8) :: b2mnc(3, 51), b2mni(3, 51)
    real(kind=8) :: b1mncm(3, 51), b1mnim(3, 51)
    real(kind=8) :: b2mncm(3, 51), b2mnim(3, 51)
    real(kind=8) :: b1mri(3, 51, 4)
    real(kind=8) :: b2mri(3, 51, 4)
    real(kind=8) :: b1mrim(3, 51, 4)
    real(kind=8) :: b2mrim(3, 51, 4)
    real(kind=8) :: dudxri(9), dudxni(9)
    real(kind=8) :: dudxrc(9), dudxnc(9)
    real(kind=8) :: dudrim(9), dudnim(9)
    real(kind=8) :: dudrcm(9), dudncm(9)
    real(kind=8) :: vecu(8, 3), vecthe(9, 3)
    real(kind=8) :: vecum(8, 3), vecthm(9, 3)
    real(kind=8) :: vecpe(51)
    real(kind=8) :: vecpem(51)
    real(kind=8) :: blam(9, 3, 3)
    real(kind=8) :: blamm(9, 3, 3)
    real(kind=8) :: theta(3), thetan
    real(kind=8) :: tmoin1(3, 3), tm1t(3, 3)
    real(kind=8) :: term(3)
    integer(kind=8), parameter :: nbv = 2
    character(len=16), parameter :: nomres(nbv) = (/'E ', 'NU'/)
    integer(kind=8) :: icodre(nbv)
    real(kind=8) :: valres(nbv)
    character(len=32) :: elasKeyword
    integer(kind=8) :: imate, icarcr, iinstm, iinstp, ivarim
    integer(kind=8) :: nbvari, itab(8), lgpg, k2, iret
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: angmas(3)
    character(len=16) :: defo_comp, rela_comp
    aster_logical :: lVect, lMatr, lVari, lSigm
    real(kind=8) :: cisail
    integer(kind=8), parameter :: ndimLdc = 2
    character(len=8), parameter :: typmod(2) = (/"C_PLAN  ", "        "/)
    type(Behaviour_Integ) :: BEHinteg
    character(len=4), parameter :: fami = "MASS"
    blas_int :: b_incx, b_incy, b_n
    character(len=16), pointer :: compor(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    codret = 0

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)

! - Get input fields
    call jevech('PNBSP_I', 'L', jnbspi)
    nbcou = zi(jnbspi-1+1)
    if (nbcou .le. 0) then
        call utmess('F', 'PLATE1_10')
    end if
    call jevech('PMATERC', 'L', imate)
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PMATERC', 'L', imate)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCARCRI', 'L', icarcr)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                itab=itab)
    if (itab(6) .le. 1) then
        lgpg = itab(7)
    else
        lgpg = itab(6)*itab(7)
    end if

! - Don"t use AFFE_CARA_ELEM/MASSIF
    angmas = r8vide()

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(ndimLdc, typmod, option, &
                              compor, zr(icarcr), &
                              zr(iinstm), zr(iinstp), &
                              fami, zi(imate), &
                              BEHinteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Properties of behaviour
    rela_comp = compor(RELA_NAME)
    defo_comp = compor(DEFO)
    read (compor(NVAR), '(I16)') nbvari
!
! - Get elastic properties
!
    call rccoma(zi(imate), 'ELAS', 1, elasKeyword, icodre(1))
    if (elasKeyword .ne. 'ELAS') then
        call utmess('F', 'PLATE1_12', sk=elasKeyword)
    end if
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
!______________________________________________________________________
!
!---- RECUPERATION DES POINTEURS ( E : ECRITURE ) SELON OPTION
!______________________________________________________________________
!

! - Get output fields
    ivarip = ivarim
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
    end if
    ndimv = lgpg*npgsn
    call jevech('PVARIMP', 'L', ivarix)
    b_n = to_blas_int(ndimv)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(ivarix), b_incx, zr(ivarip), b_incy)
!
    if (lMatr) then
!
!------- MATRICE TANGENTE DE RIGIDITE ET INITIALISATION
!
        call jevech('PMATUNS', 'E', imatun)
!
!------- INITIALISATION DES MATRICES GEOMETRIQUES
!
!------- NORMAL   COMPLET ( CONTRAINTES MEMBRANE FLEXION )
!
        call r8inir(51*51, 0.d0, vrignc, 1)
!
!------- NORMAL INCOMPLET ( CONTRAINTES MEMBRANE FLEXION )
!
        call r8inir(51*51, 0.d0, vrigni, 1)
!
!------- REDUIT INCOMPLET ( CONTRAINTES SHEAR            )
!
        call r8inir(51*51, 0.d0, vrigri, 1)
!
!------- REDUIT   COMPLET ( CONTRAINTES SHEAR            )
!
        call r8inir(51*51, 0.d0, vrigrc, 1)
!
!------- INITIALISATION DE VECZN AVANT INTEGRATION
!
        call r8inir(27, 0.d0, veczn, 1)
!
    end if
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
!---- COEFFICIENT DE CORRECTION DU SHEAR
!
    kappa = zr(jcara+3)
!
!---- COEFFICIENT DE RIGIDITE AUTOUR DE LA TRANSFORMEE DE LA NORMALE
!
    ctor = zr(jcara+4)
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
!______________________________________________________________________
!
!---- RECUPERATION DE L ADRESSE DES VARIABLES NODALES TOTALES
!     QUI NE POSSEDE PAS LE MEME SENS POUR LES DEPLACEMENTS
!     ET LES ROTATIONS
!
!---- A L INSTANT MOINS  ( PAS PRECEDENT )
!
    call jevech('PDEPLMR', 'L', ium)
    call jevech('PCONTMR', 'L', icontm)
!
!---- A L INSTANT PLUS  ( DEPUIS LE PAS PRECEDENT PAS PRECEDENT )
!
    call jevech('PDEPLPR', 'L', iup)
!
!______________________________________________________________________
!
!
!---- REPERE LOCAUX AUX NOEUDS SUR LA CONFIGURATION INITIALE
!
    call vectan(nb1, nb2, zr(igeom), zr(lzr), vecta, &
                vectn, vectpt)
!
!---- DEPLACEMENT TOTAL AUX NOEUDS DE SERENDIP
!
    call r8inir(8*3, 0.d0, vecu, 1)
    call r8inir(8*3, 0.d0, vecum, 1)
!
    do in = 1, nb1
        do ii = 1, 3
!
            vecu(in, ii) = zr(ium-1+6*(in-1)+ii)+zr(iup-1+6*(in-1)+ii)
            vecum(in, ii) = zr(ium-1+6*(in-1)+ii)
!
        end do
    end do
!
!---- ROTATION TOTALE AUX NOEUDS
!
    call r8inir(9*3, 0.d0, vecthe, 1)
    call r8inir(9*3, 0.d0, vecthm, 1)
!
    if (defo_comp .eq. 'GROT_GDEP') then
!
!------- EN ACCORD AVEC LA MISE A JOUR DES GRANDES ROTATIONS AUFAURE
!
!------- NOEUD DE SERENDIP
!
        do in = 1, nb1
            do ii = 1, 3
                vecthe(in, ii) = zr(iup-1+6*(in-1)+ii+3)
                vecthm(in, ii) = zr(ium-1+6*(in-1)+ii+3)
            end do
        end do
!
!------- SUPERNOEUD
!
        do ii = 1, 3
            vecthe(nb2, ii) = zr(iup-1+6*(nb1)+ii)
            vecthm(nb2, ii) = zr(ium-1+6*(nb1)+ii)
        end do
!
    else
!
!------- EN ACCORD AVEC LA MISE A JOUR CLASSIQUE DE STAT_NON_LINE
!
!------- NOEUDS DE SERENDIP
!
        do in = 1, nb1
            do ii = 1, 3
                vecthe(in, ii) = zr(ium-1+6*(in-1)+ii+3)+zr(iup-1+6*(in-1)+ii+3)
                vecthm(in, ii) = zr(ium-1+6*(in-1)+ii+3)
            end do
        end do
!
!--------- SUPERNOEUD
!
        do ii = 1, 3
            vecthe(nb2, ii) = zr(ium-1+6*(nb1)+ii)+zr(iup-1+6*(nb1)+ii)
            vecthm(nb2, ii) = zr(ium-1+6*(nb1)+ii)
        end do
!
    end if
!
!---- TRANSFORMEES NORMALES ET MATRICES DE ROTATION AUX NOEUDS
!
    call vectrn(nb2, vectpt, vectn, vecthe, vecnph, &
                blam)
    call vectrn(nb2, vectpt, vectn, vecthm, vecphm, &
                blamm)
!
!---- VECTEUR PE DES VARIABLES NODALES TOTALES GENERALISEES
!
    call vectpe(nb1, nb2, vecu, vectn, vecnph, &
                vecpe)
    call vectpe(nb1, nb2, vecum, vectn, vecphm, &
                vecpem)
!
!______________________________________________________________________
!
!---- INITIALISATION DES OPERATEURS DE DEFORMATION A EXTRAPOLER
!
!---- MEMBRANE REDUIT INCOMPLET
!
    call r8inir(3*51*4, 0.d0, b1mri, 1)
    call r8inir(3*51*4, 0.d0, b1mrim, 1)
!
    call r8inir(3*51*4, 0.d0, b2mri, 1)
    call r8inir(3*51*4, 0.d0, b2mrim, 1)
!
!---- SHEAR    REDUIT   COMPLET
!
    call r8inir(2*51*4, 0.d0, b1src, 1)
    call r8inir(2*51*4, 0.d0, b1srcm, 1)
!
    call r8inir(2*51*4, 0.d0, b2src, 1)
    call r8inir(2*51*4, 0.d0, b2srcm, 1)
!
!---- COMPTEUR DES POINTS D INTEGRATIONS ( EPAISSEUR * SURFACE )
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
!
                zic = zmin+(icou-1)*epais
!
                coef = 1.d0/3.d0
!
            else if (inte .eq. 2) then
!
                zic = zmin+epais/2.d0+(icou-1)*epais
!
                coef = 4.d0/3.d0
!
            else
!
                zic = zmin+epais+(icou-1)*epais
!
                coef = 1.d0/3.d0
!
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
                call promat(jdn1ri, 9, 9, 6*nb1+3, vecpem, &
                            6*nb1+3, 6*nb1+3, 1, dudrim)
!
!+++++++++++++ B1MRI ( 3 , 51 , 4 ) MEMBRANE REDUIT INCOMPLET
!              B2MRI ( 3 , 51 , 4 )
!
                call matbmr(nb1, vectt, dudxri, intsr, jdn1ri, &
                            b1mri, b2mri)
                call matbmr(nb1, vectt, dudrim, intsr, jdn1ri, &
                            b1mrim, b2mrim)
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
                call promat(jdn1rc, 9, 9, 6*nb1+3, vecpem, &
                            6*nb1+3, 6*nb1+3, 1, dudrcm)
!
!------------- J1DN2RC ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 1 COMPLET
!
                call jm1dn2(0, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsr, vecnph, jm1, &
                            jdn2rc)
                call jm1dn2(0, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsr, vecphm, jm1, &
                            jd2rcm)
!
!+++++++++++++ B1SRC ( 2 , 51 , 4 ) SHEAR REDUIT COMPLET
!              B2SRC ( 2 , 51 , 4 )
!
                call matbsr(nb1, vectt, dudxrc, intsr, jdn1rc, &
                            jdn2rc, b1src, b2src)
                call matbsr(nb1, vectt, dudrcm, intsr, jdn1rc, &
                            jd2rcm, b1srcm, b2srcm)
!
!========== FIN 1 ERE BOUCLE NPGSR
!
            end do
!---------- INITIALISATION DES CONTRAINTES A LISSER
!
            if (lMatr) then
                call r8inir(5*4, 0.d0, stlis, 1)
            end if
!
!========== BOUCLE SUR POINTS INTEGRATION NORMALE SURFACE MOYENNE
!
            do intsn = 1, npgsn
!
!C
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
                call promat(jdn1nc, 9, 9, 6*nb1+3, vecpem, &
                            6*nb1+3, 6*nb1+3, 1, dudncm)
!
!------------- J1DN2NC ( 9 , 6 * NB1 + 3 ) INDN = 1 NORMAL
!                                          INDC = 1 COMPLET
!
                call jm1dn2(1, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsn, vecnph, jm1, &
                            jdn2nc)
                call jm1dn2(1, 1, nb1, nb2, zr(lzr), &
                            epais, ksi3s2, intsn, vecphm, jm1, &
                            jd2ncm)
!
!+++++++++++++ B1MNC ( 3 , 51 ) MEMBRANE NORMAL COMPLET
!              B2MNC ( 3 , 51 )
!
                call matbmn(nb1, vectt, dudxnc, jdn1nc, jdn2nc, &
                            b1mnc, b2mnc)
                call matbmn(nb1, vectt, dudncm, jdn1nc, jd2ncm, &
                            b1mncm, b2mncm)
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
                call promat(jdn1ni, 9, 9, 6*nb1+3, vecpem, &
                            6*nb1+3, 6*nb1+3, 1, dudnim)
!
!+++++++++++++ B1MNI ( 3 , 51 ) MEMBRANE NORMAL INCOMPLET
!              B2MNI ( 3 , 51 )
!
                call matbmn(nb1, vectt, dudxni, jdn1ni, jdn1ni, &
                            b1mni, b2mni)
                call matbmn(nb1, vectt, dudnim, jdn1ni, jdn1ni, &
                            b1mnim, b2mnim)
!
!============= B1SU ( 5 , 51 ) SUBSTITUTION TOTAL
!              B2SU ( 5 , 51 ) SUBSTITUTION DIFFERENTIEL
!
                call matbsu(nb1, zr(lzr), npgsr, intsn, b1mnc, &
                            b2mnc, b1mni, b2mni, b1mri, b2mri, &
                            b1src, b2src, b1su, b2su)
                call matbsu(nb1, zr(lzr), npgsr, intsn, b1mncm, &
                            b2mncm, b1mnim, b2mnim, b1mrim, b2mrim, &
                            b1srcm, b2srcm, b1sum, b2sum)
!
!------------- LA  DEFORMATION TOTALE  DE GREEN LAGRANGE ETILD ( 5 )
!
                call promat(b1su, 5, 5, 6*nb1+3, vecpe, &
                            6*nb1+3, 6*nb1+3, 1, etild)
                call promat(b1sum, 5, 5, 6*nb1+3, vecpem, &
                            6*nb1+3, 6*nb1+3, 1, etildm)
!
!
!------------ INCREMENENT DE DEFORMATION
!
                do i = 1, 5
                    detild(i) = etild(i)-etildm(i)
                end do
!
                eps2d(1) = etildm(1)
                eps2d(2) = etildm(2)
                eps2d(3) = 0.d0
                eps2d(4) = etildm(3)/rac2
!
                deps2d(1) = detild(1)
                deps2d(2) = detild(2)
                deps2d(3) = 0.d0
                deps2d(4) = detild(3)/rac2
!
                gxz = etildm(4)+detild(4)
                gyz = etildm(5)+detild(5)
!
                k2 = lgpg*(intsn-1)+(npge*(icou-1)+inte-1)*nbvari
!
!------- CONTRAINTES DE CAUCHY = PK2 AUX POINTS DE GAUSS INSTANT MOINS
!
                k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
                do i = 1, 3
                    sign(i) = zr(icontm-1+k1+i)
                end do
                sign(4) = zr(icontm-1+k1+4)*rac2

! ------------- Index of "sub"-point
                ksp = (icou-1)*npge+inte

! ------------- Set main parameters for behaviour (on point)
                call behaviourSetParaPoin(intsn, ksp, BEHinteg)

! ------------- Integrator
                sigma = 0.d0
                call nmcomp(BEHinteg, &
                            fami, intsn, ksp, ndimLdc, typmod, &
                            zi(imate), compor, zr(icarcr), zr(iinstm), zr(iinstp), &
                            4, eps2d, deps2d, 4, sign, &
                            zr(ivarim+k2), option, angmas, &
                            sigma, zr(ivarip+k2), 36, dsidep, cod)
!
                call rcvalb(fami, intsn, ksp, '+', zi(imate), &
                            ' ', elasKeyword, 0, ' ', [0.d0], &
                            nbv, nomres, valres, icodre, 1)
!
                cisail = valres(1)/(1.d0+valres(2))
!
!           COD=1 : ECHEC INTEGRATION LOI DE COMPORTEMENT
!           COD=3 : C_PLAN DEBORST SIGZZ NON NUL
                if (cod .ne. 0) then
                    if (codret .ne. 1) then
                        codret = cod
                    end if
                end if
!
!------------- LA  CONTRAINTE TOTALE  PK2 STILD ( 5 )
!
!
                if (lMatr) then
!
                    dtild(1, 1) = dsidep(1, 1)
                    dtild(1, 2) = dsidep(1, 2)
                    dtild(1, 3) = dsidep(1, 4)/rac2
                    dtild(1, 4) = 0.d0
                    dtild(1, 5) = 0.d0
!
                    dtild(2, 1) = dsidep(2, 1)
                    dtild(2, 2) = dsidep(2, 2)
                    dtild(2, 3) = dsidep(2, 4)/rac2
                    dtild(2, 4) = 0.d0
                    dtild(2, 5) = 0.d0
!
                    dtild(3, 1) = dsidep(4, 1)/rac2
                    dtild(3, 2) = dsidep(4, 2)/rac2
                    dtild(3, 3) = dsidep(4, 4)/2.d0
                    dtild(3, 4) = 0.d0
                    dtild(3, 5) = 0.d0
!
                    dtild(4, 1) = 0.d0
                    dtild(4, 2) = 0.d0
                    dtild(4, 3) = 0.d0
                    dtild(4, 4) = cisail*kappa/2.d0
                    dtild(4, 5) = 0.d0
!
                    dtild(5, 1) = 0.d0
                    dtild(5, 2) = 0.d0
                    dtild(5, 3) = 0.d0
                    dtild(5, 4) = 0.d0
                    dtild(5, 5) = cisail*kappa/2.d0
!
                    do i = 1, 5
                        do j = 1, 5
                            matc(i, j) = dtild(i, j)
                        end do
                    end do
!
                end if
!
                if (option(1:16) .eq. 'RIGI_MECA_TANG') then
                    stild(1) = sign(1)
                    stild(2) = sign(2)
                    stild(3) = sign(4)/rac2
                else
                    stild(1) = sigma(1)
                    stild(2) = sigma(2)
                    stild(3) = sigma(4)/rac2
                end if
!
                stild(4) = cisail*kappa*gxz/2.d0
                stild(5) = cisail*kappa*gyz/2.d0
!
                if (lSigm) then
!
!------- CONTRAINTES DE CAUCHY = PK2 AUX POINTS DE GAUSS
!
                    k1 = 6*((intsn-1)*npge*nbcou+(icou-1)*npge+inte-1)
                    zr(icontp-1+k1+1) = stild(1)
                    zr(icontp-1+k1+2) = stild(2)
!
                    zr(icontp-1+k1+3) = 0.d0
!
                    zr(icontp-1+k1+4) = stild(3)
!
                    zr(icontp-1+k1+5) = stild(4)
                    zr(icontp-1+k1+6) = stild(5)
                end if
                if (lVect) then
!
!------------- FINT ( 6 * NB1 + 3 )  =     INTEGRALE  DE
!              ( B2SU ( 5 , 6 * NB1 + 3 ) ) T * STILD ( 5 ) *
!              POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btsig(6*nb1+3, 5, zr(lzr-1+127+intsn-1)*detj*coef, b2su, stild, &
                               zr(ivectu))
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!------------- VARIABLES INTERNES INACTIVES COMPORTEMENT NON PLASTIQUE
!
                end if
!
!
                if (lMatr) then
!
!------------- INTEGRATION DES CONTRAINTES LISSEES
!
                    do kntsr = 1, npgsr
!
                        do i = 1, 5
!
                            stlis(i, kntsr) = stlis(i, kntsr)+zr(lzr-1+702+4*(intsn-1)+kntsr)*sti&
                                              &ld(i)*zr(lzr-1+127+intsn-1)
!
                        end do
                    end do
!
!------------- KM ( 6 * NB1 + 3 , 6 * NB1 + 3 )  =     INTEGRALE  DE
!                ( B2SU ( 5 , 6 * NB1 + 3 ) ) T * MATC ( 5 , 5 ) *
!                  B2SU ( 5 , 6 * NB1 + 3 )
!                POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btdbma(b2su, matc, zr(lzr-1+127+intsn-1)*detj*coef, 5, 6*nb1+3, &
                                zr(imatun))
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!RR
!RR   RIGIDITE GEOMETRIQUE NON CLASSIQUE TOUT
!RR
!
!---------- POUR LE TERME NON CLASSIQUE
!           HSC ( 5 , 9 ) = H ( 5 , 6 )  * S ( 6 , 9 )
!
                    call hsaco(vectt, dudxnc, hsc)
!
!---------- CALCUL DE
!           J1DN3( 9 , 3 * NB2 )=JTILDM1( 9 , 9 )*DNDQSI3( 9 , 3 * NB2 )
!
                    call jm1dn3(nb2, zr(lzr), epais, ksi3s2, intsn, &
                                jm1, j1dn3)
!---------- CALCUL DE
!           BTILD3 ( 5 , 27 ) = HSC ( 5 , 9 ) * J1DN3 ( 9 , 3 * NB2 )
!
                    call promat(hsc, 5, 5, 9, j1dn3, &
                                9, 9, 3*nb2, btild3)
!
!---------- VECZN ( 27 )  =     INTEGRALE  DE
!           ( BTILD3 ( 5 , 27 ) ) T * STILD ( 5 ) *
!           POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btsig(3*nb2, 5, zr(lzr-1+127+intsn-1)*detj*coef, btild3, stild, &
                               veczn)
!
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!----------------------------------------------------------------------
!RR
!RR   RIGIDITE GEOMETRIQUE CLASSIQUE MEMBRANE FLEXION
!RR
!
!------------- ANNULATION DU SHEAR
!---------------------------------
!
                    call r8inir(2, 0.d0, stild(4), 1)
!
!------------- BARS ( 9 , 9 )
!
                    call tilbar(stild, vectt, bars)
!
!------------- VRIGNC  ( 6 * NB1 + 3 , 6 * NB1 + 3 )  = INTEGRALE
!              ( JDN2NC ( 9 , 6 * NB1 + 3 ) ) T * BARS   ( 9 , 9 )
!           *                               JDN2NC ( 9 , 6 * NB1 + 3 ) *
!              POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btdbma(jdn2nc, bars, zr(lzr-1+127+intsn-1)*detj*coef, 9, 6*nb1+3, &
                                vrignc)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!------------- VRIGNI  ( 6 * NB1 + 3 , 6 * NB1 + 3 )  = INTEGRALE
!              ( JDN1NI ( 9 , 6 * NB1 + 3 ) ) T * BARS   ( 9 , 9 )
!           *                               JDN1NI ( 9 , 6 * NB1 + 3 ) *
!              POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
                    call btdbma(jdn1ni, bars, zr(lzr-1+127+intsn-1)*detj*coef, 9, 6*nb1+3, &
                                vrigni)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
                end if
!
!========== FIN BOUCLE NPGSN
!
            end do
!
            if (lMatr) then
!
!========== 2 EME BOUCLE SUR POINTS INTEGRATION REDUITE SURFACE MOYENNE
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
!
                    call jm1dn1(0, 0, nb1, nb2, zr(lzr), &
                                epais, ksi3s2, intsr, jm1, jdn1ri)
!
!
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
!
!------------- RESTITUTION DES CONTRAINTES LISSEES MEMBRANE FLEXION
!
                    do i = 1, 3
                        stild(i) = stlis(i, intsr)
                    end do
!
!------------- ANNULATION DU SHEAR
!
                    call r8inir(2, 0.d0, stild(4), 1)
!
!------------- BARS ( 9 , 9 )
!
                    call tilbar(stild, vectt, bars)
!
!------------- VRIGRI  ( 6 * NB1 + 3 , 6 * NB1 + 3 )  = INTEGRALE
!              ( JDN1RI ( 9 , 6 * NB1 + 3 ) ) T * BARS   ( 9 , 9 )
!           *                               JDN1RI ( 9 , 6 * NB1 + 3 ) *
!                                      DETJ * POIDS EPAISSEUR
!DDDDDDDDDDDDD
!------------- PAS D INTEGRATION REDUITE SURFACE MOYENNE
!DDDDDDDDDDDDD
!
                    call btdbma(jdn1ri, bars, detj*coef, 9, 6*nb1+3, &
                                vrigri)
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!------------- J1DN2RC ( 9 , 6 * NB1 + 3 ) INDN = 0 REDUIT
!                                          INDC = 1 COMPLET
!
                    call jm1dn2(0, 1, nb1, nb2, zr(lzr), &
                                epais, ksi3s2, intsr, vecnph, jm1, &
                                jdn2rc)
!
!------------- ANNULATION DE MEMBRANE FLEXION
!
                    call r8inir(3, 0.d0, stild(1), 1)
!
!------------- RESTITUTION DES CONTRAINTES LISSEES DE SHEAR
!
                    do i = 4, 5
                        stild(i) = stlis(i, intsr)
                    end do
!
!------------- BARS ( 9 , 9 )
!
                    call tilbar(stild, vectt, bars)
!
!------------- VRIGRC  ( 6 * NB1 + 3 , 6 * NB1 + 3 )  = INTEGRALE
!              ( JDN2RC ( 9 , 6 * NB1 + 3 ) ) T * BARS   ( 9 , 9 )
!           *                               JDN2RC ( 9 , 6 * NB1 + 3 ) *
!              POIDS SURFACE MOYENNE * DETJ * POIDS EPAISSEUR
!
!DDDDDDDDDDDDD
!------------- PAS D INTEGRATION REDUITE SURFACE MOYENNE
!DDDDDDDDDDDDD
!
                    call btdbma(jdn2rc, bars, detj*coef, 9, 6*nb1+3, &
                                vrigrc)
!
!========== FIN 2 EME BOUCLE NPGSR
!
                end do
!
            end if
!
!-------- FIN BOUCLE NPGE
!
        end do
!
!---- FIN BOUCLE NBCOU
!
    end do
!
    if (lMatr) then
!
!------- AFFECTATION DE LA RIGIDITE GEOMETRIQUE
!
        do jd = 1, (6*nb1+3)*(6*nb1+3)
            zr(imatun-1+jd) = zr(imatun-1+jd)+vrignc(jd)-vrigni(jd)+vrigri(jd)+vrigrc(jd)
!
        end do
!
!------- AFFECTATION DE LA RIGIDITE NON CLASSIQUE RIGNC ( 3 , 3 )
!
        do in = 1, nb2
!
!---------- MATRICE ANTISYMETRIQUE    ANTZI ( 3 , 3 ) AU NOEUD
!
            call antisy(veczn((in-1)*3+1), 1.d0, antzi)
!
!---------- TRANSFOR DE NORMALE ET SA MATRICE ANTISYM AU NOEUD
!
            do ii = 1, 3
                vecni(ii) = vecnph(in, ii)
            end do
!
            call antisy(vecni, 1.d0, antni)
!
!---------- RIGIDITE NON CLASSIQUE RIGN ( 3 , 3 ) NON SYMETRIQUE
!
            call promat(antzi, 3, 3, 3, antni, &
                        3, 3, 3, rignc)
!
        end do
!
!------- ROTATION DE TOUTE LA MATRICE AU REPERE LOCAL
!
        call rogllo(nb1, nb2, zr(imatun), blam, ctor, &
                    knn)
!
    else
!
!++++ MATRICE ELASTIQUE
!
        knn = 0.d0
!
    end if
!
!++++ SECOND MEMBRE DES FORCES INTERIEURES
!
!++++++++ BOUCLE SUR LES NOEUDS DE ROTATION
!
    do in = 1, nb2
!
!+++++++++++ ROTATION AUTOUR DE LA NORMALE INITIALE
!
        do ii = 1, 3
            vecni(ii) = vectn(in, ii)
            theta(ii) = vecthe(in, ii)
        end do
!
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        thetan = ddot(b_n, theta, b_incx, vecni, b_incy)
!
!+++++++++++ MATRICE T MOIUNS 1 DE THETA
!
        call gdt(theta, tmoin1)
!
!+++++++++++ SON TRANSPOSE
!
        call transp(tmoin1, 3, 3, 3, tm1t, &
                    3)
!
!+++++++++++ PRODUIT T MOINS 1 T FOIS VECNI
!
        call promat(tm1t, 3, 3, 3, vecni, &
                    3, 3, 1, term)
!
!
!
!
        if (lMatr) then
!
            if (in .le. nb1) then
!
!-------------- AFFECTATION
!
!-------------- NOEUDS DE SERENDIP
                do jj = 1, 3
                    j = 6*(in-1)+jj+3
                    do ii = 1, 3
                        i = 6*(in-1)+ii+3
                        zr(imatun-1+(6*nb1+3)*(j-1) &
                           +i) = zr(imatun-1+(6*nb1+3)*( &
                                    j-1)+i)+knn*term(ii)*term &
                              (jj)
                    end do
                end do
!
            else
!
!-------------- SUPERNOEUD
                do jj = 1, 3
                    j = 6*nb1+jj
                    do ii = 1, 3
                        i = 6*nb1+ii
                        zr(imatun-1+(6*nb1+3)*(j-1) &
                           +i) = zr(imatun-1+(6*nb1+3)*( &
                                    j-1)+i)+knn*term(ii)*term &
                              (jj)
                    end do
                end do
!
            end if
!
        end if
!
!
!
        if (lVect) then
!
            if (in .le. nb1) then
!
                do ii = 1, 3
!
                    zr(ivectu-1+6*(in-1)+ii+3) = zr(ivectu-1+6*(in-1)+ii+3)+knn*term(ii)*thetan
!
                end do
!
            else
!
                do ii = 1, 3
                    zr(ivectu-1+6*(in-1)+ii) = zr(ivectu-1+6*(in-1)+ii)+knn*term(ii)*thetan
                end do
!
            end if
!
!
        end if
!
    end do
!
!
! FIN
!
end subroutine
