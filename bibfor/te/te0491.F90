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

subroutine te0491(option, nomte)
!.......................................................................
    implicit none
!
!
!     BUTS: .CALCUL DES INDICATEURS GLOBAUX
!            DE PERTE DE PROPORTIONNALITE DU CHARGEMENT.
!            POUR LES ELEMENTS ISOPARAMETRIQUES 3D
!           .CALCUL DES ENERGIES DE DEFORMATION ELASTIQUE ET TOTALE
!
! -----------------------------------------------------------------
!
!  OPTION INDIC_ENER : CALCUL DE  L'INDICATEUR GLOBAL
!  =================   ENERGETIQUE DETERMINE PAR L'EXPRESSION SUIVANTE:
!
!            IE = (SOMME_DOMAINE((1 - PSI(EPS)/OMEGA(EPS,VARI)).DV)/V
!
!        OU  .OMEGA EST LA DENSITE D'ENERGIE TOTALE
!            (I.E. OMEGA = SOMME_0->T(SIGMA:D(EPS)/DT).DTAU
!            .PSI EST LA DENSITE D'ENERGIE ELASTIQUE 'TOTALE'
!            (I.E. ASSOCIEE A LA COURBE DE TRACTION SI ON
!                  CONSIDERAIT LE MATERIAU ELASTIQUE NON-LINEAIRE)
!            .V EST LE VOLUME DU GROUPE DE MAILLES TRAITE
! -----------------------------------------------------------------
!
!  OPTION INDIC_SEUIL : CALCUL DE  L'INDICATEUR GLOBAL
!  ==================   DETERMINE PAR L'EXPRESSION SUIVANTE :
!
!   IS = (SOMME_DOMAINE(1 - ((SIG-X):EPS_PLAST)/((SIG_Y+R)*P)).DV)/V
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!            .X         EST LE TENSEUR DE RAPPEL
!            .EPS_PLAST EST LE TENSEUR DES DEFORMATIONS PLASTIQUES
!            .SIG_Y     EST LA LIMITE D'ELASTICITE
!            .R         EST LA FONCTION D'ECROUISSAGE
!            .P         EST LA DEFORMATION PLASTIQUE CUMULEE
!            .V         EST LE VOLUME DU GROUPE DE MAILLES TRAITE
!
!
! -----------------------------------------------------------------
!
!  OPTION ENEL_ELEM : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  ================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!  EN HPP
!   ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!
!        OU  .SIG       EST LE TENSEUR DES CONTRAINTES DE CAUCHY
!            .D         EST LE TENSEUR DE HOOKE
!
!  EN GRANDES DEFORMATIONS SIMO MIEHE POUR ELAS OU VMIS_ISOT
!   ENERLAS = ENERGIE ELASTIQUE SPECIFIQUE
!           = K/2(0.5(J^2-1)-lnJ)+0.5mu(tr(J^(-2/3)be)-3)
!           SI PRESENCE DE THERMIQUE, ON AJOUTE UNE CORRECTION
!           SPECIFIQUE PRESENTEE DANS LA DOC R
!  EN GRANDES DEFORMATIONS GDEF_LOG
!   ENERELAS = SOMME_VOLUME((T_T*(1/D)*T).DV)
!        OU  .T       EST LE TENSEUR DES CONTRAINTES DU FORMALISME
!            .D         EST LE TENSEUR DE HOOKE
! -----------------------------------------------------------------
!
!  OPTION ENTR_ELEM : CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
!  =================   MODIFIEE DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!   ENELAS =  0.5*Lame*H(tr(EPS))*tr(EPS)**2+mu*SUM(H(Ei)*Ei**2)
!
!        OU  .EPS      EST LE TENSEUR DES DEFORMATIONS ELASTIQUES
!            .Ei       SONT (pour i=1..3) LES DEFORMATIONS PROPRES
!            .H        LA FONCTION D'HEAVISIDE
!
! -----------------------------------------------------------------
!
!  OPTION ENER_TOTALE : CALCUL DE L'ENERGIE DE DEFORMATION TOTALE
!  ==================   DETERMINEE PAR L'EXPRESSION SUIVANTE :
!
!   ENER_TOTALE =  ENELAS + EPLAS
!
!          AVEC : ENELAS =  SOMME_VOLUME((SIG_T*(1/D)*SIG).DV)
!                 ENELAS EST L'ENERGIE DE DEFORMATION ELASTIQUE
!
!           OU  .SIG       EST LE TENSEUR DES CONTRAINTES
!               .D         EST LE TENSEUR DE HOOKE
!
!          ET   : EPLAS = SOMME_VOLUME((R(P))*D(P))
!                 EPLAS EST L'ENERGIE DE DEFORMATION PLASTIQUE
!
!           OU  .P         EST LA DEFORMATION PLASTIQUE CUMULEE
!           ET   R(P) EST CALCULE POUR LES COMPORTEMENTS SUIVANTS :
!                      .VMIS_ISOT_LINE
!                      .VMIS_ISOT_TRAC
!
!   REMARQUE : EN GRANDE DEFORMATION ON INTEGRE SUR LE VOLUME INITIALE
!
!         POUR LES AUTRES COMPORTEMENTS ON S'ARRETE EN ERREUR FATALE
! -----------------------------------------------------------------
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTIONS : 'INDIC_ENER'
!                    'INDIC_SEUIL'
!                    'ENEL_ELEM'
!                    'ENER_TOTALE'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/enelpg.h"
#include "asterfort/eps1mc.h"
#include "asterfort/epsvmc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nmgeom.h"
#include "asterfort/getElemOrientation.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/diago3.h"
!-----------------------------------------------------------------------
    integer(kind=8) :: idconm, idene1, idene2, idepl, ideplm, idepmm
    integer(kind=8) :: idfde, idsig, idsigm, idvari, igau, igeom, imate, icodre(5)
    integer(kind=8) :: ipoids, ivf, jgano, jprol, jvale, mxcmel, nbsgm, itemps
    integer(kind=8) :: nbsig, nbsig2, nbval, nbvari, ndim, nno, nnos, npg1
    integer(kind=8) :: iret, iret1
    integer(kind=8) :: i, jtab(7)
    parameter(mxcmel=162)
    parameter(nbsgm=6)
    real(kind=8) :: airep, c1, c2, deux, deuxmu, dsde, e
    real(kind=8) :: lame, mu, vecp(3, 3), epm(3)
    real(kind=8) :: enelas, eneldv, enelsp, enelto, eplaeq, eplast, epseq
    real(kind=8) :: enelastr, enelpart1, enelpart2, trepstraction, welastr
    real(kind=8) :: omega, p, poids, psi, rp
    real(kind=8) :: rprim, sigeq, sigy, tempg, trepsm, trois, trsig
    real(kind=8) :: un, undemi, untier, volume, welas, wtotal, zero
    real(kind=8) :: valres(5)
    real(kind=8) :: sigma(nbsgm), epsdv(nbsgm)
    real(kind=8) :: epsel(nbsgm), epspla(nbsgm), x(nbsgm)
    real(kind=8) :: epsim(nbsgm), sigmm(nbsgm), delta(nbsgm)
    real(kind=8) :: epsi(nbsgm), epssm(mxcmel), epss(mxcmel)
    real(kind=8) :: angl_naut(3), instan, nharm, integ, integ1
    real(kind=8) :: epsm(mxcmel), integ2, nu, k, indigl
    real(kind=8) :: f(3, 3), r, para_vale, epsbid(6), dfdbid(27*3)
    character(len=4) :: fami
    character(len=16) :: nomres(5)
    character(len=8) :: para_type
    character(len=16) :: nomte, option, optio2, compor(3)
    aster_logical :: grand, axi
!-----------------------------------------------------------------------
!
! ---- INITIALISATIONS :
!
    zero = 0.0d0
    undemi = 0.5d0
    un = 1.0d0
    deux = 2.0d0
    trois = 3.0d0
    untier = 1.0d0/3.0d0
    nharm = zero
    omega = zero
    enelas = zero
    eplast = zero
    welas = zero
    wtotal = zero
    welastr = zero
    psi = zero
    volume = zero
    indigl = zero
    instan = r8vide()
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    axi = .false.
!
! ----NOMBRE DE CONTRAINTES ASSOCIE A L'ELEMENT :
!
    nbsig = nbsigm()
!
! ----RECUPERATION DES COORDONNEES DES CONNECTIVITES :
!
    call jevech('PGEOMER', 'L', igeom)
!
! ----RECUPERATION DU MATERIAU :
!
    call jevech('PMATERC', 'L', imate)
!
! ----RECUPERATION  DES DONNEEES RELATIVES AU REPERE D'ORTHOTROPIE :
    call getElemOrientation(ndim, nno, igeom, angl_naut)
!
! ---- RECUPERATION DU CHAMP DE DEPLACEMENTS AUX NOEUDS  :
!
    call jevech('PDEPLR', 'L', idepl)
!
! ---- RECUPERATION DU CHAMP DE CONTRAINTES AUX POINTS D'INTEGRATION :
!
    call jevech('PCONTPR', 'L', idsig)
!
! ---- RECUPERATION DE L'INSTANT DE CALCUL
!      -----------------------------------
    call tecach('NNO', 'PINSTR', 'L', iret, iad=itemps)
    if (itemps .ne. 0) instan = zr(itemps)
!
! ----RECUPERATION DU TYPE DE COMPORTEMENT  :
!     N'EXISTE PAS EN LINEAIRE
    call tecach('NNO', 'PCOMPOR', 'L', iret, nval=7, &
                itab=jtab)
    compor(1) = 'ELAS'
    compor(2) = ' '
    compor(3) = 'PETIT'
    if (iret .eq. 0) then
        compor(1) = zk16(jtab(1))
        compor(3) = zk16(jtab(1)+2)
    end if
!
!     GRANDES DEFORMATIONS
!
    if ((compor(3) .eq. 'SIMO_MIEHE') .or. (compor(3) .eq. 'GDEF_LOG')) then
        grand = .true.
    else
        grand = .false.
    end if
!
! ---- RECUPERATION DU CHAMP DE DEPLACEMENTS AUX NOEUDS  :
!
    if (option .eq. 'ENER_TOTALE') then
        if (grand) then
            call utmess('F', 'COMPOR1_78', sk=compor(3))
        end if
        call tecach('NNO', 'PDEPLM', 'L', iret, iad=ideplm)
        if (ideplm .ne. 0) then
            call jevech('PDEPLM', 'L', idepmm)
        end if
    end if
!
! ON TESTE LA RECUPERATION DU CHAMP DE CONTRAINTES DU PAS PRECEDENT
!
    if (option .eq. 'ENER_TOTALE') then
        if ((compor(1) (1:9) .ne. 'VMIS_ISOT') .and. (compor(1) (1:4) .ne. 'ELAS')) then
            call tecach('NNO', 'PCONTMR', 'L', iret, iad=idconm)
            if (idconm .ne. 0) then
                call jevech('PCONTMR', 'L', idsigm)
            end if
        end if
    end if
!
! ----   RECUPERATION DU CHAMP DE VARIABLES INTERNES  :
!        N'EXISTE PAS EN LINEAIRE
    call tecach('ONO', 'PVARIPR', 'L', iret, nval=7, &
                itab=jtab)
    if (iret .eq. 0) then
        idvari = jtab(1)
        nbvari = max(jtab(6), 1)*jtab(7)
    else
        nbvari = 0
    end if
!
! -- CALCUL DES DEFORMATIONS TOTALES DANS LE CAS DE
! -- ENERGIE TOTALE A L INSTANT COURANT ET CELUI D AVANT
!
    if ((compor(1) (1:9) .ne. 'VMIS_ISOT') .and. (compor(1) (1:4) .ne. 'ELAS')) then
!
        if (option .eq. 'ENER_TOTALE') then
!
!  CALCUL DE B.U
!
            call eps1mc(nno, ndim, nbsig, npg1, ipoids, &
                        ivf, idfde, zr(igeom), zr(idepl), nharm, &
                        epss)
!
            if (ideplm .ne. 0) then
                call eps1mc(nno, ndim, nbsig, npg1, ipoids, &
                            ivf, idfde, zr(igeom), zr(idepmm), nharm, &
                            epssm)
            end if
        end if
!
    end if
!
! ---- CALCUL DES DEFORMATIONS HORS THERMIQUES CORRESPONDANTES AU
! ---- CHAMP DE DEPLACEMENT I.E. EPSM = EPST - EPSTH
! ---- OU EPST  SONT LES DEFORMATIONS TOTALES
! ----    EPST = B.U
! ---- ET EPSTH SONT LES DEFORMATIONS THERMIQUES
! ----    EPSTH = ALPHA*(T-TREF) :
!
    optio2 = 'EPME_ELGA'
    call epsvmc(fami, nno, ndim, nbsig, npg1, &
                ipoids, ivf, idfde, zr(igeom), zr(idepl), &
                instan, angl_naut, nharm, optio2, epsm)
!
!                      ===========================
!                      =                         =
!                      = OPTION   "INDIC_ENER"   =
!                      = OPTION   "ENEL_ELEM"    =
!                      = OPTION   "ENTR_ELEM"    =
!                      = OPTION   "ENER_TOTALE"  =
!                      =                         =
!                      ===========================
!
    if (option .eq. 'INDIC_ENER' .or. option .eq. 'ENEL_ELEM' .or. option .eq. &
        'ENER_TOTALE' .or. option .eq. 'ENTR_ELEM') then
!
! --- BOUCLE SUR LES POINTS D'INTEGRATION
!
        do igau = 1, npg1
!
            omega = zero
            psi = zero
!
! --- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT :
!
            do i = 1, nbsig
                sigma(i) = zr(idsig+(igau-1)*nbsig+i-1)
            end do
!
!
! --- CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT :
            call nmgeom(3, nno, axi, grand, zr(igeom), &
                        igau, ipoids, ivf, idfde, zr(idepl), &
                        .true._1, poids, dfdbid, f, epsbid, &
                        r)
!
!
! --- CALCUL DE L'ENERGIE ELASTIQUE AU POINT D'INTEGRATION COURANT
!
            call enelpg(fami, zi(imate), instan, igau, angl_naut, &
                        compor, f, sigma, nbvari, &
                        zr(idvari+(igau-1)*nbvari), enelas)
!
!
! --- TRAITEMENT DE L'OPTION ENEL_ELEM :
!
            if (option .eq. 'ENEL_ELEM') then
!
                welas = welas+enelas*poids
!
!  ===============================================
!  = FIN TRAITEMENT DE L'OPTION ENEL_ELEM        =
!  ===============================================
!
                goto 120
!
            end if
!
! --- RECUPERATION DES CARACTERISTIQUES DU MATERIAU :
!
            nomres(1) = 'E'
            nomres(2) = 'NU'
!
            call rcvalb(fami, igau, 1, '+', zi(imate), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        2, nomres, valres, icodre, 2)
!
            e = valres(1)
            nu = valres(2)
!
            deuxmu = e/(un+nu)
            k = untier*e/(un-deux*nu)
            lame = (e*nu)/((1+nu)*(1-2*nu))
            mu = e/(2*(1+nu))
!
! --- CALCUL DES DEFORMATIONS ELASTIQUES EN CONSIDERANT LE MATERIAU ISOTROPE :
! --- EPS_ELAS    = 1/D*SIGMA
! ---             = ((1+NU)/E)*SIGMA-(NU/E)*TRACE(SIGMA) :
!
            c1 = (un+nu)/e
            c2 = nu/e

            trsig = sigma(1)+sigma(2)+sigma(3)
            epsel(1) = c1*sigma(1)-c2*trsig
            epsel(2) = c1*sigma(2)-c2*trsig
            epsel(3) = c1*sigma(3)-c2*trsig
            epsel(4) = c1*sigma(4)
            epsel(5) = c1*sigma(5)
            epsel(6) = c1*sigma(6)

!
! --- PARTIE SPHERIQUE DE L'ENERGIE DE DEFORMATION ELASTIQUE POSITIVE OU NULLE:
!

            trepstraction = epsel(1)+epsel(2)+epsel(3)
            if (trepstraction .le. zero) then
                trepstraction = zero
            end if

            enelpart1 = undemi*lame*trepstraction*trepstraction

!
! --- CALCUL DES DEFORMATIONS ELASTIQUES PRINCIPALES POSITIVES OU NULLES:
            call diago3(epsel, vecp, epm)

            if (epm(1) .le. zero) then
                epm(1) = zero
            end if
            if (epm(2) .le. zero) then
                epm(2) = zero
            end if
            if (epm(3) .le. zero) then
                epm(3) = zero
            end if

!
! --- PARTIE DE L'ENERGIE DE DEFORMATION ELASTIQUE (DIRECTIONS PRINCIPALES) POSITIVE OU NULLE:
!

            enelpart2 = mu*(epm(1)*epm(1)+&
                     &epm(2)*epm(2)+&
                     &epm(3)*epm(3))

!
! --- ENERGIE DE DEFORMATION ELASTIQUE DE TRACTION:
!
            enelastr = enelpart1+enelpart2

! --- TRAITEMENT DE L'OPTION ENEL_ELTR :

            if (option .eq. 'ENTR_ELEM') then
!
                welastr = welastr+enelastr*poids
!
                goto 120
!
!  ===============================================
!  = FIN TRAITEMENT DE L'OPTION ENEL_ELTR        =
!  ===============================================
!
            end if

!
! --- RECUPERATION DES CARACTERISTIQUES DU MATERIAU :
!
            nomres(1) = 'E'
            nomres(2) = 'NU'
!
            call rcvalb(fami, igau, 1, '+', zi(imate), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        2, nomres, valres, icodre, 2)
!
!
            e = valres(1)
            nu = valres(2)
!
!-------------------------------------------------------
!   CACUL DU TERME OMEGA REPRESENTANT L'ENERGIE TOTALE -
!   OMEGA = SOMME_0->T(SIGMA:D(EPS)/DT).DTAU           -
!-------------------------------------------------------
! --- TRAITEMENT DU CAS DE L'ECROUISSAGE LINEAIRE ISOTROPE :
!
            if (compor(1) .eq. 'VMIS_ISOT_LINE') then
!
! --- RECUPERATION DE LA LIMITE D'ELASTICITE SY
! --- ET DE LA PENTE DE LA COURBE DE TRACTION D_SIGM_EPSI :
!
                nomres(1) = 'D_SIGM_EPSI'
                nomres(2) = 'SY'
!
                call rcvalb(fami, igau, 1, '+', zi(imate), &
                            ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                            2, nomres, valres, icodre, 2)
!
                dsde = valres(1)
                sigy = valres(2)
!
! --- RECUPERATION DE LA DEFORMATION PLASTIQUE CUMULEE :
!
                p = zr(idvari+(igau-1)*nbvari+1-1)
!
! --- PENTE DE LA COURBE DE TRACTION DANS LE DIAGRAMME 'REDRESSE' :
!
                rprim = e*dsde/(e-dsde)
!
! --- CONTRAINTE UNIAXIALE SUR LA COURBE DE TRACTION :
!
                rp = sigy+rprim*p
!
! --- TRAVAIL PLASTIQUE 'EQUIVALENT' :
!
                eplast = undemi*(sigy+rp)*p
!
! --- TRAITEMENT DU CAS DE L'ECROUISSAGE NON-LINEAIRE ISOTROPE :
!
            else if (compor(1) .eq. 'VMIS_ISOT_TRAC') then
!
! --- RECUPERATION DE LA COURBE DE TRACTION :
!
                call rcvarc(' ', 'TEMP', '+', fami, igau, &
                            1, tempg, iret1)
                call rctype(zi(imate), 1, 'TEMP', [tempg], para_vale, &
                            para_type)
                if ((para_type(1:4) .eq. 'TEMP') .and. (iret1 .eq. 1)) then
                    call utmess('F', 'COMPOR5_5', sk=para_type)
                end if
                call rctrac(zi(imate), 1, 'SIGM', tempg, jprol, &
                            jvale, nbval, e)
!
! --- RECUPERATION DE LA DEFORMATION PLASTIQUE CUMULEE :
!
                p = zr(idvari+(igau-1)*nbvari+1-1)
!
! --- TRAVAIL PLASTIQUE 'EQUIVALENT' :
!
                call rcfonc('V', 1, jprol, jvale, nbval, &
                            p=p, rp=rp, rprim=rprim, airerp=airep)
!
                eplast = airep
            end if
!
! --- AFFECTATION A OMEGA QUI EST L'ENERGIE TOTALE
! --- DE LA CONTRIBUTION AU POINT D'INTEGRATION DE L'ENERGIE
! --- ELASTIQUE ET DE L'ENERGIE PLASTIQUE :
!
            omega = omega+enelas+eplast
!
! --- TRAITEMENT DE L'OPTION ENER_TOTALE :
!
            if (option .eq. 'ENER_TOTALE') then
!
                if ((compor(1) (1:9) .ne. 'VMIS_ISOT') .and. (compor(1) (1:4) .ne. 'ELAS')) then
!
!             TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION PRECEDENT
!             ON LE CALCULE SEULEMENT DANS LE CAS DE LOI DE COMPORTEMENT
!             NI VMIS_ISOT_ NI ELAS
!
                    if (idconm .ne. 0) then
                        do i = 1, nbsig
                            sigmm(i) = zr(idsigm+(igau-1)*nbsig+i-1)
                        end do
                    end if
!
! ---         TENSEUR DES DEFORMATIONS AU POINT D'INTEGRATION COURANT
!             ON LE CALCULE SEULEMENT DANS LE CAS DE LOI DE COMPORTEMENT
!             NI VMIS_ISOT_ NI ELAS
!
                    do i = 1, nbsig
                        epsi(i) = epss(i+(igau-1)*nbsig)
                    end do
!
! ---         TENSEUR DES DEFORMATIONS AU POINT D'INTEGRATION PRECEDENT
!             ON LE CALCULE SEULEMENT DANS LE CAS DE LOI DE COMPORTEMENT
!             NI VMIS_ISOT_ NI ELAS
!
                    if (ideplm .ne. 0) then
                        do i = 1, nbsig
                            epsim(i) = epssm(i+(igau-1)*nbsig)
                        end do
                    end if
!
                    if ((idconm .ne. 0) .and. (ideplm .ne. 0)) then
                        do i = 1, nbsig
                            delta(i) = epsi(i)-epsim(i)
                        end do
!
!---             CALCUL DES TERMES A SOMMER
!
                        integ1 = sigmm(1)*delta(1)+sigmm(2)*delta(2)+sigmm(3)*delta(3)+2.0d&
                                 &0*sigmm(4)*delta(4)+2.0d0*sigmm(5)*delta(5)+2.0d0*sigmm(6)*&
                                 & delta(6)
!
                        integ2 = sigma(1)*delta(1)+sigma(2)*delta(2)+sigma(3)*delta(3)+2.0d&
                                 &0*sigma(4)*delta(4)+2.0d0*sigma(5)*delta(5)+2.0d0*sigma(6)*&
                                 & delta(6)
!
                        integ = undemi*(integ1+integ2)*poids
!
                    else
!
!---             CAS OU LE NUMERO D ORDRE EST UN
!
                        integ = sigma(1)*epsi(1)+sigma(2)*epsi(2)+sigma(3)*epsi(3)+2.0d0*si&
                                &gma(4)*epsi(4)+2.0d0*sigma(5)*epsi(5)+2.0d0*sigma(6)*epsi(6&
                                &)
!
                        integ = undemi*integ*poids
!
                    end if
!
                    wtotal = wtotal+integ
!
                else
!
!-- DANS LE CAS VMIS_ISOT ON CALCULE L ENERGIE TOTALE
!
                    wtotal = wtotal+(enelas+eplast)*poids
!
                end if
!
!  ===============================================
!  = FIN TRAITEMENT DE L'OPTION ENER_TOTALE      =
!  ===============================================
!
                goto 120
!
            end if
!
!---------------------------------------------------------
!   CACUL DU TERME PSI REPRESENTANT L'ENERGIE ELASTIQUE  -
!   NON-LINEAIRE TOTALE ASSOCIEE A LA COURBE DE TRACTION -
!---------------------------------------------------------
!
! --- CALCUL DE LA DILATATION VOLUMIQUE AU POINT D'INTEGRATION COURANT
!
            trepsm = epsm(1+(igau-1)*nbsig)+epsm(2+(igau-1)*nbsig)+epsm(3+(igau-1)*nbsig)
!
! --- DEVIATEUR DES DEFORMATIONS AU POINT D'INTEGRATION COURANT
!
            epsdv(1) = epsm(1+(igau-1)*nbsig)-untier*trepsm
            epsdv(2) = epsm(2+(igau-1)*nbsig)-untier*trepsm
            epsdv(3) = epsm(3+(igau-1)*nbsig)-untier*trepsm
            epsdv(4) = epsm(4+(igau-1)*nbsig)
            epsdv(5) = epsm(5+(igau-1)*nbsig)
            epsdv(6) = epsm(6+(igau-1)*nbsig)
!
! --- CALCUL DE LA DEFORMATION ELASTIQUE EQUIVALENTE AU
! --- POINT D'INTEGRATION COURANT :
!
            epseq = sqrt( &
                    trois/deux*( &
                    epsdv(1)*epsdv(1)+epsdv(2)*epsdv(2)+epsdv(3)*epsdv(3)+deux*(epsdv(4)*epsdv&
                    &(4)+epsdv(5)*epsdv(5)+epsdv(6)*epsdv(6)) &
                    ) &
                    )
!
! --- COEFFICIENTS DU MATERIAU (CONSTANTE ELASTIQUE DE CISAILLEMENT
! --- DE LAME : MU ET  MODULE ELASTIQUE DE DILATATION : K) :
!
            deuxmu = e/(un+nu)
            k = untier*e/(un-deux*nu)
!
! --- CALCUL DE LA CONTRAINTE ELASTIQUE EQUIVALENTE AU
! --- POINT D'INTEGRATION COURANT :
!
            sigeq = deuxmu*epseq
!
! --- PARTIE SPHERIQUE DE L'ENERGIE DE DEFORMATION ELASTIQUE :
!
            enelsp = undemi*k*trepsm*trepsm
!
! --- TRAITEMENT DU CAS DE L'ECROUISSAGE LINEAIRE ISOTROPE :
!
            if (compor(1) .eq. 'VMIS_ISOT_LINE') then
!
! --- RECUPERATION DE LA LIMITE D'ELASTICITE SY
! --- ET DE LA PENTE DE LA COURBE DE TRACTION D_SIGM_EPSI :
!
                nomres(1) = 'D_SIGM_EPSI'
                nomres(2) = 'SY'
!
                call rcvalb(fami, igau, 1, '+', zi(imate), &
                            ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                            2, nomres, valres, icodre, 2)
!
                dsde = valres(1)
                sigy = valres(2)
!
! --- PENTE DE LA COURBE DE TRACTION DANS LE DIAGRAMME 'REDRESSE'
!
                rprim = e*dsde/(e-dsde)
!
! --- DEFORMATION NON-LINEAIRE CUMULEE EQUIVALENTE :
!
                p = (sigeq-sigy)/(rprim+trois/deux*deuxmu)
                if (p .le. r8prem()) p = zero
!
! --- CONTRAINTE UNIAXIALE SUR LA COURBE DE TRACTION :
!
                rp = sigy+rprim*p
!
! --- TRAVAIL ELASTIQUE NON-LINEAIRE 'EQUIVALENT' :
!
                eplast = undemi*(sigy+rp)*p
!
! --- TRAITEMENT DU CAS DE L'ECROUISSAGE NON-LINEAIRE ISOTROPE :
!
            else if (compor(1) .eq. 'VMIS_ISOT_TRAC') then
!
! --- RECUPERATION DE LA COURBE DE TRACTION :
!
                call rctrac(zi(imate), 1, 'SIGM', tempg, jprol, &
                            jvale, nbval, e)
!
! --- CALCUL DE LA LIMITE ELASTIQUE SIGY :
!
                call rcfonc('S', 1, jprol, jvale, nbval, &
                            sigy=sigy)
!
! --- CALCUL DU TRAVAIL ELASTIQUE NON-LINEAIRE ET DE LA
! --- CONTRAINTE EQUIVALENTE :
!
                call rcfonc('E', 1, jprol, jvale, nbval, &
                            e=e, nu=nu, p=zero, rp=rp, rprim=rprim, &
                            airerp=airep, sieleq=sigeq, dp=p)
!
! --- TRAVAIL ELASTIQUE NON-LINEAIRE 'EQUIVALENT' :
!
                eplast = airep
            else
                call utmess('F', 'ELEMENTS4_2')
            end if
!
! --- PARTIE DEVIATORIQUE DE L'ENERGIE DE DEFORMATION ELASTIQUE
! --- TOTALE 'EQUIVALENTE' (I.E. ASSOCIEE A LA COURBE DE
! --- TRACTION SI ON CONSIDERAIT LE MATERIAU ELASTIQUE
! --- NON-LINEAIRE :
!
            if (p .le. r8prem()) then
                eneldv = epseq*epseq*deuxmu/trois
            else
                eneldv = rp*rp/deuxmu/trois
            end if
!
! --- ENERGIE DE DEFORMATION ELASTIQUE TOTALE AU POINT
! --- D'INTEGRATION COURANT :
!
            enelto = enelsp+eneldv+eplast
!
! --- AFFECTATION A PSI QUI EST L'ENERGIE ELASTIQUE TOTALE
! --- DE LA CONTRIBUTION AU POINT D'INTEGRATION DE CETTE ENERGIE :
!
            psi = psi+enelto
!
! --- VOLUME DE L'ELEMENT :
!
            volume = volume+poids
!
! --- INDICATEUR GLOBAL ENERGETIQUE (NON NORMALISE) :
!
            if (omega .ge. 1d04*r8prem()) then
                indigl = indigl+(un-psi/omega)*poids
            end if
!
120         continue
        end do
!
! ----  RECUPERATION ET AFFECTATION DES GRANDEURS EN SORTIE
! ----  AVEC RESPECTIVEMENT LA VALEUR DE L'INDICATEUR GLOBAL SUR
! ----  L'ELEMENT ET LE VOLUME DE L'ELEMENT POUR L'OPTION
! ----  INDIC_ENER
! ----  AFFECTATION DE L'ENERGIE DE DEFORMATION ELASTIQUE
! ----  ET DE L'ENERGIE DE DEFORMATION TOTALE RESPECTIVEMENT
! ----  POUR LES OPTIONS ENEL_ELEM ET ENER_TOTALE :
!
        if (option .eq. 'INDIC_ENER') then
            call jevech('PENERD1', 'E', idene1)
            zr(idene1) = indigl
            call jevech('PENERD2', 'E', idene2)
            zr(idene2) = volume
        else if (option .eq. 'ENEL_ELEM') then
            call jevech('PENERD1', 'E', idene1)
            zr(idene1) = welas
        else if (option .eq. 'ENTR_ELEM') then
            call jevech('PENTRD1', 'E', idene1)
            zr(idene1) = welastr
        else if (option .eq. 'ENER_TOTALE') then
            call jevech('PENERD1', 'E', idene1)
            zr(idene1) = wtotal
        end if
!
!
    else if (option .eq. 'INDIC_SEUIL') then
!
! ---    BOUCLE SUR LES POINTS D'INTEGRATION
!
        do igau = 1, npg1
!
! ---      RECUPERATION DES CARACTERISTIQUES DU MATERIAU :
!
            nomres(1) = 'E'
            nomres(2) = 'NU'
!
            call rcvalb(fami, igau, 1, '+', zi(imate), &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        2, nomres, valres, icodre, 2)
!
!
            e = valres(1)
            nu = valres(2)
!
! ---      TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT :
!
            do i = 1, nbsig
                sigma(i) = zr(idsig+(igau-1)*nbsig+i-1)
            end do
!
! ---      CALCUL DES DEFORMATIONS ELASTIQUES AU POINT
! ---      D'INTEGRATION COURANT EN CONSIDERANT LE MATERIAU ISOTROPE :
! ---       EPS_ELAS    = 1/D*SIGMA
! ---                   = ((1+NU)/E)*SIGMA-(NU/E)*TRACE(SIGMA) :
!
            c1 = (un+nu)/e
            c2 = nu/e
            trsig = sigma(1)+sigma(2)+sigma(3)
            epsel(1) = c1*sigma(1)-c2*trsig
            epsel(2) = c1*sigma(2)-c2*trsig
            epsel(3) = c1*sigma(3)-c2*trsig
            epsel(4) = c1*sigma(4)
            epsel(5) = c1*sigma(5)
            epsel(6) = c1*sigma(6)
!
! ---      CALCUL DES DEFORMATIONS PLASTIQUES AU POINT
! ---      D'INTEGRATION COURANT
! ---      EPS_PLAST = EPS_TOT - EPS_ELAS - EPSTH
! ---      EPS_PLAST = EPS_HORS_THERMIQUE - EPS_ELAS :
!
            epspla(1) = epsm(1+(igau-1)*nbsig)-epsel(1)
            epspla(2) = epsm(2+(igau-1)*nbsig)-epsel(2)
            epspla(3) = epsm(3+(igau-1)*nbsig)-epsel(3)
            epspla(4) = epsm(4+(igau-1)*nbsig)-epsel(4)
            epspla(5) = epsm(5+(igau-1)*nbsig)-epsel(5)
            epspla(6) = epsm(6+(igau-1)*nbsig)-epsel(6)
!
! ---      CAS DE L'ECROUISSAGE CINEMATIQUE :
! ---      LE TENSEUR A PRENDRE EN CONSIDERATION POUR LE CALCUL
! ---      DU TRAVAIL PLASTIQUE EST SIGMA - X OU X EST LE TENSEUR
! ---      DE RAPPEL :
!
            if (compor(1) .eq. 'VMIS_CINE_LINE') then
                nbsig2 = 7
                do i = 1, nbsig
                    x(i) = zr(idvari+(igau-1)*nbsig2+i-1)
                end do
                do i = 1, nbsig
                    sigma(i) = sigma(i)-x(i)
                end do
            end if
!
! ---      CALCUL DU TRAVAIL PLASTIQUE AU POINT D'INTEGRATION COURANT :
!
            eplast = sigma(1)*epspla(1)+sigma(2)*epspla(2)+sigma(3)*epspla(3)+deux*(sigma(&
                     &4)*epspla(4)+sigma(5)*epspla(5)+sigma(6)*epspla(6))
!
! ---      CALCUL DU TRAVAIL PLASTIQUE EQUIVALENT AU POINT
! ---      D'INTEGRATION COURANT :
!
! ---      TRAITEMENT DU CAS DE L'ECROUISSAGE LINEAIRE  :
!
            if (compor(1) .eq. 'VMIS_ISOT_LINE' .or. compor(1) .eq. 'VMIS_CINE_LINE') then
!
! ---          RECUPERATION DE LA LIMITE D'ELASTICITE SY
! ---          ET DE LA PENTE DE LA COURBE DE TRACTION D_SIGM_EPSI :
!
                nomres(1) = 'D_SIGM_EPSI'
                nomres(2) = 'SY'
!
                call rcvalb(fami, igau, 1, '+', zi(imate), &
                            ' ', 'ECRO_LINE', 0, ' ', [0.d0], &
                            2, nomres, valres, icodre, 2)
!
                dsde = valres(1)
                sigy = valres(2)
!
! ---          CALCUL DE LA DEFORMATION PLASTIQUE EQUIVALENTE :
!
                epseq = sqrt( &
                        trois/deux*( &
                        epspla(1)*epspla(1)+epspla(2)*epspla(2)+epspla(3)*epspla(3)+deux*(eps&
                        &pla(4)*epspla(4)+epspla(5)*epspla(5)+epspla(6)*epspla(6)) &
                        ) &
                        )
!
! ---          DEFORMATION PLASTIQUE CUMULEE :
!
! ---         (TEMPORAIRE POUR L'ECROUISSAGE CINEMATIQUE)
                if (compor(1) .eq. 'VMIS_CINE_LINE') then
                    p = epseq
                else if (compor(1) .eq. 'VMIS_ISOT_LINE') then
                    p = zr(idvari+(igau-1)*nbvari+1-1)
                end if
!
! ---          PENTE DE LA COURBE DE TRACTION DANS LE DIAGRAMME
! ---          'REDRESSE' :
!
                rprim = e*dsde/(e-dsde)
!
! ---          CONTRAINTE UNIAXIALE SUR LA COURBE DE TRACTION :
!
                rp = sigy+rprim*p
!
! ---          TRAVAIL PLASTIQUE 'EQUIVALENT' :
!
                eplaeq = rp*p
!
! ---      TRAITEMENT DU CAS DE L'ECROUISSAGE NON-LINEAIRE ISOTROPE :
!
            else if (compor(1) .eq. 'VMIS_ISOT_TRAC') then
!
! ---          RECUPERATION DE LA COURBE DE TRACTION :
!
                call rcvarc(' ', 'TEMP', '+', fami, igau, &
                            1, tempg, iret1)
                call rctrac(zi(imate), 1, 'SIGM', tempg, jprol, &
                            jvale, nbval, e)
!
! ---      RECUPERATION DE LA DEFORMATION PLASTIQUE CUMULEE :
!
                p = zr(idvari+(igau-1)*nbvari+1-1)
!
! ---          TRAVAIL PLASTIQUE 'EQUIVALENT' :
!
                call rcfonc('V', 1, jprol, jvale, nbval, &
                            p=p, rp=rp, rprim=rprim, airerp=airep)
!
                eplaeq = rp*p
            else
                call utmess('F', 'ELEMENTS4_3')
            end if
!
! ---      CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT :
            call dfdm3d(nno, igau, ipoids, idfde, zr(igeom), &
                        poids)
!
! ---      VOLUME DE L'ELEMENT :
!
            volume = volume+poids
!
! ---      INDICATEUR GLOBAL ENERGETIQUE (NON NORMALISE) :
!
            if (eplaeq .ge. 1d04*r8prem()) then
                indigl = indigl+(un-eplast/eplaeq)*poids
            end if
!
        end do
!
! ----   RECUPERATION ET AFFECTATION DES GRANDEURS EN SORTIE
! ----   AVEC RESPECTIVEMENT LA VALEUR DE L'INDICATEUR GLOBAL SUR
! ----   L'ELEMENT ET LE VOLUME DE L'ELEMENT :
!
        call jevech('PENERD1', 'E', idene1)
        zr(idene1) = indigl
        call jevech('PENERD2', 'E', idene2)
        zr(idene2) = volume
!
    end if
!
end subroutine
