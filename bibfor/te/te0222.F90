! --------------------------------------------------------------------
! Copyright (C) 1991 - 2021 - EDF R&D - www.code-aster.org
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

subroutine te0222(option, nomte)
!
!------------------------------------------------------------------------------------------
! FONCTION REALISEE:
!
!      CALCUL DU TAUX DE RESTITUTION D'ENERGIE ET DES FACTEURS
!      D'INTENSITE DES CONTRAINTES EN ELASTICITE LINEAIRE ET NON LINEAIRE
!
!      ELEMENTS ISOPARAMETRIQUES 2D/3D
!
!      OPTION : 'CALCH_G'          (LOCAL,CHARGES REELLES)
!               'CALCH_G_F'        (LOCAL,CHARGES FONCTIONS)
!               'CALCH_K_G'        (LOCAL,CHARGES REELLES)
!               'CALCH_K_G_F'      (LOCAL,CHARGES FONCTIONS)
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!
!------------------------------------------------------------------------------------------
!
use Behaviour_type
use Behaviour_module
!
implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/chauxi.h"
#include "asterfort/cgverho.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/gbilin.h"
#include "asterfort/gbil3d.h"
#include "asterfort/jevech.h"
#include "jeveux.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmelnl.h"
#include "asterfort/nmgeom.h"
#include "asterfort/nmplru.h"
#include "asterc/r8prem.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/teattr.h"
#include "asterfort/thetapdg.h"
#include "asterfort/utmess.h"
#include "asterfort/vecini.h"
!
! =====================================================================
!                       DECLARATION DES VARIABLES
! =====================================================================
!
    type(Behaviour_Integ) :: BEHinteg
!
    integer           :: i, j, k, j1, j2, m, kk, l
    integer           :: ier, icodre(4), kp, ncmp
    integer           :: ipoids, ivf, idfde, jgano
    integer           :: ndim, nno, nnos, npg, compt
    integer           :: ivites, iaccel, ithet, igthet, icour, ipuls
    integer           :: icomp, igeom, idepl, imate, matcod
    integer           :: iforf, itemps, iepsf, iforc, iepsr
    integer           :: nbvari, irota, ipesa, isigi, isigm
    integer           :: iret, ireth, ibalo, ideg, ilag
    real(kind=8)      :: tcla, tthe, tfor, tini, thet, poids, f(3,3)
    real(kind=8)      :: der(4), energi(2), divt, dsidep(6,6), mu
    real(kind=8)      :: epsi, valpar(4), accele(3), dsigin(6,3)
    real(kind=8)      :: crit(13), e, ecin, tpg(27), tref, absno
    real(kind=8)      :: prod, prod1, prod2, prod3, prod4, puls
    real(kind=8)      :: dtdm(3, 4), dfdm(3,4), dudm(3,4), dvdm(3,4)
    real(kind=8)      :: rbid, rho, om, omo, epsref(6), depsin(6,3), u1(2), u2(2)
    real(kind=8)      :: tgd(20), nu, rac2, eps(6), epsin(6), epsp(6)
    real(kind=8)      :: tgdm(3), sr(3,3), sigl(6), sigin(6), r_axi
    real(kind=8)      :: c1, c2, c3, k1, k2, k3, guv1, guv2, guv3
    real(kind=8)      :: p(3,3), invp(3,3), du1dm(3,4), du2dm(3,4), du3dm(3,4)
    real(kind=8)      :: courb(3, 3, 3), valres(4), alpha, coeff_K1K2, coeff_K3
    real(kind=8)      :: dfvdm(3,4), k3a, ka, lambda, phig, r8bid, rg, th
    real(kind=8)      :: ttrgv, ttrg, u1l(3), u2l(3), u3l(3), tgvdm(3)
    character(len=8)  :: typmod(2), nompar(4), discr
    character(len=4)  :: fami
    character(len=16) :: nomte, option, compor(4), phenom, oprupt, nomres(4)
!
    aster_logical :: axi, cp, fonc, epsini, grand, incr, notelas, lcour, l_not_zero
!
    real(kind=8), pointer :: fno(:) => null()
    real(kind=8), pointer :: epsino(:) => null()
    real(kind=8), pointer :: dfdi(:) => null()
    real(kind=8), pointer :: ffp(:) => null()
!
! =====================================================================
!                       INITIALISATION PARAMETRES
!                       Cas 2D/3D
!                       Option calc_G/calc_K_G
! =====================================================================
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg,&
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!-- Initialisation des champs et paramètres
    call behaviourInit(BEHinteg)
!
!-- Initialisation des paramètres
    epsi = r8prem()
    tcla = 0.d0
    tthe = 0.d0
    tfor = 0.d0
    tini = 0.d0
    rbid = 0.d0
    r_axi = 0.d0
    der(:) = 0.d0
    valpar(:) = 0.d0
    k1 = 0.d0
    k2 = 0.d0
    k3 = 0.d0
    rac2 = sqrt(2.d0)
    ivites = 0
    iaccel = 0
    nompar(:) = ' '
    epsini = .false.
    axi = .false.
    cp = .false.
    if (ndim .eq. 2) then
        lcour = .false.
    else
        lcour = .true.
    endif
    typmod(1) = '3D'
    typmod(2) = ' '
    oprupt    = 'RUPTURE'
    nomres(1) = 'E'
    nomres(2) = 'NU'
    nomres(3) = 'ALPHA'
    nomres(4) = 'RHO'
    discr = ' '
!
!-- Maillage quadratique ou linéaire
!    call teattr('S','TYPMA',typma, iret)
!
!-- Nombre de composantes des tenseurs
    ncmp = 2*ndim
!
!-- Allocation dynamique : vecteurs 2D/3D
    AS_ALLOCATE(vr=epsino, size=ncmp*nno)
    AS_ALLOCATE(vr=fno, size=ndim*nno)
    AS_ALLOCATE(vr=dfdi, size=ndim*nno)
    AS_ALLOCATE(vr=ffp, size= nno)
!
!-- Cas 2D
    if (ndim .eq. 2) then
        if (lteatt('AXIS','OUI')) then
            typmod(1) = 'AXIS'
            axi = .true.
        else if (lteatt('C_PLAN','OUI')) then
            typmod(1) = 'C_PLAN'
            cp = .true.
        else if (lteatt('D_PLAN','OUI')) then
            typmod(1) = 'D_PLAN'
        endif
    endif
!
! =====================================================================
!                       CALCUL DE THETA
! =====================================================================
!
!   Recuperation des parametres pour créer theta
    call jevech('PTHETAR', 'L', ithet)  ! champ d'entree
!
!   Recuperation du champ de sortie
    call jevech('PGTHETA', 'E', igthet) ! champ de sortie
!
!   Recuperation de DEG si LEGENDRE ou abscisse si LAGRANGE
    if (ndim .eq. 3) then
!
        call tecach('ONO', 'PDEG', 'L', iret, iad=ideg)
!
        if (iret .eq. 0 ) then
            discr = "LEGENDRE"
        else
            discr = "LINEAIRE"
            call jevech('PLAG', 'L', ilag)
        endif
!
    else
        discr = "2D"
    endif
!
!-- Valeur de theta0(r) nulle sur element : pas de calcul
!-- si LAGRANGE, test sur les abscisses curvilignes
    compt = 0
    do i = 1, nno
        thet = zr(ithet-1+6*(i-1)+1)
        if (thet .lt. epsi) compt = compt+1
    end do
    if (compt .eq. nno) goto 999
!
    compt = 0
    if (discr == "LINEAIRE" ) then
        do i = 1, nno
            absno= zr(ithet-1+6*(i-1)+5)
!---------- Cas fond fermé
            if (zr(ilag) .gt. zr(ilag+1) ) then 
                if (absno .ge. zr(ilag) .and. absno .le. zr(ithet-1+6*(i-1)+6)) then
                    compt = compt+1
                elseif (absno .ge. zr(ilag+1) .and. absno .le. zr(ilag+2)) then
                    compt = compt+1
                endif
            else
                if (absno .ge. zr(ilag) .and. absno .lt. zr(ilag+2)) compt = compt+1
            endif 
        end do
        if (compt .eq. 0) goto 999
    endif
!
! =====================================================================
!                  RECUPERATION DES CHAMPS LOCAUX
! =====================================================================
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCOMPOR', 'L', icomp)
    if (option .eq. 'CALCH_K_G' .or. option .eq. 'CALCH_K_G_F') then
        call jevech('PBASLOR', 'L', ibalo)
        if (ndim .eq. 3) then
            call jevech('PCOURB', 'L', icour)
        endif
!------ Verification coherence RHO <-> PESANTEUR, ROTATION, PULSATION
        if (.not. cgverho(imate)) call utmess('F', 'RUPTURE1_26')
    endif
!
    matcod = zi(imate)
!
!-- Recuperation du champ local (carte) associe au pre-epsi
!-- Ce champ est issu d un chargement pre epsi
    if ( option .eq. 'CALCH_G_F' .or. option .eq. 'CALCH_K_G_F' ) then
        fonc = .true.
        call jevech('PFFVOLU', 'L', iforf)
        call jevech('PTEMPSR', 'L', itemps)
        if (ndim .eq. 3) then
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'Z'
            nompar(4) = 'INST'
            valpar(4) = zr(itemps)
        else
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'INST'
            valpar(3) = zr(itemps)
        endif
        call tecach('ONO', 'PEPSINF', 'L', iret, iad=iepsf)
        if (iepsf .ne. 0) then
            epsini = .true.
            do i = 1, ncmp*nno
               epsino(i) = 0.d0
            end do
        endif
    else
        fonc = .false.
        call jevech('PFRVOLU', 'L', iforc)
        call tecach('ONO', 'PEPSINR', 'L', iret, iad=iepsr)
        if (iepsr .ne. 0) then
            epsini = .true.
            do i = 1, ncmp*nno
               epsino(i) = 0.d0
            end do
        endif
    endif
!
!-- Loi de comportement
    do i = 1, 4
        compor(i)= zk16(icomp+i-1)
    end do
    grand = compor(3).eq.'GROT_GDEP'
    incr = compor(4)(1:9).eq.'COMP_INCR'
    notelas = compor(1).ne.'ELAS'
    read(zk16(icomp+1),'(I16)') nbvari
    call tecach('ONO', 'PPESANR', 'L', iret, iad=ipesa)
    call tecach('ONO', 'PROTATR', 'L', iret, iad=irota)
    call tecach('ONO', 'PSIGINR', 'L', iret, iad=isigi)
!
    if (option .eq. 'CALCH_G' .or. option .eq. 'CALCH_G_F') then
        call tecach('ONO', 'PVITESS', 'L', iret, iad=ivites)
        call tecach('ONO', 'PACCELE', 'L', iret, iad=iaccel)
        if (incr) then
            call jevech('PCONTRR', 'L', isigm)
        endif
    else
!------ Restriction avec option calc_k_g
        if (grand) then
            call utmess('F', 'RUPTURE1_24')
        end if
        if (notelas .or.  incr) then
            if (notelas) then
                call utmess('F', 'RUPTURE1_24')
            end if
        endif
!       Pas de comportement incrémentale avec option K
        incr = .FALSE.
!
!------ Récupération de la pulsation
        call tecach('ONO', 'PPULPRO', 'L', iret, iad=ipuls)
        if (ipuls .ne. 0) then
            puls = zr(ipuls)
        else
            puls = 0.d0
        endif
    endif
!
! =====================================================================
!                     MESSAGES D'ERREURS
! =====================================================================
!
!-- On ne peut avoir simultanement pre-deformations et contraintes init
    if ((isigi.ne.0) .and. epsini) then
        call utmess('F', 'RUPTURE1_20')
    endif
!
! =====================================================================
!   RECUPERATION DES CHARGES ET PRE-DEFORMATIONS (CHARGEMENT PRE-EPSI)
! =====================================================================
!
    if (fonc) then
        do i = 1, nno
            do j = 1, ndim
                valpar(j) = zr(igeom+ndim* (i-1)+j-1)
            end do
            do j = 1, ndim
                kk = ndim* (i-1) + j
                call fointe('FM', zk8(iforf+j-1), ndim + 1, nompar, valpar, fno(kk), ier)
            end do
            if (epsini) then
                do j = 1, ncmp
                    kk = ncmp* (i-1) + j
                    call fointe('FM', zk8(iepsf+j-1), ndim + 1, nompar, valpar, epsino(kk), ier)
                end do
            endif
        end do
    else
        do i = 1, nno
            do j = 1, ndim
                fno(ndim* (i-1)+j) = zr(iforc+ndim* (i-1)+j-1)
            end do
            if (epsini) then
                do j = 1, 3
                    epsino(ncmp* (i-1)+j) = zr(iepsr+ncmp* (i-1)+j-1)
                end do
                if (ndim .eq. 3) then
                    do j = 1, 3
                        epsino(ncmp* (i-1)+j+3) = zr(iepsr+ncmp* (i-1)+j-1+3)*rac2
                    enddo
                else
                    epsino(ncmp* (i-1)+4) = zr(iepsr+ncmp*(i-1)-1+4)*rac2
                endif
            endif
        end do
    endif
!
    if (ivites .ne. 0) then
        call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
        call rcvalb(fami, 1, 1, '+', matcod,&
                    ' ', phenom, 0, ' ', [0.d0],&
                    4, nomres, valres, icodre, 0)
!
        nu = valres(2)
        rho = valres(4)
!
    endif
!
!-- Correction des forces volumiques
    if ((ipesa.ne.0) .or. (irota.ne.0)) then
        call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
        call rcvalb(fami, 1, 1, '+', matcod,&
                    ' ', phenom, 0, ' ', [0.d0],&
                    4, nomres, valres, icodre, 0)
!
        rho = valres(4)
!
        if (ipesa .ne. 0) then
            do i = 1, nno
                do j = 1, ndim
                    kk = ndim* (i-1) + j
                    fno(kk) = fno(kk) + rho*zr(ipesa)*zr(ipesa+j)
                end do
            end do
        endif
        if (irota .ne. 0) then
            om = zr(irota)
            do i = 1, nno
                omo = 0.d0
                do j = 1, ndim
                    omo = omo + zr(irota+j)*zr(igeom+ndim* (i-1)+j-1)
                end do
                do j = 1, ndim
                    kk = ndim* (i-1) + j
                    fno(kk) = fno(kk) + rho*om*om* (zr(igeom+kk-1)- omo*zr(irota+j))
                end do
            end do
        endif
    endif
!
! =====================================================================
!     CALCUL DE LA TEMPERATURE AUX NOEUDS, MISE A 0 ZERO SINON
! =====================================================================
!
!-- Recuperation aux noeuds
    do kp = 1, nno
        call rcvarc(' ', 'TEMP', '+', 'NOEU', kp, 1, tgd(kp), ireth)
        if (ireth .eq. 1) tgd(kp) = 0.d0
    end do
!
!-- Recuperation aux PDG
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, 1, tref, iret)
    if (iret .ne. 0) tref = 0.d0
    do kp = 1, npg
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, 1, tpg(kp), iret)
        if (iret .ne. 0) tpg(kp) = 0.d0
    enddo
!
! =====================================================================
!           BOUCLE PRINCIPALE SUR LES POINTS DE GAUSS
! =====================================================================
!
    do kp = 1, npg
!
        ! ===========================================
        !              INITIALISATION
        ! ===========================================
        l = (kp-1)*nno
        do i = 1, 3
            tgdm(i)   = 0.d0
            accele(i) = 0.d0
            tgvdm(i) = 0.d0
            do j = 1, 3
                sr(i,j) = 0.d0
            end do
            do j = 1, 4
                dudm(i,j) = 0.d0
                du1dm(i,j)= 0.d0
                du2dm(i,j)= 0.d0
                du3dm(i,j)= 0.d0
                dvdm(i,j) = 0.d0
                dtdm(i,j) = 0.d0
                dfdm(i,j) = 0.d0
                dfvdm(i,j) = 0.d0
            end do
        end do
        do i = 1, 6
            sigl (i) = 0.d0
            sigin(i) = 0.d0
            epsin(i) = 0.d0
            epsp(i)  = 0.d0
            eps(i)   = 0.d0
            epsref(i)= 0.d0
            do j = 1, 3
                dsigin(i,j) = 0.d0
                depsin(i,j) = 0.d0
            end do
        end do
        p(:,:) = 0.d0
        invp(:,:) = 0.d0
        courb (:,:,:) = 0.d0
        guv1 = 0.d0
        guv2 = 0.d0
        guv3 = 0.d0
        ttrgv = 0.d0
!
        ! ===========================================
        !       CALCUL DES ELEMENTS GEOMETRIQUES
        ! ===========================================
!
        if (ndim .eq. 2) then
            call nmgeom(ndim, nno, axi, grand, zr(igeom),&
                        kp, ipoids, ivf, idfde, zr(idepl),&
                        .true._1, poids, dfdi, f, eps, r_axi)
        else
            call nmgeom(ndim, nno, .false._1, grand, zr(igeom),&
                        kp, ipoids, ivf, idfde, zr(idepl),&
                        .true._1, poids, dfdi, f, eps, rbid)
        endif
        ! ===========================================
        !      CALCULS AUX POINTS DE GAUSS
        !        --    GRADIENT DE U (DUDM)
        !        --    GRADIENT THETA (DTDM)
        !        --    GRADIENT FORCE(DFDM)
        !        --    GRADIENT TEMPERATURE
        ! ===========================================
!
        do i = 1, nno
            der(1) = dfdi(i)
            der(2) = dfdi(i+nno)
            if (ndim .eq. 3) der(3) = dfdi(i+2*nno)
            der(4) = zr(ivf+l+i-1)
            do j = 1, ndim
                tgdm(j) = tgdm(j) + tgd(i)*der(j)
                do k = 1, ndim
                    dudm(j,k) = dudm(j,k) + zr(idepl+ndim* (i-1)+j-1)* der(k)
                    dfdm(j,k) = dfdm(j,k) + fno(ndim* (i-1)+j)*der(k)
                end do
                if (ivites .ne. 0) then
                    do k = 1, ndim
                        dvdm(j,k) = dvdm(j,k) + zr(ivites+ndim*(i-1)+ j-1)*der(k)
                    end do
                    dvdm(j,4) = dvdm(j,4) + zr(ivites+ndim*(i-1)+j-1)* der(4)
                    accele(j) = accele(j) + zr(iaccel+ndim*(i-1)+j-1)* der(4)
                endif
                dudm(j,4) = dudm(j,4) + zr(idepl+ndim* (i-1)+j-1)*der(4)
                dfdm(j,4) = dfdm(j,4) + fno(ndim* (i-1)+j)*der(4)
            end do
        end do
!
!------ Calcul de theta et de son gadient aux points de gauss : DTDM
        if (ndim .eq. 2) then
!
            call thetapdg(ndim, nno, discr, &
                          zr(ivf + l - 1 + 1:ivf + l - 1 + nno), &
                          dfdi, 1, 1, ithet, dtdm)
!
            if (cp) then
                dudm(3,3)= eps(3)
                if (ivites .ne. 0) then
                    dvdm(3,3)= -nu/(1.d0-nu)*(dvdm(1,1)+dvdm(2,2))
                endif
            endif
!
            if (axi) then
                dudm(3,3)= dudm(1,4)/r_axi
                dtdm(3,3)= dtdm(1,4)/r_axi
                dfdm(3,3)= dfdm(1,4)/r_axi
                if (ivites .ne. 0) then
                    dvdm(3,3)= dvdm(1,4)/r_axi
                endif
            endif
!
!---------- Calcul de la divergence 2D
            divt = dtdm(1,1) + dtdm(2,2) + dtdm(3,3)
        else
!
!           Calcul de DTDM pour LEGENDRE ou LAGRANGE
            call thetapdg(ndim, nno, discr , &
                          zr(ivf + l - 1 + 1:ivf + l - 1 + nno), &
                          dfdi, ideg, ilag, ithet, dtdm)
!
!---------- Calcul de la divergence 3D
            divt = dtdm(1,1) + dtdm(2,2) + dtdm(3,3)
        endif
!
        ! ===========================================
        ! PRE-DEFORMATIONS ET GRADIENT DEPSIN
        ! (seule intervenant dans le calcul de G)
        ! ===========================================
!
        if (epsini) then
            do i = 1, nno
                der(1) = dfdi(i)
                der(2) = dfdi(i+nno)
                if (ndim .eq. 3) der(3) = dfdi(i+2*nno)
                der(4) = zr(ivf+l+i-1)
                do j = 1, ncmp
                    epsin(j) = epsin(j) + epsino(ncmp* (i-1)+j)*der(4)
                end do
                do j = 1, ncmp
                    do k = 1, ndim
                        depsin(j,k) = depsin(j,k) + epsino(ncmp* (i-1) +j)*der(k)
                    end do
                end do
            end do
            do i = 1, ncmp
                eps(i) = eps(i) - epsin(i)
            end do
        endif
!
        ! ===========================================
        ! CALCUL DES CONTRAINTES LAGRANGIENNES SIGL
        ! ET DE L'ENERGIE LIBRE
        ! ===========================================
!
        if (incr) then
!---------- Besoin de garder l'appel nmplru pour le calcul de l'energie
!---------- en presence d'etat initial (nettoyage à completer)
            call nmplru(fami, kp, 1, '+', ndim,&
                        typmod, matcod, compor, rbid, eps,&
                        epsp, rbid, energi)
            do i = 1, 3
                sigl(i) = zr(isigm+ncmp* (kp-1)+i-1)
            end do
            if (ndim .eq. 3) then
                do i = 1, 3
                    sigl(i+3) = zr(isigm+ncmp* (kp-1)+i-1+3)*rac2
                enddo
            else
                sigl(4)= zr(isigm+ncmp*(kp-1)+3)*rac2
            endif
        else
!
            crit(1) = 300
            crit(2) = 0.d0
            crit(3) = 1.d-3
            crit(9) = 300
            crit(8) = 1.d-3
!
            call nmelnl(BEHinteg,&
                        fami, kp, 1, '+',&
                        ndim, typmod, matcod, compor, crit,&
                        oprupt, eps, sigl, rbid, dsidep,&
                        energi)
            call tecach('NNO', 'PCONTGR', 'L', iret, iad=isigm)
            if (iret .eq. 0) then
                call jevech('PCONTGR', 'L', isigm)
                do i = 1, 3
                    sigl(i) = zr(isigm+ncmp* (kp-1)+i-1)
                end do
                if (ndim .eq. 3) then
                    do i = 1, 3
                        sigl(i+3) = zr(isigm+ncmp* (kp-1)+i-1+3)*rac2
                    enddo
                else
                    sigl(4)= zr(isigm+ncmp*(kp-1)+3)*rac2
                endif
            endif
        endif
!
!------ Stockage de sigma et des termes croises
        sr(1,1) = sigl(1)
        sr(2,2) = sigl(2)
        sr(3,3) = sigl(3)
        sr(1,2) = sigl(4)/rac2
        sr(1,3) = sigl(5)/rac2
        sr(2,3) = sigl(6)/rac2
        sr(2,1) = sr(1,2)
        sr(3,1) = sr(1,3)
        sr(3,2) = sr(2,3)
!
        ! ===========================================================
        ! RECUPERATION DE E, NU, ALPHA (dilation thermique), RHO
        ! CALCUL DES PARAMETRES UTILES
        ! ===========================================================
!
        call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
!
        call rcvarc(' ', 'TEMP', '+', 'RIGI', kp,&
                    1, r8bid, iret)
!
        call rcvalb(fami, kp, 1, '+', matcod,&
                    ' ', phenom, 0, ' ', [0.d0],&
                    4, nomres, valres, icodre, 0)
!
        ASSERT(icodre(1)+icodre(2).eq.0)
!
        if (icodre(3) .ne. 0) then
            ASSERT(iret.ne.0)
            valres(3) = 0.d0
        endif
!
        if (icodre(4) .ne. 0) then
            valres(4) = 0.d0
        endif
!
        e = valres(1)
        nu = valres(2)
        alpha = valres(3)
        rho = valres(4)
        mu = e/(2.d0*(1.d0+nu))
        lambda = nu*e/((1.d0+nu)*(1.d0-2.d0*nu))
!------ Terme du a la dialation thermique et temperature
        k3a = alpha * e / (1.d0-2.d0*nu)
        ttrg = tpg(kp) - tref
!
        if (ndim .eq.3 .or. lteatt('D_PLAN','OUI') .or. lteatt('AXIS','OUI') ) then
            ka = 3.d0-4.d0*nu
            coeff_K1K2 = e/(1.d0-nu*nu)
            coeff_K3 = 2.d0 * mu
            c1 = lambda + 2.d0 * mu
            c2 = lambda
            c3 = mu
            th = 1.d0
        else
!------ Contrainte plane
            ka = (3.d0-nu)/(1.d0+nu)
            coeff_K1K2 = e
            coeff_K3 = 0.d0
            c1 = e/(1.d0-nu*nu)
            c2 = nu*c1
            c3 = mu
            th = (1.d0-2.d0*nu)/(1.d0-nu)
        endif
!
        ! ===========================================================
        ! CORRECTIONS LIEES A LA CONTRAINTE INITIALE (SIGM DE CALC_G)
        ! CONTRAINTE, DEFORMATION DE REFERENCE, ENERGIE ELASTIQUE
        ! ===========================================================
!
        if (isigi .ne. 0) then
            do i = 1, nno
                der(1) = dfdi(i)
                der(2) = dfdi(i+nno)
                if (ndim .eq. 3) der(3) = dfdi(i+2*nno)
                der(4) = zr(ivf+l+i-1)
!
!-------------- Calcul de sigma initial
                do j = 1, ncmp
                    sigin(j) = sigin(j) + zr(isigi+ncmp* (i-1)+j-1)* der(4)
                end do
!
!-------------- Calcul du gradient de sigma initiale
                do j = 1, ncmp
                    do k = 1, ndim
                        dsigin(j,k) = dsigin(j,k) + zr(isigi+ncmp* (i- 1)+j-1)*der(k)
                    end do
                end do
            end do
!
!---------- Traitement particulier des termes croises
            do i = 4, ncmp
                sigin(i) = sigin(i)*rac2
                do j = 1, ndim
                    dsigin(i,j) = dsigin(i,j)*rac2
                end do
            end do
!
            epsref(1)=-(1.d0/e)*(sigin(1)-(nu*(sigin(2)+sigin(3))))
            epsref(2)=-(1.d0/e)*(sigin(2)-(nu*(sigin(3)+sigin(1))))
            epsref(3)=-(1.d0/e)*(sigin(3)-(nu*(sigin(1)+sigin(2))))
            epsref(4)=-(1.d0/mu)*sigin(4)
            if (ndim .eq. 3) then
                epsref(5)=-(1.d0/mu)*sigin(5)
                epsref(6)=-(1.d0/mu)*sigin(6)
            endif
!
!---------- Energie elastique (expression wadier)
            do i = 1, ncmp
                energi(1) = energi(1) + (eps(i)-0.5d0*epsref(i))* sigin(i)
            end do
        endif
!
        ! ===========================================================
        !                  CALCUL DE G(THETA)
        ! ===========================================================
!
        ! ===========================================================
        !            TERME THERMOELASTIQUE CLASSIQUE
        !          F.SIG:(GRAD(U).GRAD(THET))-ENER*DIVT
        ! ===========================================================
!
        ecin = 0.d0
        prod3 = 0.d0
        prod4 = 0.d0
        if (ivites .ne. 0) then
            do j1 = 1, ndim
                ecin = ecin + dvdm(j1,4)*dvdm(j1,4)
            end do
            do j1 = 1, ndim
                do j2 = 1, ndim
                    prod3 = prod3 + accele(j1)*dudm(j1,j2)*dtdm(j2,4)
                    prod4 = prod4 + dvdm(j1,4)*dvdm(j1,j2)*dtdm(j2,4)
                end do
            end do
            ecin = 0.5d0*rho*ecin
            prod3 = rho*prod3
            prod4 = rho*prod4
        endif
!
        prod = 0.d0
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do m = 1, 3
                        prod = prod + f(i,j)*sr(j,k)*dudm(i,m)*dtdm(m, k)
                    end do
                end do
            end do
        end do
        prod = prod - ecin*divt + prod3 - prod4
        tcla = tcla + poids* (prod-energi(1)*divt)
!
        ! ===========================================================
        !                     TERME THERMIQUE
        !               -(D(ENER)/DT)(GRAD(T).THETA)
        ! ===========================================================
!
        if (ireth .eq. 0) then
            prod = 0.d0
            do i = 1, ndim
                prod = prod + tgdm(i)*dtdm(i,4)
            end do
            tthe = tthe - poids*prod*energi(2)
        endif
!
        ! ===========================================================
        !                  TERME FORCE VOLUMIQUE
        ! ===========================================================
!
        do i = 1, ndim
            prod = 0.d0
            do j = 1, ndim
                prod = prod + dfdm(i,j)*dtdm(j,4)
            end do
            tfor = tfor + dudm(i,4)* (prod+dfdm(i,4)*divt)*poids
        end do
        ! ===========================================================
        !                     TERME INITIAL
        ! PROD1 LIE A LA CONTRAINTE (EPS-EPSREF):GRAD(SIGIN).THETA
        ! PROD2 LIE A LA PREDEFORMATION SIG:GRAD(EPSIN).THETA
        ! ===========================================================
        if ((isigi.ne.0) .or. epsini) then
            prod1 = 0.d0
            prod2 = 0.d0
            if (isigi .ne. 0) then
                do i = 1, ncmp
                    do j = 1, ndim
                        prod1=prod1-(eps(i)-epsref(i))*dsigin(i,j)*&
                        dtdm(j,4)
                    end do
                end do
            else if (epsini) then
                do i = 1, ncmp
                    do j = 1, ndim
                        prod2=prod2+sigl(i)*depsin(i,j)*dtdm(j,4)
                    end do
                end do
            endif
            tini = tini + (prod1+prod2)*poids
        endif
!
        ! ===========================================================
        !         CALCUL DE K1, K2, K3 : OPTION K_G
        !         UTILISATION DE LA FORME BILINEAIRE
        ! ===========================================================
!
        if ( option .eq. 'CALCH_K_G' .or. option .eq. 'CALCH_K_G_F' ) then
!
!---------- Base locale associée au PDG KP
            do i = 1, nno
                ffp(i)=zr(ivf-1+nno*(kp-1)+i)
            enddo
!
            call coor_cyl(ndim, nno, zr(ibalo), zr(igeom), ffp,&
                          p, invp, rg, phig, l_not_zero)
!
!---------- Recuparation du tenseur de courbure
            if (lcour) then
                call jevech('PCOURB', 'L', icour)
                do  i = 1, ndim
                    do  j = 1, ndim
                        courb(i,1,j)=zr(icour-1+ndim*(i-1)+j)
                        courb(i,2,j)=zr(icour-1+ndim*(i+3-1)+j)
                        courb(i,3,j)=zr(icour-1+ndim*(i+6-1)+j)
                    enddo
                enddo
            endif
!
!---------- Calcul des Champs auxiliaires
            call chauxi(ndim, mu, ka, rg, phig,&
                        invp, lcour, courb, du1dm, du2dm,&
                        du3dm, u1l, u2l, u3l)
!
            if (axi) then
!---------------Champs singuliers dans la base locale
                call vecini(ndim, 0.d0, u1)
                call vecini(ndim, 0.d0, u2)
                do i = 1, ndim
                    do j = 1, ndim
                        u1(i) = u1(i) + p(i,j) * u1l(j)
                        u2(i) = u2(i) + p(i,j) * u2l(j)
                    end do
                end do
!
                du1dm(3,3)= u1(1)/r_axi
                du2dm(3,3)= u2(1)/r_axi
            endif
!
            ! ======================
            !         CAS 3D
            ! ======================
            if (ndim .eq. 3) then
!
                call gbil3d(dudm, du1dm, dtdm, dfdm, dfvdm,&
                            tgdm, tgvdm, ttrg, ttrgv, poids, sigin,&
                            dsigin, epsref, c1, c2, c3, k3a, alpha,&
                            1.d0, rho, puls, guv1)
        !
                call gbil3d(dudm, du2dm, dtdm, dfdm, dfvdm,&
                            tgdm, tgvdm, ttrg, ttrgv, poids, sigin,&
                            dsigin, epsref, c1, c2, c3, k3a, alpha,&
                            1.d0, rho, puls, guv2)
        !
                call gbil3d(dudm, du3dm, dtdm, dfdm, dfvdm,&
                            tgdm, tgvdm, ttrg, ttrgv, poids, sigin,&
                            dsigin, epsref, c1, c2, c3, k3a, alpha,&
                            1.d0, rho, puls, guv3)
!
                k1 = k1 + guv1
                k2 = k2 + guv2
                k3 = k3 + guv3
!
            ! ======================
            !         CAS 2D
            ! ======================
            else
!
                call gbilin(fami, kp, zi(imate), dudm, du1dm,&
                            dtdm, dfdm, tgdm, poids, sigin,&
                            dsigin, epsref, c1,c2, c3, 0.5d0, &
                            th, 1.d0, rho, puls, axi, guv1)

        !
                call gbilin(fami, kp, zi(imate), dudm, du2dm,&
                            dtdm, dfdm, tgdm, poids, sigin,&
                            dsigin, epsref, c1, c2, c3, 0.5d0, &
                            th, 1.d0, rho, puls, axi, guv2)
!
                k1 = k1 + guv1
                k2 = k2 + guv2
!
            endif
        endif
!
!-- Fin boucle points de Gauss
    enddo
!
!-- Assemblage final des termes de G, K* et des K* réduits pour le
!-- calcul de G_IRWIN
    zr(igthet)   = tthe + tcla + tfor + tini
    zr(igthet+1) = k1 * coeff_K1K2
    zr(igthet+2) = k2 * coeff_K1K2
    zr(igthet+3) = k3 * coeff_K3
!
    zr(igthet+4) = k1 * sqrt(coeff_K1K2)
    zr(igthet+5) = k2 * sqrt(coeff_K1K2)
    zr(igthet+6) = k3 * sqrt(coeff_K3)
!
!-- Exit sur valeur de theta nulle sur element
    999 continue
!
    AS_DEALLOCATE(vr=epsino)
    AS_DEALLOCATE(vr=fno)
    AS_DEALLOCATE(vr=dfdi)
    AS_DEALLOCATE(vr=ffp)

end subroutine
