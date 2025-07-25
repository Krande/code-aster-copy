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

subroutine te0027(option, nomte)
! FONCTION REALISEE:
!
!      CALCUL DU TAUX DE RESTITUTION D'ENERGIE ELEMENTAIRE
!      EN ELASTICITE LINEAIRE ET NON LINEAIRE
!      ELEMENTS ISOPARAMETRIQUES 3D
!
!      OPTION : 'CALC_G'          (LOCAL,CHARGES REELLES)
!               'CALC_G_F'        (LOCAL,CHARGES FONCTIONS)
!
! ENTREES  ---> OPTION : OPTION DE CALCUL
!          ---> NOMTE  : NOM DU TYPE ELEMENT
!
! ----------------------------------------------------------------------
!
    use Behaviour_type
    use Behaviour_module
!
    implicit none
!
!
! DECLARATION PARAMETRES D'APPELS
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/nmelnl.h"
#include "asterfort/nmgeom.h"
#include "asterfort/nmplru.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
!
!
! DECLARATION VARIABLES LOCALES
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: igeom, itemps, idepl, imate
    integer(kind=8) :: iepsr, iepsf, isigi, isigm
    integer(kind=8) :: iforc, iforf, ithet, igthet, irota, ipesa, ier
    integer(kind=8) :: jgano, nno, nnos, npg, ncmp
    integer(kind=8) :: i, j, k, kk, l, m, kp, ndim, compt, iret
    integer(kind=8) :: ivites, iaccel, j1, j2, ireth, matcod
!
    real(kind=8) :: epsref(6), epsi, rac2
    real(kind=8) :: dfdi(81), f(3, 3), sr(3, 3)
    real(kind=8) :: eps(6), epsin(6), depsin(6, 3), epsp(6)
    real(kind=8) :: epsino(162), fno(81)
    real(kind=8) :: sigl(6), sigin(6), dsigin(6, 3)
    real(kind=8) :: thet, tgdm(3), tgd(20)
    real(kind=8) :: prod, prod1, prod2, divt, valpar(4)
    real(kind=8) :: tcla, tthe, tfor, tini, poids, rbid
    real(kind=8) :: dudm(3, 4), dfdm(3, 4), dtdm(3, 4), der(4), dvdm(3, 4)
    real(kind=8) :: energi(2), rho(1), om, omo
    real(kind=8) :: ecin, prod3, prod4, accele(3), e(1), nu(1), mu
    type(Behaviour_Integ) :: BEHinteg
!
    aster_logical :: grand, fonc, incr, epsini
!
    integer(kind=8) :: icodre(1)
    character(len=4) :: fami
    character(len=8) :: nompar(4), typmod(2)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: phenom
!
! DEB ------------------------------------------------------------------
!
!
    epsi = r8prem()
    rac2 = sqrt(2.d0)
    epsini = .false.
    typmod(1) = '3D      '
    typmod(2) = ' '

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHinteg)
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    ncmp = 2*ndim
!
! - PAS DE CALCUL DE G POUR LES ELEMENTS OU LA VALEUR DE THETA EST NULLE
    call jevech('PTHETAR', 'L', ithet)
    call jevech('PGTHETA', 'E', igthet)
    tcla = 0.d0
    tthe = 0.d0
    tfor = 0.d0
    tini = 0.d0
    compt = 0
    do i = 1, nno
        thet = 0.d0
        do j = 1, ndim
            thet = thet+abs(zr(ithet+ndim*(i-1)+j-1))
        end do
        if (thet .lt. epsi) compt = compt+1
    end do
!
    if (compt .eq. nno) goto 999
!
    ivites = 0
    iaccel = 0
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PDEPLAR', 'L', idepl)
    call jevech('PMATERC', 'L', imate)
    call jevech('PCOMPOR', 'L', vk16=compor)
    matcod = zi(imate)
! RECUPERATION DU CHAMP LOCAL (CARTE) ASSOCIE AU PRE-EPSI
! CE CHAMP EST ISSU D UN CHARGEMENT PRE-EPSI
    if (option .eq. 'CALC_G_XFEM_F') then
        fonc = .true.
        call jevech('PFFVOLU', 'L', iforf)
        call jevech('PINSTR', 'L', itemps)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
        valpar(4) = zr(itemps)
        call tecach('ONO', 'PEPSINF', 'L', iret, iad=iepsf)
        if (iepsf .ne. 0) epsini = .true.
    else
        fonc = .false.
        call jevech('PFRVOLU', 'L', iforc)
        call tecach('ONO', 'PEPSINR', 'L', iret, iad=iepsr)
        if (iepsr .ne. 0) epsini = .true.
    end if
!
    grand = ASTER_FALSE
    incr = compor(4) (1:9) .eq. 'COMP_INCR'
    if (incr) then
        call jevech('PCONTRR', 'L', isigm)
    end if
    call tecach('ONO', 'PPESANR', 'L', iret, iad=ipesa)
    call tecach('ONO', 'PROTATR', 'L', iret, iad=irota)
    call tecach('ONO', 'PSIGINR', 'L', iret, iad=isigi)
    if (option .eq. 'CALC_G_XFEM' .or. option .eq. 'CALC_G_XFEM_F') then
        call tecach('ONO', 'PVITESS', 'L', iret, iad=ivites)
        call tecach('ONO', 'PACCELE', 'L', iret, iad=iaccel)
    end if
!
    do i = 1, ncmp*nno
        epsino(i) = 0.d0
    end do
!
! =====================================================================
! MESSAGES D'ERREURS
! =====================================================================
!
! ON NE PEUT AVOIR SIMULTANEMENT PRE-DEFORMATIONS ET CONTRAINTES INIT.
    if ((isigi .ne. 0) .and. epsini) then
        call utmess('F', 'RUPTURE1_20')
    end if
!
! =====================================================================
! RECUPERATION DES CHARGES ET PRE-DEFORMATIONS (CHARGEMENT PRE-EPSI)
! =====================================================================
    if (fonc) then
        do i = 1, nno
            do j = 1, ndim
                valpar(j) = zr(igeom+ndim*(i-1)+j-1)
            end do
            do j = 1, ndim
                kk = ndim*(i-1)+j
                call fointe('FM', zk8(iforf+j-1), 4, nompar, valpar, fno(kk), ier)
            end do
            if (epsini) then
                do j = 1, ncmp
                    kk = ncmp*(i-1)+j
                    call fointe('FM', zk8(iepsf+j-1), 4, nompar, valpar, epsino(kk), ier)
                end do
            end if
        end do
    else
        do i = 1, nno
            do j = 1, ndim
                fno(ndim*(i-1)+j) = zr(iforc+ndim*(i-1)+j-1)
            end do
            if (epsini) then
                do j = 1, 3
                    epsino(ncmp*(i-1)+j) = zr(iepsr+ncmp*(i-1)+j-1)
                    epsino(ncmp*(i-1)+j+3) = zr(iepsr+ncmp*(i-1)+j-1+3)*rac2
                end do
            end if
        end do
    end if
!
    if (ivites .ne. 0) then
        call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
        call rcvalb('RIGI', 1, 1, '+', matcod, &
                    ' ', phenom, 1, ' ', [rbid], &
                    1, 'RHO', rho, icodre(1), 1)
    end if
!
!
! CORRECTION DES FORCES VOLUMIQUES
    if ((ipesa .ne. 0) .or. (irota .ne. 0)) then
        call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
        call rcvalb('RIGI', 1, 1, '+', matcod, &
                    ' ', phenom, 1, ' ', [rbid], &
                    1, 'RHO', rho, icodre(1), 1)
        if (ipesa .ne. 0) then
            do i = 1, nno
                do j = 1, ndim
                    kk = ndim*(i-1)+j
                    fno(kk) = fno(kk)+rho(1)*zr(ipesa)*zr(ipesa+j)
                end do
            end do
        end if
        if (irota .ne. 0) then
            om = zr(irota)
            do i = 1, nno
                omo = 0.d0
                do j = 1, ndim
                    omo = omo+zr(irota+j)*zr(igeom+ndim*(i-1)+j-1)
                end do
                do j = 1, ndim
                    kk = ndim*(i-1)+j
                    fno(kk) = fno(kk)+rho(1)*om*om*(zr(igeom+kk-1)-omo*zr(irota+j))
                end do
            end do
        end if
    end if
!
!
!
! ======================================================================
!
! CALCUL DU GRADIENT DE TEMPERATURE :
! SI LA TEMPERATURE N'EXISTE PAS, ON LUI IMPOSE UNE VALEUR NULLE
    do kp = 1, nno
        call rcvarc(' ', 'TEMP', '+', 'NOEU', kp, 1, tgd(kp), ireth)
        if (ireth .eq. 1) tgd(kp) = 0.d0
    end do
!
!
! ======================================================================
! BOUCLE PRINCIPALE SUR LES POINTS DE GAUSS
! ======================================================================
!
    do kp = 1, npg
!
! INITIALISATIONS
        l = (kp-1)*nno
        do i = 1, 3
            tgdm(i) = 0.d0
            accele(i) = 0.d0
            do j = 1, 3
                sr(i, j) = 0.d0
            end do
            do j = 1, 4
                dudm(i, j) = 0.d0
                dvdm(i, j) = 0.d0
                dtdm(i, j) = 0.d0
                dfdm(i, j) = 0.d0
            end do
        end do
        do i = 1, 6
            sigl(i) = 0.d0
            sigin(i) = 0.d0
            epsin(i) = 0.d0
            epsp(i) = 0.d0
            eps(i) = 0.d0
            epsref(i) = 0.d0
            do j = 1, 3
                dsigin(i, j) = 0.d0
                depsin(i, j) = 0.d0
            end do
        end do
!
! - CALCUL DES ELEMENTS GEOMETRIQUES
!
        call nmgeom(ndim, nno, .false._1, grand, zr(igeom), &
                    kp, ipoids, ivf, idfde, zr(idepl), &
                    .true._1, poids, dfdi, f, eps, &
                    rbid)
!
! - CALCULS DES GRADIENTS DE U (DUDM),THETA (DTDM) ET FORCE(DFDM)
!   DU GRADIENT (TGDM) DE LA TEMPERATURE AUX POINTS DE GAUSS (TG)
!
        do i = 1, nno
            der(1) = dfdi(i)
            der(2) = dfdi(i+nno)
            der(3) = dfdi(i+2*nno)
            der(4) = zr(ivf+l+i-1)
            do j = 1, ndim
                tgdm(j) = tgdm(j)+tgd(i)*der(j)
                do k = 1, ndim
                    dudm(j, k) = dudm(j, k)+zr(idepl+ndim*(i-1)+j-1)*der(k)
                    dtdm(j, k) = dtdm(j, k)+zr(ithet+ndim*(i-1)+j-1)*der(k)
                    dfdm(j, k) = dfdm(j, k)+fno(ndim*(i-1)+j)*der(k)
                end do
                if (ivites .ne. 0) then
                    do k = 1, ndim
                        dvdm(j, k) = dvdm(j, k)+zr(ivites+ndim*(i-1)+j-1)*der(k)
                    end do
                    dvdm(j, 4) = dvdm(j, 4)+zr(ivites+ndim*(i-1)+j-1)*der(4)
                    accele(j) = accele(j)+zr(iaccel+ndim*(i-1)+j-1)*der(4)
                end if
                dudm(j, 4) = dudm(j, 4)+zr(idepl+ndim*(i-1)+j-1)*der(4)
                dtdm(j, 4) = dtdm(j, 4)+zr(ithet+ndim*(i-1)+j-1)*der(4)
                dfdm(j, 4) = dfdm(j, 4)+fno(ndim*(i-1)+j)*der(4)
            end do
        end do
!
! =======================================================
! PRE DEFORMATIONS ET LEUR GRADIENT DEPSIN
! (seule intervenant dans le calcul de G)
! =======================================================
!
        if (epsini) then
            do i = 1, nno
                der(1) = dfdi(i)
                der(2) = dfdi(i+nno)
                der(3) = dfdi(i+2*nno)
                der(4) = zr(ivf+l+i-1)
                do j = 1, ncmp
                    epsin(j) = epsin(j)+epsino(ncmp*(i-1)+j)*der(4)
                end do
                do j = 1, ncmp
                    do k = 1, ndim
                        depsin(j, k) = depsin(j, k)+epsino(ncmp*(i-1)+j)*der(k)
                    end do
                end do
            end do
            do i = 1, ncmp
                eps(i) = eps(i)-epsin(i)
            end do
        end if
!
! =======================================================
! CALCUL DES CONTRAINTES LAGRANDIENNES SIGL ET DE L'ENERGIE LIBRE
! =======================================================
!
        if (incr) then
! BESOIN GARDER APPEL NMPLRU POUR LE CALCUL DE L'ENERGIE EN PRESENCE
! D'ETAT INITIAL MAIS NETTOYAGE A COMPLETER
            call nmplru(fami, kp, 1, '+', ndim, &
                        typmod, matcod, compor, rbid, eps, &
                        epsp, rbid, energi)
            do i = 1, 3
                sigl(i) = zr(isigm+ncmp*(kp-1)+i-1)
                sigl(i+3) = zr(isigm+ncmp*(kp-1)+i-1+3)*rac2
            end do
        else
!
            call nmelnl(BEHinteg, &
                        fami, kp, 1, &
                        ndim, typmod, matcod, compor, &
                        eps, 0.d0, 0.d0, sigl, energi)
            call tecach('NNO', 'PCONTGR', 'L', iret, iad=isigm)
            if (iret .eq. 0) then
                call jevech('PCONTGR', 'L', isigm)
                do i = 1, 3
                    sigl(i) = zr(isigm+ncmp*(kp-1)+i-1)
                    sigl(i+3) = zr(isigm+ncmp*(kp-1)+i-1+3)*rac2
                end do
            end if
        end if
        divt = dtdm(1, 1)+dtdm(2, 2)+dtdm(3, 3)
!
! =======================================================
! CORRECTIONS LIEES A LA CONTRAINTE INITIALE (SIGM DE CALC_G)
! CONTRAINTE, DEFORMATION DE REFERENCE, ENERGIE ELASTIQUE
! =======================================================
!
        if (isigi .ne. 0) then
            do i = 1, nno
                der(1) = dfdi(i)
                der(2) = dfdi(i+nno)
                der(3) = dfdi(i+2*nno)
                der(4) = zr(ivf+l+i-1)
! CALCUL DE SIGMA INITIAL
                do j = 1, ncmp
                    sigin(j) = sigin(j)+zr(isigi+ncmp*(i-1)+j-1)*der(4)
                end do
! CALCUL DU GRADIENT DE SIGMA INITIAL
                do j = 1, ncmp
                    do k = 1, ndim
                        dsigin(j, k) = dsigin(j, k)+zr(isigi+ncmp*(i-1)+j-1)*der(k)
                    end do
                end do
            end do
!
! TRAITEMENTS PARTICULIERS DES TERMES CROISES
            do i = 4, ncmp
                sigin(i) = sigin(i)*rac2
                do j = 1, ndim
                    dsigin(i, j) = dsigin(i, j)*rac2
                end do
            end do
!
! CALCUL DE LA DEFORMATION DE REFERENCE
            call rccoma(matcod, 'ELAS', 1, phenom, icodre(1))
            call rcvala(matcod, ' ', phenom, 1, ' ', &
                        [rbid], 1, 'NU', nu(1), icodre(1), &
                        1)
            call rcvala(matcod, ' ', phenom, 1, ' ', &
                        [rbid], 1, 'E', e(1), icodre(1), &
                        1)
!
            mu = e(1)/(2.d0*(1.d0+nu(1)))
!
            epsref(1) = -(1.d0/e(1))*(sigin(1)-(nu(1)*(sigin(2)+sigin(3))))
            epsref(2) = -(1.d0/e(1))*(sigin(2)-(nu(1)*(sigin(3)+sigin(1))))
            epsref(3) = -(1.d0/e(1))*(sigin(3)-(nu(1)*(sigin(1)+sigin(2))))
            epsref(4) = -(1.d0/mu)*sigin(4)
            epsref(5) = -(1.d0/mu)*sigin(5)
            epsref(6) = -(1.d0/mu)*sigin(6)
!
! ENERGIE ELASTIQUE (expression WADIER)
            do i = 1, ncmp
                energi(1) = energi(1)+(eps(i)-0.5d0*epsref(i))*sigin(i)
            end do
        end if
!
! =======================================================
! STOCKAGE DE SIGMA ET TRAITEMENTS DES TERMES CROISES
! =======================================================
        sr(1, 1) = sigl(1)
        sr(2, 2) = sigl(2)
        sr(3, 3) = sigl(3)
        sr(1, 2) = sigl(4)/rac2
        sr(1, 3) = sigl(5)/rac2
        sr(2, 3) = sigl(6)/rac2
        sr(2, 1) = sr(1, 2)
        sr(3, 1) = sr(1, 3)
        sr(3, 2) = sr(2, 3)
!
! - CALCUL DE G
!
! =======================================================
! - TERME THERMOELASTIQUE CLASSIQUE F.SIG:(GRAD(U).GRAD(THET))-ENER*DIVT
! =======================================================
        ecin = 0.d0
        prod3 = 0.d0
        prod4 = 0.d0
        if (ivites .ne. 0) then
            do j1 = 1, ndim
                ecin = ecin+dvdm(j1, 4)*dvdm(j1, 4)
            end do
            do j1 = 1, ndim
                do j2 = 1, ndim
                    prod3 = prod3+accele(j1)*dudm(j1, j2)*dtdm(j2, 4)
                    prod4 = prod4+dvdm(j1, 4)*dvdm(j1, j2)*dtdm(j2, 4)
                end do
            end do
            ecin = 0.5d0*rho(1)*ecin
            prod3 = rho(1)*prod3
            prod4 = rho(1)*prod4
        end if
!
        prod = 0.d0
        do i = 1, 3
            do j = 1, 3
                do k = 1, 3
                    do m = 1, 3
                        prod = prod+f(i, j)*sr(j, k)*dudm(i, m)*dtdm(m, k)
                    end do
                end do
            end do
        end do
        prod = prod-ecin*divt+prod3-prod4
        tcla = tcla+poids*(prod-energi(1)*divt)
!
!
! =======================================================
! - TERME THERMIQUE :   -(D(ENER)/DT)(GRAD(T).THETA)
! =======================================================
        if (ireth .eq. 0) then
            prod = 0.d0
            do i = 1, ndim
                prod = prod+tgdm(i)*dtdm(i, 4)
            end do
            tthe = tthe-poids*prod*energi(2)
        end if
! =======================================================
! - TERME FORCE VOLUMIQUE
! =======================================================
        do i = 1, ndim
            prod = 0.d0
            do j = 1, ndim
                prod = prod+dfdm(i, j)*dtdm(j, 4)
            end do
            tfor = tfor+dudm(i, 4)*(prod+dfdm(i, 4)*divt)*poids
        end do
!
! =======================================================
! TERME INITIAL:PROD1 LIE A LA CONTRAINTE (EPS-EPSREF):GRAD(SIGIN).THETA
!               PROD2 LIE A LA PREDEFORMATION SIG:GRAD(EPSIN).THETA
! =======================================================
!
        if ((isigi .ne. 0) .or. epsini) then
            prod1 = 0.d0
            prod2 = 0.d0
            if (isigi .ne. 0) then
                do i = 1, ncmp
                    do j = 1, ndim
                        prod1 = prod1-(eps(i)-epsref(i))*dsigin(i, j)* &
                                dtdm(j, 4)
                    end do
                end do
            else if (epsini) then
                do i = 1, ncmp
                    do j = 1, ndim
                        prod2 = prod2+sigl(i)*depsin(i, j)*dtdm(j, 4)
                    end do
                end do
            end if
            tini = tini+(prod1+prod2)*poids
        end if
! POUR VERIFIER QUE SIGL EST EGAL A E*DUDM+SIGIN (TEST A NU=0)
!        write(6,*)'sigin',sigin(1)
!        write(6,*)'dudm',dudm(1,1
!        write(6,*)'sigl',sigl(1)
!        write(6,*)'E',e(1)
!        write(6,*)'******************'

! ==================================================================
! FIN DE BOUCLE PRINCIPALE SUR LES POINTS DE GAUSS
! ==================================================================
    end do
! EXIT EN CAS DE THETA FISSURE NUL PARTOUT
999 continue
! ASSEMBLAGE FINAL DES TERMES DE G OU DG
    zr(igthet) = tthe+tcla+tfor+tini
!
end subroutine
