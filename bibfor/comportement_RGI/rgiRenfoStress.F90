! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
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
subroutine rgiRenfoStress(xmat, iadrmat, sigmf6, epstf6, epspt6,&
                          teta1, teta2, dt, ppas, theta, fl3d,&
                          end3d, wpl3, vwpl33, vwpl33t, dt3, dr3, ipzero, &
                          nvarbe, ngf, rc00, var0, varf, sigf6d, ierr1)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!
    implicit none
#include "asterf_types.h"
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/x6x33.h"
#include "asterfort/x33x6.h"
#include "asterfort/b3d_valp33.h"
#include "asterfort/chrep3d.h"
#include "asterfort/goujon3d.h"
#include "asterfort/indice0.h"
#include "asterfort/renfort3d.h"
#include "asterfort/getValVect.h"
#include "asterfort/utmess.h"
    integer, intent(in) :: nvarbe, ngf, iadrmat, ipzero(ngf)
    real(kind=8), intent(in) :: xmat(*), sigmf6(6), epstf6(6), epspt6(6)
    real(kind=8), intent(in) :: var0(*), rc00, teta1, teta2, dt, theta
    real(kind=8), intent(in) :: wpl3(3), vwpl33(3,3), vwpl33t(3,3), dt3(3), dr3(3)
    aster_logical, intent(in) :: end3d, fl3d, ppas
    real(kind=8), intent(out) :: varf(*), sigf6d(6)
    integer, intent(out) :: ierr1
!-----------------------------------------------------------------------
!   local variables
    integer :: nbrenf, nbparr,  ngfba
    parameter(nbrenf=5, nbparr=21, ngfba=65)
    integer :: i, j, k, l, istepbid, iadr, numf, numr, nrenf00
    aster_logical :: errr, plast_seule
    real(kind=8) :: XN(ngfba), BN(ngfba), ANN(ngfba, ngfba+1)
    real(kind=8) :: rhor(nbrenf), deqr(nbrenf), yor(nbrenf), syr(nbrenf), taur(nbrenf)
    real(kind=8) :: vecr(nbrenf, 3), hplr(nbrenf), tor(nbrenf), ekr(nbrenf), skr(nbrenf)
    real(kind=8) :: ATRR(nbrenf), gamr(nbrenf), khir(nbrenf), sprec(nbrenf), spref(nbrenf)
    real(kind=8) :: ttaref(nbrenf), xnr(nbrenf), xmuthr(nbrenf), tokr(nbrenf), yksyr(nbrenf)
    real(kind=8) :: sigrm33(3, 3), sigrf33(3, 3), epsr0(nbrenf), vnorm
    real(kind=8) :: eplr0(nbrenf), sigr0(nbrenf), eprm0(nbrenf), mu_r0(nbrenf)
    real(kind=8) :: spre0(nbrenf), sigrfissp(nbrenf, 3, 3), sigrd33p(3, 3), sigrf33p(3, 3)
    real(kind=8) :: sigrh6p(6), dri, epspmf33(3, 3), epstf33(3, 3)
    real(kind=8) :: rhov, sigrf6p(6), sigrh33(3, 3), sigrh33p(3, 3), sigrm33p(3, 3), sigrm6p(6)
    real(kind=8) :: epsrf(nbrenf), epspmf(nbrenf),sigrf(nbrenf), eprk0(nbrenf), eplrf(nbrenf)
    real(kind=8) :: eps_nl(nbrenf), eprkf(nbrenf), eprmf(nbrenf), sigrh6(6), sigrm6(6)
    real(kind=8) :: sig133(3,3), vsig133(3,3), sig13(3)
    
!-----------------------------------------------------------------------
!   pour ne pas avoir de fluage mettre plast_seule à vrai
    plast_seule = ASTER_FALSE
    
    ASSERT(ngf.eq.ngfba)
    sigrm33(:,:)=0.d0
    sigrf33(:,:)=0.d0
    errr = ASTER_FALSE
    ierr1 = 0
    
!   nombre de renforts anisotropes reparties (min 0, max 5)
    nrenf00=int(xmat(iadrmat))
    ASSERT(nrenf00.le.nbrenf)   
!
!   récupération des données matériaux concernant les armatures
    do i = 1, nrenf00
!       rhor(i)=taux de renforts élastoplastiques dans 3 directions 
!       deqr(i)=diamètre équivalent des renforts           
!       yor(i)=module d Young des renforts passifs
!       syr(i)=limite élastique des renforts
!       taur(i)=cisaillement maxi renfort matrice
!       vecr(i,1:3)=direction projetée axe base fixe 1             
        iadr = iadrmat + nbparr*(i-1)+1
        call getValVect(xmat, rhor(i), deqr(i), yor(i), syr(i),&
                    taur(i), vecr(i,1), vecr(i,2),  vecr(i,3), ind1=iadr)

!       vérif / normalisation vecteur direction             
        vnorm=sqrt(vecr(i,1)**2+vecr(i,2)**2+vecr(i,3)**2)  
        if (vnorm .lt. r8miem()) call utmess('F','COMPOR3_1')
        vecr(i,:)=vecr(i,:)/vnorm

!       hplr(i)=module d ecrouissage cinematique des armatures
!       tor(i)=temps caracteristique de fluage/relaxation             
!       ekr(i)=deformation caracteristique de fluage relaxation             
!       skr(i)=contrainte caracteristique de la loi de relaxation/ fluage             
!       ATRR(i)=xenergie d'activation de reference pour la relaxation             
!       gamr(i)=coeff de couplage thermo-mecanique pour la relaxation             
!       khir(i)=coeff de non linearite mecanique de la relaxation             
!       sprec(i)=precontrainte imposee             
!       spref(i)=precontrainte fin de pas                       
!       ttaref(i)=temperature de reference armature
!       xnr(i)=exposant de loi d'activation thermique
!       xmuthr(i)=taux de chargement à partir duquel l activation therm dépend du chargement
!       tokr(i)=temps caractéristique de kelvin pour les renforts             
!       yksyr(i)=rapport Ekelvin/Eélastique pour les renforts             
        iadr = iadrmat + nbparr*(i-1)+9
        call getValVect(xmat, hplr(i), tor(i), ekr(i), skr(i),&
                    ATRR(i), gamr(i), khir(i),  sprec(i), ttaref(i),&
                    xnr(i), xmuthr(i), tokr(i), yksyr(i), ind1=iadr)
        spref(i)= sprec(i)
        
!       calcul du taux volumique d armature
        rhov= sum(rhor(1:nrenf00))
!       do j = 1, nrenf00
!           rhov=rhov+rhor(j)
!       end do
        if (rhov .gt. 1.d0) then
            call utmess('F', 'COMPOR3_10')
        end if
    end do 

!   Contrainte princales du beton
    call x6x33(sigmf6, sig133)
    call b3d_valp33(sig133, sig13, vsig133)
    varf(156)=MAX(sig13(1),sig13(2),sig13(3))
    varf(158)=MIN(sig13(1),sig13(2),sig13(3))
    varf(157)=(sig13(1)+sig13(2)+sig13(3))-(varf(156)+varf(158))

! - Calcul des contraintes dans les renforts
    if (nrenf00 .ne. 0) then
!       J'enlève le istep dans la ligne de dessous car pas de non local pour le moment
!       .and.((istep.eq.0).or.(istep.eq.2))) then

!       recuperation du tenseur des deformations finales
        call x6x33(epstf6, epstf33)
!       recuperation du tenseur des deformations localisees dans la matrice
        call x6x33(epspt6, epspmf33)
!       calcul des contraintes axiales dans chaque renfort            
        do i = 1, nrenf00
            epsr0(i)=var0(nvarbe+7*(i-1)+1)
            eplr0(i)=var0(nvarbe+7*(i-1)+2)
            sigr0(i)=var0(nvarbe+7*(i-1)+3)
            eprm0(i)=var0(nvarbe+7*(i-1)+4)
            mu_r0(i)=var0(nvarbe+7*(i-1)+5)
            eprk0(i)=var0(nvarbe+7*(i-1)+6)
            spre0(i)=var0(nvarbe+7*(i-1)+7)
!             if(istep.eq.2) then
!              recuperation de la deformation finale non locale
!               eps_nl(i)=var0(NVARFLU3D+NVARSUP3D+(i-1)*NVARENF3D+2)
!              sinon elle est estimée dans renfort3d
!            end if
!           contrainte axiale dans les renforts                 
            call renfort3d(istepbid, nbrenf, i, epstf33, vecr,&
                           epsr0(i), epsrf(i), eplr0(i), eplrf(i), yor(i),&
                           syr(i), sigr0(i), sigrf(i), hplr(i), tor(i),&
                           ekr(i), skr(i), ATRR(i), khir(i), gamr(i),&
                           sprec(i), teta1, teta2, dt, ppas,&
                           theta, eprm0(i), eprmf(i), ttaref(i), rhor(i),&
                           mu_r0(i), fl3d, errr, xnr(i), xmuthr(i),&
                           eprk0(i), eprkf(i), tokr(i), yksyr(i), plast_seule,&
                           ann, xn, bn, ngf, ipzero,&
                           epspmf33, epspmf(i), eps_nl(i), spre0(i), spref(i))
            if (errr) then
!                print*,'erreur dans renfort3d'
                ierr1=1
                goto 999
            end if 
            varf(nvarbe+7*(i-1)+1)=epsrf(i)
            varf(nvarbe+7*(i-1)+2)=eplrf(i)
            varf(nvarbe+7*(i-1)+3)=sigrf(i)
            varf(nvarbe+7*(i-1)+4)=eprmf(i)
            varf(nvarbe+7*(i-1)+5)=mu_r0(i)
            varf(nvarbe+7*(i-1)+6)=eprkf(i)
            varf(nvarbe+7*(i-1)+7)=spref(i)                
!
!            modification de la direction par effet goujon, calcul des      
!            vecteurs forces sur chaque fissure                         
            do j = 1, 3
!                on ne calcule l effet goujon que si la fissure
!                localisee existe
                if (dt3(j) .gt. 0.) then
!                   cas du renfort i traversant la fissure j
                    call goujon3d(end3d, nbrenf, i, j, vecr,&
                                  deqr, rhor, wpl3, vwpl33, vwpl33t,&
                                  rc00, sigrf(i), sigrfissp)
                else 
                    do k = 1, 3
                        sigrfissp(i,j,k)=0.d0
                    end do
                end if 
!                     do k=1,3
!                        print*,'sigfissp','i',i,'j',j,'=',sigrfissp(i,j,k)
!                     end do
!                    read*
            end do
!                 end if                        
!             end do                 
!            tenseur des contraintes dans la matrice dues au renfort
            do j = 1, 3
                do k = 1, 3
                    sigrm33(j,k)=sigrm33(j,k)+rhor(i)*sigrf(i)* &
                    vecr(i,j)*vecr(i,k)
                end do
            end do
        end do 
        if (end3d .and. (.not.plast_seule)) then
!          calcul de la resultante sur chaque fissure en cas 
!          d endommagement (on est en base principale des endommagements)
            do numf = 1, 3
                do k = 1, 3
                    sigrd33p(k,numf)=0.d0
                    if (dt3(numf) .gt. 0.) then
                        do numr = 1, nrenf00
                            sigrd33p(k,numf)=sigrd33p(k,numf) &
                     +sigrfissp(numr,numf,k)*dr3(numf)
                        end do
                    end if
                end do
            end do
!          partie symetrique de la contribution des renforts 
!          dans les fissures localisees
            do k = 1, 3
                do l = k, 3
                    if (k .eq. l) then
                        sigrf33p(k,l)=sigrd33p(k,l)
                    else if (k.ne.l) then
                        sigrf33p(k,l)=0.5d0*(sigrd33p(k,l)+sigrd33p(l,k))
                        sigrf33p(l,k)=sigrf33p(k,l)
                    end if
                end do
            end do
            call x33x6(sigrf33p, sigrf6p)
            call chrep3d(sigrf33, sigrf33p, vwpl33t)
!
!          passage en base fixe pour verif             
!          combinaison des contraintes dans les renforts et la matrice
!          en presence d endommagement, reduction de la contribution des
!          renforts dans la partie non endommagee
!          on se place dans la base principale des endommagements de traction
            call chrep3d(sigrm33p, sigrm33, vwpl33)
            call x33x6(sigrm33p, sigrm6p)
            do i = 1, 6
                call indice0(i, k, l)
                dri=max(dr3(k),dr3(l))
                sigrh6p(i)=sigrm6p(i)*(1.d0-dri)+sigrf6p(i)
            end do
            call x6x33(sigrh6p, sigrh33p)
            call chrep3d(sigrh33, sigrh33p, vwpl33t)
            call x33x6(sigrh33, sigrh6)
!          combinaison avec les contraintes dans la matrice et les 
!          contraintes dans les renforts traversant les fissures localisees
            sigf6d(:)=(1.d0-rhov)*sigmf6(:)+sigrh6(:)
        else
!          combinaison des contraintes dans les renforts et la matrice
!          en l abscence d endommagement
            call x33x6(sigrm33, sigrm6)
            sigf6d(:)=(1.d0-rhov)*sigmf6(:)+sigrm6(:)
        end if 
    else
!        pas de renforts
        sigf6d(:)=sigmf6(:)
    end if 
!   tenseur des contrainte matrice seule et homogeneise  Sbei 
    do i = 1, 6
        varf(nvarbe+7*5+i)= sigmf6(i)
    end do
999 continue
        
end subroutine
