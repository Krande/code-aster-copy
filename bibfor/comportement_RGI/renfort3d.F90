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

subroutine renfort3d(istep,nbrenf,numr,epstf33,vecr,&
                     epsr0,epsrf,eprp0,eprpf,your,&
                     syr,sigr0,sigrf,hplr,tomr,&
                     ekr,skr,ATRref,khir,gamr,&
                     sprec,teta1,teta2,dt,ppas,&
                     theta,eprm0,eprmf,ttaref,rhor,&
                     mu00,fl3d,errr,xnr,xmuthr,&
                     eprk0,eprkf,tokr,kr,plast_seule,&
                     ann,xn,bn,ngf,ipzero,&
                     epspmf33,epspmf,eps_nl,spre0,spref)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
 
implicit none
#include "asterfort/gauss3d.h"   
      
!     calcul de la deformation et de la contrainte axiale dans un renfort
      
      integer i,j
      
      integer istep,nbrenf,numr
      aster_logical :: plast_seule
!     def finale, orientation      
      real(kind=8) :: epstf33(3,3),vecr(nbrenf,3),epspmf33(3,3)
!     def totale, plastique...      
      real(kind=8) :: epsr0,epsrf,epspmf,eprp0,eprpf,your,syr,sigr0,sigrf
      real(kind=8) :: hplr,tomr,ekr,skr,spre0
      real(kind=8) :: ATRref,khir,gamr,sprec
      real(kind=8) :: teta1,teta2,dt,theta
      aster_logical  :: ppas
      real(kind=8) :: eprm0,eprmf,ttaref,rhor,mu00,eprk0,eprkf,tokr,kr,K
      aster_logical :: fl3d,errr
      real(kind=8) :: xnr,xmuthr
      integer ngf,errgauss
      real(kind=8) :: ann(ngf,ngf+1),bn(ngf),xn(ngf),spref
      integer  ipzero(ngf) 

!     variables locales      
      real(kind=8) :: epsrf3(3),gamr3(3)
      real(kind=8) :: epsr,ynorm,Teta,Tetaref
      real(kind=8) :: kt,km,kf,mucr,mu0,ATR1,epsmk,cc0,dsigr,f1,sigeq
      real(kind=8) :: epsemin,dtk,xkk,ake,akk,bk,xmm,ame,amm,bm,dep,deps
      real(kind=8) :: depse1,depse2,depsk1,depsk2,depsm1,depsm2,depsp,ep,epse0
      real(kind=8) :: tol_syr,epsek,hp,tauk,taum,x1
      real(kind=8) :: eps_nl
      parameter (epsemin=1.0d-20,tol_syr=1.00001)  

!     intialisation      
      errr=.false.
     
!   if (istep.eq.0) then
!     calcul des deformations axiales (projection) 
      epsr=0.d0
      do i=1,3
          epsrf3(i)=0.d0
          do j=1,3
              epsrf3(i)=epsrf3(i)+epstf33(i,j)*vecr(numr,j)
          end do
          epsr=epsr+epsrf3(i)*vecr(numr,i)
      end do
!     deformation orthogonale 
      ynorm=0.d0
      do i=1,3
          gamr3(i)=epsrf3(i)-epsr*vecr(numr,i)
          ynorm=ynorm+gamr3(i)**2
      end do
!     tension du renfort en l abscence de localisation eps+gama2/2
      epsrf=epsr+0.5d0*ynorm
!     on met 0 ds la deformation eps_nl car le calcul est local
      eps_nl=0.d0         
!         print*,'renfort3d',numr,epsr
!         call affiche33(epstf33)
!         read*  
!      else if (istep.eq.1) then
!c        1ere etape non local on prepare la source de deformation
!c        non locale (la def d ouverture de fissure plastique de la
!c        matrice)  
!         epspm=0.d0       
!         do i=1,3
!             epspm3(i)=0.d0
!             do j=1,3
!                 epspm3(i)=epspm3(i)+epspmf33(i,j)*vecr(numr,j)
!             end do
!             epspm=epspm+epspm3(i)*vecr(numr,i)
!         end do
!c        deformation plastique orthogonale 
!         ynorm=0.d0
!         do i=1,3
!             gamp3(i)=epspm3(i)-epspm*vecr(numr,i)
!             ynorm=ynorm+gamp3(i)**2
!         end do
!c        tension du renfort en l abscence de localisation eps+gama2/2
!         epspmf=epspm+0.5d0*ynorm
!c        au cas ou toutes les armatures ne seraient pas non locales
!c        la variable sera récupéré telle quelle         
!      else if(istep.eq.2) then
!c        2eme etape non locale on reconstruit la deformation dans le
!c        renfort à partir des deformation locales et non locales
!         epsr=0.d0       
!         do i=1,3
!             epsrf3(i)=0.d0
!             do j=1,3
!c                on retranche la def plastique locale             
!                 epsrf3(i)=epsrf3(i)+(epstf33(i,j)-epspmf33(i,j))
!     #           *vecr(numr,j)
!             end do
!             epsr=epsr+epsrf3(i)*vecr(numr,i)
!         end do
!c        deformation orthogonale 
!         ynorm=0.d0
!         do i=1,3
!             gamr3(i)=epsrf3(i)-epsr*vecr(numr,i)
!             ynorm=ynorm+gamr3(i)**2
!         end do
!c        tension du renfort en l abscence de localisation eps+gama2/2
!         epsrf=epsr+0.5d0*ynorm 
!c        on rajoute la def plastique non locale (supplement de eps du
!c        au glissement renfort - matrice)
!         epsrf=epsrf+eps_nl       
!c         print*,'renfort3d',numr,epstf33(numr,numr),epsr
!c         read* 
!      else
!         print*,'istep imprevu ds renfort 3d'
!         errr=.true.
!         return         
!      end if
         
      
!      if((istep.eq.0).or.(istep.eq.2)) then
!       istep: 0 calcul local
!       istep: 2 2eme etape non local
!       calcul de la contrainte dans les renforts 
!       increment de deformation
        deps=epsrf-epsr0
        if ((spre0.ne.spref).and.(.not.ppas)) then
!         on est en train de faire varier la precontrainte
!         on ne calcule que la variation de deformation elastique
!         pour qu elle soit compatible avec la precontrainte imposee
          sigr0=spref
!         taux de chargement maxi        
          mu0=abs(sigr0/syr)
          mu00=max(mu0,mu00)
!         contrainte finale        
          sigrf=sigr0
!         def elastique
          epse0=sigrf/your          
!         actualisation de la def plastique         
          eprpf=eprp0
!         actualisation de la def visqueuse  de Maxwell       
          eprmf=eprm0 
!         actualisation de la def visqueuse  de Kelvin       
          eprkf=eprk0          
!         fin de traitement si rho=0
        else
           if (ppas) then
!            initialisation de la variable interne si premier pas          
             sigr0=sprec
!            initialisation des autres variables internes du renfort
             eprp0=0.d0
             eprm0=0.d0
             eprk0=0.d0
             mu0=abs(sigr0/syr)
             mu00=max(mu0,mu00)        
           end if
!          la precontrainte n a pas varie donc on calcule normalement        
!          on teste si le calcul de relaxation est nécessaire ou pas
           if (rhor.gt.0.) then       
             if((dt.eq.0.).or.(tomr.eq.0.).or.(ekr.eq.0.) &
                .or.(.not.fl3d).or.(plast_seule)) then   
!                calcul elastoplastique : tir elastique
                 sigrf=sigr0+your*deps
!                test du critere de plasticite
                 if((sigrf-hplr*eprp0).gt.syr) then
                    f1=(sigrf-hplr*eprp0)-syr
!                   franchissement du critere en traction, on plastifie
                    dep=f1/(your+hplr)
                    sigrf=sigrf-your*dep
                 else if ((sigrf-hplr*eprp0).lt.(-syr)) then
                   f1=-(sigrf-hplr*eprp0)-syr
!                  franchissement du critere en compression, on plastifie
                   dep=-f1/(your+hplr) 
                   sigrf=sigrf-your*dep
                 else
!                  pas d increment plastique
                  dep=0.d0
                 end if
!                actualisation de la def plastique         
                 eprpf=eprp0+dep
!                actualisation de la def visqueuse         
                 eprmf=eprm0
                 eprkf=eprk0
!                coeff de consolidation debut de pas
!                prise en compt du taux de chargement maxi dans le temps
                 mu0=abs(sigr0/syr)
                 mu00=max(mu0,mu00)         
             else
!                calcul visco elastique : tir visco elastique
!                deformation elastique debut de pas
                 epse0=sigr0/your
!                coeff de consolidation debut de pas
!                prise en compt du taux de chargement
                 mu0=min(abs((sigr0-hplr*eprp0)/syr),1.d0)
                 mu00=max(mu0,mu00)
                 if(khir.gt.1.) then
                     mucr=0.6667d0*khir/(khir-1.d0)
                     km=mucr/(mucr-mu0)
                     if(km.lt.0.) then
!                          print*,'km<0 ds renfort3d',km
!                          print*,mu0,sigr0,syr
                          errr=.true.
                          go to 999
                     end if
                 else
                     km=1.d0
                 end if
!                prise en compte de la temperature
!                energie d activation fonction du chargement
                 ATR1=ATRref*exp(gamr*max(mu00,xmuthr))
                 Teta=0.5d0*(teta1+teta2)+273.15d0
                 Tetaref=ttaref+273.15d0
!                coeff  d amplification de la def differee         
!                kt=exp(-(ATR1/8.31d0)*(1.d0/Teta-1.d0/Tetaref)) 
                 if(Teta.gt.Tetaref) then 
                     kt=exp(ATR1*((Teta-Tetaref)**xnr))
                 else
                     if(Teta.lt.Tetaref) then
                        kt=exp(-ATR1*((Tetaref-Teta)**xnr))
                     else
                        kt=1.d0
                     end if
                 end if
                 if(kt.gt.200.) then
!                       print*,'Amplification thermique renfort3d :',kt
!                       print*,'Valeur tres grande'
!                       print*,'verifier les donnees ATR et XNR'
!                       print*,ATRref,Teta,xnr,Tetaref,mu00
                    errr=.true.
                    go to 999
                end if
!               potentiel de deformation differee
                epsmk=ekr*km*kt
!               deformation elastique caracteristique
                epsek=skr/your
!               coefficient de deformation differee
                kf=epsmk/epsek        
!               coeff de consolidation
                if (abs(epse0).gt.epsemin) then
                     x1=eprm0/epse0
                else
                     x1=eprm0/sign(epsemin,epse0)
               end if
               if(x1.gt.0.d0) then
!               deformation elastique dans le sens de la def de consolidation : ok         
                cc0=(exp(x1/kf))/kf
               else
!               les deformations epse et epsm sont de signe opposée: 
!               il faut attendre que la def visqueuse revienne du bon signe
!               en bloquant la consolidation à sa valeur minimale         
                cc0=1.d0/kf
               end if
!              resolution du tir visco elastique
!              depsm=dt*(epse0+theta*deps)/(tomr*cc0*(1.d0+theta*dt/(tomr*cc0)))
               taum=tomr*cc0
               tauk=tokr/kt
               K=kr/kt
!              increment de deformation totale deps
!              *** tir : calcul des increments sans deformation plastique ****             
!              preparation de la matrice de couplage 
               do i=1,3
                    do j=1,3           
                        if(i.eq.1) then
                            ann(i,j)=1.d0
                        else
                            ann(i,j)=0.d0
                        end if
                    end do
                    bn(i)=0.d0
               end do
!              1ere ligne deformation totale         
               bn(1)=deps
!              2eme ligne Kelvin
               dtk=min(dt,tauk)
               if(tauk.gt.0.) then
                  xkk=1.d0-exp(-(dt/tauk))
               else
                 xkk=0.d0
                 dtk=0.d0
               end if
               ake=-(xkk*theta)/kr
               akk=1.d0
               bk=xkk*(epse0/kr-eprk0)
!              affectation dans la matrice de couplage         
               ann(2,1)=ake
               ann(2,2)=akk
               bn(2)=bk
!              3eme ligne Maxwell
               xmm=kf*log(1.d0+dt/(kf*taum))
               ame=-xmm*theta
               amm=1.d0
               bm=epse0*xmm
!             affectation dans la matrice de couplage
              ann(3,1)=ame
              ann(3,3)=amm
              bn(3)=bm         
!             resolution du systeme lineaire                   
              call gauss3d(3,ann,xn,bn,ngf,errgauss,ipzero)
              if(errgauss.eq.1) then
               errr=.true.
!               print*,'Pb avec gauss3d tir viscoelastique ds renfort3d'
!               print*,tauk,taum,kt,kf,cc0,Teta,Tetaref,xnr,gamr,mu00
!               print*,'renfort3d theta',theta,'ttaref',ttaref
!               print*,'hplr,tomr,ekr,skr,ATRrf,khir,gamr,sprec,teta2,dt' 
!               print*, hplr,tomr,ekr,skr,ATRref,khir,gamr,sprec,teta2,dt
!               print*,'theta,eprm0,ttaref,rhor'
!               print*, theta,eprm0,ttaref,rhor
!               print*,'eprk0,eprkf,tokr,K'
!               print*, eprk0,eprkf,tokr,K
               go to 999
              end if          
!             recuperation des solutions
              depse1 = xn(1)         
              depsk1 = xn(2) 
              depsm1 = xn(3)         
!             actualisation des variables internes
              eprmf=eprm0+depsm1
              eprkf=eprk0+depsk1         
              eprpf=eprp0
              dsigr=your*depse1
              sigrf=sigr0+dsigr        
!             test de plasticité
              sigeq=sigrf-hplr*eprpf
              if(sigeq.gt.0.d0) then
                  f1=sigeq-syr
                  Ep=your
                  Hp=-hplr
              else
                  f1=-sigeq-syr
                  Ep=-your
                  Hp=hplr
              end if
              if(f1.gt.0.d0) then        
!              ecoulement plastique couple a la relaxation
!              preparation de la matrice de couplage 
               do i=1,4
                do j=1,4           
                    if(i.eq.1) then
                        ann(i,j)=1.d0
                    else
                        ann(i,j)=0.d0
                    end if
                end do
                bn(i)=0.d0
               end do
!              1ere ligne deformation totale nulle pour le retour radial        
               bn(1)=0.d0
!              2eme ligne Kelvin
!              affectation dans la matrice de couplage         
               ann(2,1)=ake
               ann(2,2)=akk
!              3eme ligne Maxwell
               ame=-xmm*theta
               amm=1.d0
!              affectation dans la matrice de couplage
               ann(3,1)=ame
               ann(3,3)=amm 
!              4eme ligne plasticite
               ann(4,1)=Ep
               ann(4,4)=Hp
               bn(4)=-f1         
!              resolution du systeme lineaire                   
               call gauss3d(4,ann,xn,bn,ngf,errgauss,ipzero) 
               if(errgauss.eq.1) then
                 errr=.true.
!                 print*,'Pb avec gauss3d viscoplasticite ds renfort3d'
                 go to 999
               end if 
!              recuperation des solutions
               depse2=xn(1)
               depsk2=xn(2)
               depsm2=xn(3)
               depsp =xn(4)
!              actialisation des variables internes 
               eprkf=eprkf+depsk2    
               eprmf=eprmf+depsm2
               eprpf=eprpf+depsp
               dsigr=your*depse2
               sigrf=sigrf+dsigr    
!              verif criteres de plasticite
               sigeq=sigrf-hplr*eprpf
               if(sigeq.gt.0.d0) then
                if(sigeq.gt.syr*tol_syr) then
!                    print*,'Pd de plasticite ds renfort3d'
!                    print*,sigeq,'>',syr
                    errr=.true.
                end if
               else
                if(sigeq.lt.-syr*tol_syr) then
!                    print*,'Pd de plasticite ds renfort3d'
!                    print*,sigeq,'<',-syr
                    errr=.true.
                end if             
               end if
               if(errr) then
!                do i=1,4
!                  do j=1,4
!                    print*,'ann(',i,j,')=',ann(i,j)
!                  end do
!                  print*,'bn(',i,')=',bn(i)
!                end do
                go to 999
               end if
             end if   
            end if 
           else
!           initialisation de la contrainte à la precontrainte pour CI
!           si rhor = 0 (au cas ou il devient non nul par la suite)
            sigr0=sprec
!           taux de chargement maxi        
            mu00=abs(sigr0/syr)
!           contrainte finale        
            sigrf=sigr0
!           actualisation de la def plastique         
            eprpf=eprp0
!           actualisation de la def visqueuse  de Maxwell       
            eprmf=eprm0 
!           actualisation de la def visqueuse  de Kelvin       
            eprkf=eprk0          
!           fin de traitement si rho=0     
           end if       
        end if                     
999    continue         
end subroutine
