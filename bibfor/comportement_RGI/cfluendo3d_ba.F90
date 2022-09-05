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

subroutine cfluendo3d_ba(fami, kpg, ksp, ndim, imate,&
                  compor, carcri, instam, instap, epsm,&
                  deps, sigm, vim, option,&
                  sigp, vip, typmod,&
                  dsidep, codret)
! person_in_charge: etienne.grimal@edf.fr
!=====================================================================
!      Ficher de base de FLUAGE ET ENDOMMAGEMENT
!=====================================================================
    implicit none
#include "asterf_types.h"
#include "asterfort/rcvarc.h"
#include "asterfort/rcvalb.h"
#include "asterfort/fluendo3d_ba.h"
#include "asterc/r8prem.h"
#include "asterfort/utmess.h"
!
!
    integer :: imate, ndim, kpg, ksp, codret, iret
    integer :: ifour,ifour11,kerre,nmatflu
    real(kind=8) :: instam, instap
    real(kind=8) :: epsm(6)
    real(kind=8) :: epsmc(6), depsc(6), epsc(6)
    real(kind=8) :: sigm(6), sigp(6)
    real(kind=8) :: vim(*), vip(*), tm, tp, tref
    real(kind=8) :: dsidep(6, 6)
    character(len=16) :: compor(*), option
    character(len=8) :: typmod(*)
    character(len=*) :: fami
    real(kind=8), intent(in) :: carcri(*)
!
! DECLARATIONS LOCALES
    integer :: nvarbe, nxmat,nvarac
!   Nombre de paramtre relatif au beton uniquement
    parameter (nvarbe=59)
!   Nombre de paramtre relatif a l'acier uniquement(1+21*5)
    parameter (nvarac=106)
!   Nombre para total
    parameter (nmatflu=nvarbe+nvarac)
!   Nombre parametre pour xmat = 4+ nmatflu
    parameter (nxmat=4+nmatflu)


    character(len=8) :: nomres(nmatflu),nomres1(4)
    real(kind=8) :: valres(nmatflu), xmat(nxmat), rbid,valres1(4)
    integer :: nvari, nstrs, mfr, i, j
    integer :: retour(nmatflu),retour1(4),iteflumax
    real(kind=8) :: d(6, 6), e, nu, coef, coef1, coef2, coef3
    real(kind=8) :: zero, un, deux, rac2,var0(6),sig0(6)
!
    real(kind=8) :: hydrm, hydrp, sechp, sechm, sref, vgm, vgp
    real(kind=8) :: alpham, alphap, deps(6),epstf3d(6),teta13d,teta23d,sech


!   variables de transfert de donnees ( a declarer suivant idvar4 et idvisc)
      integer nmat3d,nstrs3d,nvari3d,ierr1,mfr11
      integer NMATAILX,NVARFLU,NMATRAG,NVARRAG

    integer :: nbelas3d
      parameter (nbelas3d=4)
      parameter (NMATAILX=0)
!     Nombre des variable totale du modele
      parameter (NVARFLU=163)

!   taille du pseudo vecteur des contraintes pour fluendo3d (tjrs 6
!   en raison de son utilisation dans fludes3d qui la suppose à 6)
      parameter (nstrs3d=6)

!   mettre a jour ici le nbre de parametres materiaux et variables
!   interne de l option du modele
      parameter (NMATRAG=0)
      parameter (NVARRAG=0)


!   nombre totale de parametres et variables internes option comprise
      parameter (nmat3d=nbelas3d+NMATFLU+NMATRAG+NMATAILX)
      parameter (nvari3d=NVARFLU+NVARRAG)

!   nbre de parametres materiaux sans les tailles des elements
      real(kind=8) :: xmat3d(nmat3d),sig03d(6),sigf3d(6),depst3d(6)
      real(kind=8) :: var03d(nvari3d),varf3d(nvari3d),varf(nvari3d),sigf(6)
!   indicateur d isotropie initiale
      aster_logical :: iso1,local11,end3d,fl3d
!   temperatures debut et fin de pas , moyenne, pas de temps, volule rgi
      real(kind=8) :: dt3d,phig3d
     
    parameter       (nvari=NVARFLU)
    integer, dimension(2) :: vali


    sig03d(:) = 0.d0
    sigf3d(:) = 0.d0
    var03d(:) = 0.d0
    varf3d(:) = 0.d0
    varf(:)   = 0.d0
    depsc(:)  = 0.d0
    epsmc(:)  = 0.d0
    xmat(:)   = 0.d0
    var0(:)   = 0.d0
    sig0(:)   = 0.d0
    sigf(:)   = 0.d0
    depst3d(:) = 0.d0
    valres1(:)=0.d0
    valres(:)=0.d0

    iteflumax = int(carcri(1))

!
! APPEL DE RCVARC POUR LE CALCUL DE LA TEMPERATURE
!
    call rcvarc('f', 'TEMP', '-', fami, kpg,&
                ksp, tm, iret)
    call rcvarc('f', 'TEMP', '+', fami, kpg,&
                ksp, tp, iret)
    call rcvarc('f', 'TEMP', 'REF', fami, kpg,&
                ksp, tref, iret)
!
! ------------------------------------------------
!!     RECUPERATION DE L HYDRATATION DEBUT DE PAS
    call rcvarc(' ', 'HYDR', '-', fami, kpg,&
                ksp, hydrm, codret)
    if (codret .ne. 0) then
        hydrm=0.d0
        codret = 0
    endif
!!
!! ------------------------------------------------
!!     RECUPERATION DE L HYDRATATION FIN DE PAS
    call rcvarc(' ', 'HYDR', '+', fami, kpg,&
                ksp, hydrp, codret)
    if (codret .ne. 0) then
        hydrp=0.d0
        codret = 0
    endif
!
! ------------------------------------------------
!     RECUPERATION DU SECHAGE
    call rcvarc(' ', 'SECH', '+', fami, kpg,&
                ksp, sechp, iret)
    if (iret .ne. 0) sechp=0.d0
!    call rcvarc(' ', 'SECH', '-', fami, kpg,&
!                ksp, sechm, iret)
!    if (iret .ne. 0) sechm=0.d0
    sechm=0.d0
    call rcvarc(' ', 'SECH', 'REF', fami, kpg,&
                ksp, sref, iret)
    if (iret .ne. 0) sref=0.d0

    sech = sechp
!
!   le séchage de référence doit être nul
!
    if (sref .ne. 0.d0) call utmess('F', 'COMPOR3_9', sk=compor(1))
!
! -----------------------------------------------
!     RECUPERATION DU VOLUME DE GEL DEBUT DE PAS
    call rcvarc(' ', 'X1', '-', fami, kpg,&
                ksp, vgm, codret)
    if (codret .ne. 0) then
        vgm=0.d0
        codret = 0
    endif
!
! ------------------------------------------------
!     RECUPERATION DU VOLUME DE GEL FIN DE PAS
    call rcvarc(' ', 'X1', '+', fami, kpg,&
                ksp, vgp, codret)
    if (codret .ne. 0) then
        vgp=0.d0
        codret = 0
    endif
!
! ------------------------------------------------
!
!   vérification que B_ENDOGE et K_DESSIC sont nuls
    nomres1(1)='B_ENDOGE'
    nomres1(2)='K_DESSIC'

    call rcvalb(fami, kpg, ksp, '+', imate,&
                ' ', 'ELAS', 0, ' ', [0.d0],&
                2, nomres1, valres1, retour1, 0)
    do i=1,2
        if (retour1(i) .eq. 0)then
            if (valres1(i) .ne. 0.d0)then
                call utmess('F', 'COMPOR3_40', nk=2, valk=[compor(1), nomres1(i)])
            endif
        endif
    enddo

    nomres1(1)='E'
    nomres1(2)='NU'
    nomres1(3)='RHO'
    nomres1(4)='ALPHA'
!
    call rcvalb(fami, kpg, ksp, '-', imate,&
                ' ', 'ELAS', 0, ' ', [0.d0],&
                4, nomres1, valres1, retour1, 2)
!
!
!        MODULES INSTANTANES ISOTROPES
    xmat(1) = valres1(1)
    xmat(2) = valres1(2)
    xmat(3)  = valres1(3)
    xmat(4) = valres1(4)
    alpham = valres1(4)
!
! --- EVALUATION PARAMETERES MATERIAU ELASTIQUES A T+
    call rcvalb(fami, kpg, ksp, '+', imate,&
                ' ', 'ELAS', 0, ' ', [0.d0],&
                4, nomres1, valres1, retour1, 2)

    alphap = valres1(4)


!
! ------------------------------------------------------------------
! --  RETRAIT INCREMENT DE DEFORMATION DUE A LA DILATATION THERMIQUE
! ------------------------------------------------------------------
    if (ndim .eq. 2) then
        nstrs = 4
    else
        nstrs = 6
    endif

    if ((option(1:9).eq.'RAPH_MECA') .or. (option(1:9).eq.'FULL_MECA')) then
        do i = 1, 3
            depsc(i) = deps(i) - (alphap*(tp-tref)-alpham*(tm-tref))
            epsmc(i) = epsm(i) - alpham*(tm-tref)
            epsc(i) = epsmc(i) + depsc(i)
        end do
        do i = 4, nstrs
            depsc(i) = deps(i)
            epsc(i) = epsm(i) + depsc(i)
        end do
    endif
!
! ------------------------------------------------------------------
! --  RECUPERATION PARAMETRES MATERIAU DU MODELE FLUENDO3D
! ------------------------------------------------------------------
!
    nomres(1) = 'HYDR'
    nomres(2) = 'HYDS'
    nomres(3) = 'RT'
    nomres(4) = 'REF'
    nomres(5) = 'RC'
    nomres(6) = 'DELT'
    nomres(7) = 'BETA'
    nomres(8) = 'EPT'
    nomres(9) = 'HRGI'
    nomres(10)= 'VVRG'
    nomres(11)= 'KGEL'
    nomres(12)= 'GFT'
    nomres(13)= 'EKFL'
    nomres(14)= 'YKSY'
    nomres(15)= 'XFLU'
    nomres(16)= 'TAUK'
    nomres(17)= 'TAUM'
    nomres(18)= 'NRJM'
    nomres(19)= 'DT80'
    nomres(20)= 'BSHR'
    nomres(21)= 'MSHR'
    nomres(22)= 'PORO'
    nomres(23)= 'VRAG'
    nomres(24)= 'NRJF'
    nomres(25)= 'SFLD'
    nomres(26)= 'MVGN'
    nomres(27)= 'EPC'
    nomres(28)= 'EKDC'
    nomres(29)= 'EKRG'
    nomres(30)= 'GFR'
    nomres(31)= 'ALAT'
    nomres(32)= 'KRGI'
    nomres(33)= 'TREF'
    nomres(34)= 'TSTH'
    nomres(35)= 'DFMX'
    nomres(36)= 'TAUG'
    nomres(37)= 'NRJG'
    nomres(38)= 'SRSG'
    nomres(39)= 'TRAG'
    nomres(40)= 'DIM3'
    nomres(41)= 'TDEF'
    nomres(42)= 'NRJP'
    nomres(43)= 'SRSD'
    nomres(44)= 'VDEF'
    nomres(45)= 'CNAD'
    nomres(46)= 'SSAD'
    nomres(47)= 'CNAK'
    nomres(48)= 'CNAB'
    nomres(49)= 'EXND'
    nomres(50)= 'EXMD'
    nomres(51)= 'TTDD'
    nomres(52)= 'TDID'
    nomres(53)= 'TFID'
    nomres(54)= 'NRJD'
    nomres(55)= 'TTRD'
    nomres(56)= 'TTKF'
    nomres(57)= 'HPEV'

    nomres(58)= 'YOUM'
    nomres(59)= 'NUM'

    nomres(nvarbe+1)= 'NREN'

    nomres(nvarbe+2)= 'ROA1'
    nomres(nvarbe+3)= 'DEQ1'
    nomres(nvarbe+4)= 'YOR1'
    nomres(nvarbe+5)= 'SYR1'
    nomres(nvarbe+6)= 'TYR1'
    nomres(nvarbe+7)= 'VR11'
    nomres(nvarbe+8)= 'VR12'
    nomres(nvarbe+9)= 'VR13'
    nomres(nvarbe+10)= 'HPL1'
    nomres(nvarbe+11)= 'TMR1'
    nomres(nvarbe+12)= 'EKR1'
    nomres(nvarbe+13)= 'SKR1'
    nomres(nvarbe+14)= 'ATR1'
    nomres(nvarbe+15)= 'CTM1'
    nomres(nvarbe+16)= 'XFL1'
    nomres(nvarbe+17)= 'PRE1'
    nomres(nvarbe+18)= 'TTR1'
    nomres(nvarbe+19)= 'XNR1'
    nomres(nvarbe+20)= 'MUS1'
    nomres(nvarbe+21)= 'TKR1'
    nomres(nvarbe+22)= 'YKY1'

    nomres(nvarbe+23)= 'ROA2'
    nomres(nvarbe+24)= 'DEQ2'
    nomres(nvarbe+25)= 'YOR2'
    nomres(nvarbe+26)= 'SYR2'
    nomres(nvarbe+27)= 'TYR2'
    nomres(nvarbe+28)= 'VR21'
    nomres(nvarbe+29)= 'VR22'
    nomres(nvarbe+30)= 'VR23'
    nomres(nvarbe+31)= 'HPL2'
    nomres(nvarbe+32)= 'TMR2'
    nomres(nvarbe+33)= 'EKR2'
    nomres(nvarbe+34)= 'SKR2'
    nomres(nvarbe+35)= 'ATR2'
    nomres(nvarbe+36)= 'CTM2'
    nomres(nvarbe+37)= 'XFL2'
    nomres(nvarbe+38)= 'PRE2'
    nomres(nvarbe+39)= 'TTR2'
    nomres(nvarbe+40)= 'XNR2'
    nomres(nvarbe+41)= 'MUS2'
    nomres(nvarbe+42)= 'TKR2'
    nomres(nvarbe+43)= 'YKY2'

    nomres(nvarbe+44)= 'ROA3'
    nomres(nvarbe+45)= 'DEQ3'
    nomres(nvarbe+46)= 'YOR3'
    nomres(nvarbe+47)= 'SYR3'
    nomres(nvarbe+48)= 'TYR3'
    nomres(nvarbe+49)= 'VR31'
    nomres(nvarbe+50)= 'VR32'
    nomres(nvarbe+51)= 'VR33'
    nomres(nvarbe+52)= 'HPL3'
    nomres(nvarbe+53)= 'TMR3'



    nomres(nvarbe+54)= 'EKR3'
    nomres(nvarbe+55)= 'SKR3'
    nomres(nvarbe+56)= 'ATR3'
    nomres(nvarbe+57)= 'CTM3'
    nomres(nvarbe+58)= 'XFL3'
    nomres(nvarbe+59)= 'PRE3'
    nomres(nvarbe+60)= 'TTR3'
    nomres(nvarbe+61)= 'XNR3'
    nomres(nvarbe+62)= 'MUS3'
    nomres(nvarbe+63)= 'TKR3'
    nomres(nvarbe+64)= 'YKY3'

    nomres(nvarbe+65)='ROA4'
    nomres(nvarbe+66)='DEQ4'
    nomres(nvarbe+67)='YOR4'
    nomres(nvarbe+68)='SYR4'
    nomres(nvarbe+69)='TYR4'
    nomres(nvarbe+70)='VR41'
    nomres(nvarbe+71)='VR42'
    nomres(nvarbe+72)='VR43'
    nomres(nvarbe+73)='HPL4'
    nomres(nvarbe+74)='TMR4'
    nomres(nvarbe+75)='EKR4'
    nomres(nvarbe+76)='SKR4'
    nomres(nvarbe+77)='ATR4'
    nomres(nvarbe+78)='CTM4'
    nomres(nvarbe+79)='XFL4'
    nomres(nvarbe+80)='PRE4'
    nomres(nvarbe+81)='TTR4'
    nomres(nvarbe+82)='XNR4'
    nomres(nvarbe+83)='MUS4'
    nomres(nvarbe+84)='TKR4'
    nomres(nvarbe+85)='YKY4'

    nomres(nvarbe+86)='ROA5'
    nomres(nvarbe+87)='DEQ5'
    nomres(nvarbe+88)='YOR5'
    nomres(nvarbe+89)='SYR5'
    nomres(nvarbe+90)='TYR5'
    nomres(nvarbe+91)='VR51'
    nomres(nvarbe+92)='VR52'
    nomres(nvarbe+93)='VR53'
    nomres(nvarbe+94)='HPL5'
    nomres(nvarbe+95)='TMR5'
    nomres(nvarbe+96)='EKR5'
    nomres(nvarbe+97)='SKR5'
    nomres(nvarbe+98)='ATR5'
    nomres(nvarbe+99)='CTM5'
    nomres(nvarbe+100)='XFL5'
    nomres(nvarbe+101)='PRE5'
    nomres(nvarbe+102)='TTR5'
    nomres(nvarbe+103)='XNR5'
    nomres(nvarbe+104)='MUS5'
    nomres(nvarbe+105)='TKR5'
    nomres(nvarbe+106)='YKY5'
!
    rbid = 0.d0

!
    call rcvalb(fami, kpg, ksp, '-', imate,&
                ' ', compor(1), 0, ' ', [rbid],&
                nmatflu, nomres, valres, retour, 0, nan='NON')
!
!
! --- ON REMPLIT XMAT DE nbelas+1 A nbelas+nmatflu
    do i = (nbelas3d+1), (nbelas3d+nmatflu)
        xmat(i) = valres(i-nbelas3d)
    end do
!
!-----VALEUR FIXEE PROVISOIREMENT POUR MFR
        mfr = 1
!-----------------------------------------

    if (typmod(1)(1:2) .eq. '3D') then
        ifour = 2
    else
        if (typmod(1)(1:6) .eq. 'D_PLAN') then
            ifour = -1
        else
               ifour = 0
        endif
    endif
   if (compor(1) .eq. 'FLUA_PORO_BETON') then
        fl3d = .true.
           end3d= .false.
   else
              if (compor(1) .eq. 'ENDO_PORO_BETON') then
            fl3d = .false.
                  end3d= .true.
              else
                  fl3d = .true.
                  end3d= .true.
              endif
   endif

!    pas de temps
       dt3d=instap-instam

          do i=1,nmat3d
             xmat3d(i)=xmat(i)
          end do

!    variables internes
       if(nvari3d.ne.nvari) then
            vali(1) = nvari3d
            vali(2) = nvari
            call utmess('E', 'COMPOR3_16', ni=2, vali=vali)
            kerre=1
       end if
       do i=1,nvari3d
             var03d(i)=vim(i)
             varf3d(i)=vip(i)
       end do

!   initialisation des contraintes effectives si premier pas
      if (abs(var03d(64)-1.d0).ge.r8prem()) then

         do i=1,6
!         initialisation de l etage elastique
            if(i.le.3) then
!              on retire la pression intra poreuse au cas où
!              elle aurait été initialisée dans les vari
                 var03d(18+i)=sig03d(i)-var03d(61)
            else
                 var03d(18+i)=sig03d(i)
            end if
!         idem pour l etage de Kelvin
            var03d(49+i)=var03d(18+i)
         end do
!      on met l indicateur de 1er passage a 1
         varf3d(64)=1.
       else
!     on reconduit l indicateur de 1er passage
        varf3d(64)=var03d(64)
       end if

!    autres parametres a renseigner
!    temperature moyenne sur le pas
       teta13d=tm
       teta23d=tp
!    initialisation indicateur d erreur
       ierr1=0
!    indicateur isostropie elastique et de resistance
       iso1=.true.
!    numero de la formulation (33 pour poreux)
       mfr11=mfr
!    type de formulation
       ifour11=ifour


!    controle de regularisation en cas d endommagement
       local11=.true.
       rac2 = sqrt(2.d0)
         do 70 i = 1, 3
            depst3d(i) = depsc(i)
            epstf3d(i) = epsc(i)
 70     continue
         do 71 i = 4, nstrs
           depst3d(i) = depsc(i) * rac2
           epstf3d(i) = epsc(i) * rac2
 71     continue

         do 72 i = 1, 3
           sig03d(i) = sigm(i)
 72     continue

         do 73 i = 4, nstrs
           sig03d(i) = sigm(i) / rac2
 73     continue


       phig3d=0.d0
       call fluendo3d_ba(xmat3d,sig03d,sigf3d,depst3d,&
                      nstrs3d,var03d,varf3d,nvari3d,nbelas3d,&
                      teta13d,teta23d,dt3d,phig3d,epstf3d,ierr1,&
                      iso1, mfr11,end3d,fl3d,local11,&
                      ndim,nvarbe,iteflumax,sech)

       do i=1,3
            sigp(i) = sigf3d(i)
        end do
        do 80 i = 4, nstrs
            sigp(i) = sigf3d(i) * rac2
 80     continue

        if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
            do i=1,nvari
               vip(i) = varf3d(i)
            end do
        endif


!**********************************************************************
        if(ierr1.eq.0)then
            kerre=0
        else
            kerre=1
        end if


!*********************************************************************!

    if ((option(1:9).eq.'RIGI_MECA') .or. (option(1:9).eq.'FULL_MECA')) then
!
        zero = 0.d0
        un = 1.d0
        deux = 2.d0
!
        d(:,:) = zero
!
        e = xmat(1)
        nu = xmat(2)

        coef = un/ ((un+nu)* (un-deux*nu))
        coef1 = e* (un-nu)*coef
        coef2 = e*nu*coef
        coef3 = e/ (un+nu)

        d(1,1) = coef1
        d(1,2) = coef2
        d(1,3) = coef2

        d(2,1) = coef2
        d(2,2) = coef1
        d(2,3) = coef2

        d(3,1) = coef2
        d(3,2) = coef2
        d(3,3) = coef1

        d(4,4) = 0.5d0*coef3
        d(5,5) = 0.5d0*coef3
        d(6,6) = 0.5d0*coef3

        do 90 i = 1, nstrs
            do 100 j = 1, nstrs
               dsidep(i,j) = d(i,j)
100         continue
 90     continue
!
    endif
!
end subroutine
