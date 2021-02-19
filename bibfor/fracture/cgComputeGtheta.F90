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
!
! person_in_charge: matthieu-m.le-cren at edf.fr
!
subroutine cgComputeGtheta(cgField, cgTheta, cgStudy, cgTable)
!
use calcG_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cescre.h"
#include "asterfort/cescel.h"
#include "asterfort/cesexi.h"
#include "asterfort/chpchd.h"
#include "asterfort/chpver.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcharg.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbelem.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/mesomm.h"
#include "asterfort/rsexch.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "jeveux.h"

    type(CalcG_field), intent(in) :: cgField
    type(CalcG_theta), intent(inout) :: cgTheta
    type(CalcG_Study), intent(inout) :: cgStudy
    type(CalcG_Table), intent(inout) :: cgTable
! --------------------------------------------------------------------------------------------------
!
!     CALC_G --- Utilities
!
!    Compute G(Theta) in 2D and 3D
!
!----------------------------------------------
    integer :: iret, nsig, ino1, ino2, inga, ibid, i
    integer :: nchin, nbgrel, nel, nute, nncp
    integer :: jcesd, jcesl, jliel, ibasf
    integer :: igr, iel, ima, iad
    real(kind=8) :: gth(7)
    character(len=2)  :: codret
    character(len=8)  :: k8b, lpain(50), lpaout(1), typmo, ma
    character(len=16) :: opti, nomte
    character(len=19) :: chrota, chpesa, cf2d3d, chpres, chvolu, cf1d2d, chepsi
    character(len=19) :: chvarc, chvref, chsdeg, chdeg, ligrmo, chslag, chlag
    character(len=24) :: chsigi, celmod, sigelno, chtime, chpuls
    character(len=24) :: chgeom, type, chsig, chepsp, chvari
    character(len=24) :: pavolu, papres, pa2d3d, pepsin, pa1d2d
    character(len=24) :: lchin(50), lchout(1)
    aster_logical     :: lfonc
    integer, pointer  :: cesv(:) => null()
    real(kind=8), pointer  :: absc(:) => null()
!    integer, pointer  :: absc(:) => null()
    real(kind=8), pointer :: v_base(:) => null()
!----------------------------------------------
!
    call jemarq()
!
!   Initialisation des champs et des paramètres
    chvarc  = '&&cgtheta.VARC'
    chvref  = '&&cgtheta.VARC.REF'
    chsigi  = '&&cgtheta.CHSIGI'
    celmod  = '&&cgtheta.CELMOD'
    sigelno = '&&cgtheta.SIGELNO'
    chvolu  = '&&cgtheta.VOLU'
    cf1d2d  = '&&cgtheta.1D2D'
    cf2d3d  = '&&cgtheta.2D3D'
    chpres  = '&&cgtheta.PRES'
    chepsi  = '&&cgtheta.EPSI'
    chpesa  = '&&cgtheta.PESA'
    chrota  = '&&cgtheta.ROTA'
    chtime  = '&&cgtheta.CH_INST_R'
    chpuls  = '&&cgtheta.PULS'
    chsdeg  = '&&cgtheta.CHSDEG'
    chslag  = '&&cgtheta.CHSLAG'
    chdeg   = '&&cgtheta.CHDEG'
    chlag   = '&&cgtheta.CHLAG'
!
    k8b = '        '

!   Recuperation du champ geometrique
    call megeom(cgStudy%model, chgeom)
!
!   Recuperation du LIGREL
    ligrmo = cgStudy%model//'.MODELE'
!
!   Recuperation du comportement
    if (cgField%l_incr) then
!
        call dismoi('TYPE_RESU', cgField%result_in, 'RESULTAT', repk=type)
!
        if (type .ne. 'EVOL_NOLI') then
            call utmess('F', 'RUPTURE1_15')
        endif
!
        call rsexch('F', cgField%result_in, 'SIEF_ELGA', cgStudy%nume_ordre, chsig, iret)
        call rsexch('F', cgField%result_out, 'EPSP_ELNO', cgStudy%nume_ordre, chepsp, iret)
        call rsexch('F', cgField%result_out, 'VARI_ELNO', cgStudy%nume_ordre, chvari, iret)
    endif
!
!   Recuperation de l'etat initial
    if (cgField%l_incr) then
!
        call getvid('ETAT_INIT', 'SIGM', iocc=1, scal=chsigi, nbret=nsig)
!
!       Verification du type de champ + transfo, si necessaire en champ elno
        if (nsig .ne. 0) then
!
!           chpver renvoit 0 si OK et 1 si PB
            call chpver('C', chsigi(1:19), 'ELNO', 'SIEF_R', ino1)
            call chpver('C', chsigi(1:19), 'NOEU', 'SIEF_R', ino2)
            call chpver('C', chsigi(1:19), 'ELGA', 'SIEF_R', inga)

!           Verification du type de champ
            if (ino1.eq.1 .and. ino2.eq.1 .and. inga.eq.1) then
                call utmess('F', 'RUPTURE1_12')
            endif

!           Transformation si champ ELGA
            if (inga.eq.0) then

!               Traitement du champ pour les elements finis classiques
                call detrsd('CHAMP', celmod)
                call alchml(ligrmo, 'CALCH_G', 'PSIGINR', 'V', celmod, iret, ' ')
                call chpchd(chsigi(1:19), 'ELNO', celmod, 'OUI', 'V', sigelno)
                call chpver('F', sigelno(1:19), 'ELNO', 'SIEF_R', ibid)
            endif
        endif
    else
        nsig = 0
    endif
!
!   Recuperation des champs de temperature (T,TREF)
    call vrcins(cgStudy%model, cgStudy%material, k8b, cgStudy%time, chvarc, codret)
    call vrcref(cgStudy%model, cgStudy%material(1:8), k8b, chvref(1:19))
!
!   Traitement des charges
    call gcharg(cgStudy%model, cgStudy%loading, chvolu, cf1d2d, cf2d3d, chpres,&
                chepsi, chpesa, chrota, lfonc, cgStudy%time, cgStudy%nume_ordre)
!
!   Select name of option
    if (lfonc) then
        if (cgField%ndim .eq. 2) then
            pavolu = 'PFFVOLU'
            pa1d2d = 'PFF1D2D'
            papres = 'PPRESSF'
            pepsin = 'PEPSINF'
        else
            pavolu = 'PFFVOLU'
            pa2d3d = 'PFF2D3D'
            papres = 'PPRESSF'
            pepsin = 'PEPSINF'
        endif
!
        if (cgStudy%option .eq. 'K') then
            opti   = 'CALCH_K_G_F'
        else if (cgStudy%option .eq. 'G'.or.cgStudy%option .eq. 'G_EPSI') then
            opti   = 'CALCH_G_F'
        else
            ASSERT(ASTER_FALSE)
        endif
    else
        if (cgField%ndim .eq. 2) then
            pavolu = 'PFRVOLU'
            pa1d2d = 'PFR1D2D'
            papres = 'PPRESSR'
            pepsin = 'PEPSINR'
        else
            pavolu = 'PFRVOLU'
            pa2d3d = 'PFR2D3D'
            papres = 'PPRESSR'
            pepsin = 'PEPSINR'
        endif
!
        if (cgStudy%option .eq. 'K') then
            opti   = 'CALCH_K_G'
        else if (cgStudy%option .eq. 'G'.or.cgStudy%option .eq. 'G_EPSI') then
            opti   = 'CALCH_G'
        else
            ASSERT(ASTER_FALSE)
        endif
    endif
!
!   Création d'un champ constant par élément : LEGENDRE, LAGRANGE
    call dismoi('NOM_MAILLA', cgStudy%model, 'MODELE', repk=ma)
!
!   Récupération du nombre de groupes d'éléments du ligrel
    call jelira(ligrmo//'.LIEL', 'NMAXOC', nbgrel)
!
    if (cgtheta%discretization .eq. 'LEGENDRE') then
!
!       Champ scalaire constant par élement :
!            --  degré du polynome
        call cescre('V', chsdeg, 'ELEM', ma, 'NEUT_I',&
                    0, ' ', [-1], [-1], [-1])
!
        call jeveuo(chsdeg//'.CESD', 'L', jcesd)
        call jeveuo(chsdeg//'.CESL', 'E', jcesl)
        call jeveuo(chsdeg//'.CESV', 'E', vi=cesv)
!
    elseif (cgtheta%discretization .eq.'LINEAIRE') then
!
        call jeveuo(cgTheta%absfond, 'L', ibasf)
!
!       Champ scalaire constant par élement :
!            -- abscisse curviligne s(i-1), s(i), s(i+1)
!               pour la fonction de forme courante   
        call cescre('V', chslag, 'ELEM', ma, 'NEUT_R',&
                    0, '', [-1], [-1], [-3])
!
        call jeveuo(chslag//'.CESD', 'L', jcesd)
        call jeveuo(chslag//'.CESL', 'E', jcesl)
        call jeveuo(chslag//'.CESV', 'E', vr=absc)
!
    endif
!
!************************************************
!   Boucle sur le nombre de champ theta
!    --- 1 si 2D
!    --- ndeg ou nnof si 3D
!************************************************
    do i = 1, cgTheta%nb_theta_field
!
!       Declaration des entrees/sorties pour l'appel à calcul
        lpaout(1) = 'PGTHETA'
        lchout(1) = '&&cgtheta.CH_G'
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PDEPLAR'
        lchin(2) = cgStudy%depl
        lpain(3) = 'PTHETAR'
        lchin(3) = cgTheta%theta_factors
        lpain(4) = 'PMATERC'
        lchin(4) = cgStudy%mateco
        lpain(5) = 'PVARCPR'
        lchin(5) = chvarc
        lpain(6) = 'PVARCRR'
        lchin(6) = chvref
        lpain(7) = pavolu(1:8)
        lchin(7) = chvolu
        if (cgField%ndim .eq. 2) then
            lpain(8) = pa1d2d(1:8)
            lchin(8) = cf1d2d
        else
            lpain(8) = pa2d3d(1:8)
            lchin(8) = cf2d3d
        endif
        lpain(9) = papres(1:8)
        lchin(9) = chpres
        lpain(10) = 'PPESANR'
        lchin(10) = chpesa
        lpain(11) = 'PROTATR'
        lchin(11) = chrota
        lpain(12) = pepsin(1:8)
        lchin(12) = chepsi
        lpain(13) = 'PCOMPOR'
        lchin(13) = cgField%compor
        lpain(14) = 'PBASLOR'
        lchin(14) = cgTheta%crack//'.BASLOC'
!
        nchin = 14
!
!       Champ constant pour construire theta_i dans le te
        if (cgField%ndim .eq. 3) then
!
            do igr = 1, nbgrel
!
!               Récupération nombre d'éléments dans le groupe
                nel = nbelem(ligrmo,igr)
                call jeveuo(jexnum(ligrmo//'.LIEL', igr), 'L', jliel)
!
!               Récupération du nom du type d'éléments du groupe
                nute = typele(ligrmo,igr)
                call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
!
!               Boucle sur les éléments du groupe
                do iel = 1, nel
!
                    ima=zi(jliel-1+iel)
                    if (ima .lt. 0) cycle
!
                    if (cgTheta%discretization.eq.'LEGENDRE') then       
!            
                        call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
                        iad = abs(iad)
                        zl(jcesl-1+iad)=.true.
                        cesv(iad)= i - 1
!
                    elseif (cgtheta%discretization .eq.'LINEAIRE') then
!
                        call cesexi('C', jcesd, jcesl, ima, 1, 1, 1, iad)
                        iad = abs(iad)
                        zl(jcesl-1+iad)=.true.
                        if (i .eq. 1 ) then 
                            absc(iad)= zr(ibasf-1+i)
                        else 
                            absc(iad)= zr(ibasf-1+i-1)
                        endif
!
                        call cesexi('C', jcesd, jcesl, ima, 1, 1, 2, iad)
                        iad = abs(iad)
                        zl(jcesl-1+iad)=.true.
                        absc(iad)= zr(ibasf-1+i)
!
                        call cesexi('C', jcesd, jcesl, ima, 1, 1, 3, iad)
                        iad = abs(iad)
                        zl(jcesl-1+iad)=.true.
                        if (i .eq. cgTheta%nb_theta_field ) then     
                            absc(iad)= zr(ibasf-1+i)
                        else
                            absc(iad)= zr(ibasf-1+i+1)
                        endif
!
                    endif
                enddo
            enddo
!
            if (cgTheta%discretization.eq.'LEGENDRE') then                  
!               Conversion du champ par élément simple en champ par éléments
                call cescel(chsdeg, ligrmo, 'CALCH_G', 'PDEG', 'NON',&
                            nncp, 'V', chdeg, 'F', iret)
!
                lpain(nchin+1) = 'PDEG'
                lchin(nchin+1) = chdeg
!
            elseif (cgtheta%discretization .eq.'LINEAIRE') then
!               Conversion du champ par élément simple en champ par éléments
                call cescel(chslag, ligrmo, 'CALCH_G', 'PLAG', 'NON',&
                            nncp, 'V', chlag, 'F', iret)
!
                lpain(nchin+1) = 'PLAG'
                lchin(nchin+1) = chlag
!
            endif
!         
            nchin = nchin + 1
!
        endif

        if (cgStudy%option .eq. 'K') then
            lpain(nchin+1) = 'PCOURB'
            lchin(nchin+1) = cgTheta%curvature
            lpain(nchin+2) = 'PLSN'
            lchin(nchin+2) = cgTheta%crack//'.LNNO'
            lpain(nchin+3) = 'PLST'
            lchin(nchin+3) = cgTheta%crack//'.LTNO'
            nchin = nchin + 3

        endif
!
        if (opti .eq. 'CALCH_G_F' .or. opti .eq. 'CALCH_K_G_F') then
            call mecact('V', chtime, 'MODELE', ligrmo, 'INST_R',&
                        ncmp=1, nomcmp='INST', sr=cgStudy%time)
            lpain(nchin+1) = 'PTEMPSR'
            lchin(nchin+1) = chtime
            nchin = nchin + 1
        endif
!
        if (cgStudy%l_modal) then
            call mecact('V', chpuls, 'MODELE', ligrmo, 'FREQ_R  ',&
                        ncmp=1, nomcmp='FREQ   ', sr=cgStudy%pulse)
            nchin = nchin + 1
            lpain(nchin) = 'PPULPRO'
            lchin(nchin) = chpuls
        endif
!
        if (cgField%l_incr ) then
            if (opti .eq. 'CALCH_G_F' .or. opti .eq. 'CALCH_G') then
                lpain(nchin+1) = 'PCONTRR'
                lchin(nchin+1) = chsig
                lpain(nchin+2) = 'PDEFOPL'
                lchin(nchin+2) = chepsp
                lpain(nchin+3) = 'PVARIPR'
                lchin(nchin+3) = chvari
                nchin = nchin + 3
            endif
        endif
!
!       Champ de contrainte initiale
        if (nsig .ne. 0) then
          if (inga .eq. 0) then
!           Champ de contrainte initiale transforme en ELNO
            lpain(nchin+1) = 'PSIGINR'
            lchin(nchin+1) = sigelno
            nchin = nchin + 1
          else
!           Champ de contrainte initiale donne par l'utilisateur (NOEUD ou ELNO)
            lpain(nchin+1) = 'PSIGINR'
            lchin(nchin+1) = chsigi
            nchin = nchin + 1
          endif
        endif
!
        if (cgStudy%option .eq. 'G'.or.cgStudy%option .eq. 'G_EPSI') then
            if (cgStudy%vitesse .ne. ' ') then
                lpain(nchin+1) = 'PVITESS'
                lchin(nchin+1) = cgStudy%vitesse
                lpain(nchin+2) = 'PACCELE'
                lchin(nchin+2) = cgStudy%acce
                nchin = nchin + 2
            endif
        endif
!
!       Recuperationdes contraintes du resultat pour option G
        if (cgStudy%option .eq. 'G') then
            call rsexch(' ', cgField%result_in, 'SIEF_ELGA', cgStudy%nume_ordre, chsig, iret)
!
            if (iret .ne. 0) then
!               Probleme pour recuperer SIEF_ELGA avec ce numero d'ordre
                call utmess('F', 'RUPTURE0_94', si=cgStudy%nume_ordre)
            endif
!
            lpain(nchin+1) = 'PCONTGR'
            lchin(nchin+1) = chsig
            nchin = nchin + 1
        endif
!
!       Sommation des G élémentaires
        call calcul('S', opti, ligrmo, nchin, lchin,&
                    lpain, 1, lchout, lpaout, 'V', 'OUI')
!
!       G, K1, K2, K3, FIC1, FIC2, FIC3 en 2D
        call mesomm(lchout(1), 7, vr=gth)
!
        print* , ' G = ' , gth(1) , 'K = ' , i
!
    end do
!
!    Cas axis, on normalise par 1/R
     call dismoi('MODELISATION', cgStudy%model, 'MODELE', repk=typmo)
     if (typmo(1:4) .eq. 'AXIS') then
         do i = 1, 7
             call cgTheta%getBaseLoc(v_base)             
             gth(i) = gth(i)/v_base(1)
         end do
     endif
!
    cgStudy%gth = 0.d0
    if (cgTheta%symech .eq. 'OUI') then
        cgStudy%gth(1:7) = [ 2.d0*gth(1), 2.d0*gth(2), 0.d0, 0.d0,&
                                          2.d0*gth(5), 0.d0, 0.d0]
    else
        cgStudy%gth(1:7) = gth(1:7)
    endif
!
!   Ajout des valeurs dans la table de G
    call cgTable%addValues(cgField, cgStudy)
!
    call detrsd('CHAMP_GD', chvarc)
    call detrsd('CHAMP_GD', chvref)
    call detrsd('CHAMP_GD', cf1d2d)
    call detrsd('CHAMP_GD', cf2d3d)
    call detrsd('CHAMP_GD', chepsi)
    call detrsd('CHAMP_GD', chpesa)
    call detrsd('CHAMP_GD', chpres)
    call detrsd('CHAMP_GD', chrota)
    call detrsd('CHAMP_GD', chtime)
    call detrsd('CHAMP_GD', chvolu)
!
    call jedema()
end subroutine
