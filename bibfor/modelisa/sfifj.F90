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

subroutine sfifj(nomres)
    implicit none
!     CALCUL DE LA FONCTION ACCEPTANCE
!     TURBULENCE DE COUCHE LIMITE
!     AUTEUR : G. ROUSSEAU
!-----------------------------------------------------------------------
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/accep1.h"
#include "asterfort/accep2.h"
#include "asterfort/accept.h"
#include "asterfort/chpver.h"
#include "asterfort/dspprs.h"
#include "asterfort/evalis.h"
#include "asterfort/fointe.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/rsadpa.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    integer(kind=8) :: nfinit, nfin, nbm, nbpoin, nbid
    integer(kind=8) :: npoin, iff, lvale, ibid, in
    integer(kind=8) :: im1, im2, nvecx, nvecy
    integer(kind=8) :: nveco, ier, ncham, jpara
    integer(kind=8) :: lnumi, lnumj, lfreq, mxval, nbabs, ij
    real(kind=8) :: fmin, fmax, finit, ffin, df, f, prs
    real(kind=8) :: kste, uflui, dhyd, rho, jc, fcoupu, fmodel
    real(kind=8) :: dir(3, 3), fcoup, fcoup_red
    real(kind=8) :: deuxpi, puls, uc, ut, long1, long2
    real(kind=8) :: valr
    character(len=8) :: k8b, nomres, is
    character(len=8) :: spectr, method
    character(len=19) :: base, fonct, chamno, pg, phi, sphi
    character(len=19) :: long1f, long2f
    character(len=24) :: ligrmo, val_spec
    character(len=24) :: chnumi, chnumj, chfreq, chvale
    aster_logical :: yang
    real(kind=8), pointer :: vecx(:) => null()
    real(kind=8), pointer :: vecy(:) => null()
    real(kind=8), pointer :: vecz(:) => null()
    integer(kind=8), pointer :: ordr(:) => null()
    character(len=16), pointer :: vate(:) => null()
    real(kind=8), pointer :: vare(:) => null()
!
    data deuxpi/6.28318530718d0/, yang/.false./
!
!-----------------------------------------------------------------------
    call jemarq()
!
! RECHERCHE DE LA PRESENCE D'UN CHAMNO
    call getvid(' ', 'CHAM_NO', nbval=0, nbret=ncham)
!
    if (ncham .eq. 0) then
!
! RECUPERATION DE LA BASE MODALE
        call getvid(' ', 'MODE_MECA', scal=base, nbret=ibid)
!
        call jeveuo(base//'.ORDR', 'L', vi=ordr)
        call jelira(base//'.ORDR', 'LONUTI', nbm)
!
! RECUPERATION DE LA FREQUENCE MINIMALE ET MAX DES MODES
!
        call rsadpa(base, 'L', 1, 'FREQ', ordr(1), &
                    0, sjv=jpara, styp=k8b)
        fmin = zr(jpara)
        call rsadpa(base, 'L', 1, 'FREQ', ordr(nbm), &
                    0, sjv=jpara, styp=k8b)
        fmax = zr(jpara)
    else
        call getvid(' ', 'CHAM_NO', scal=chamno, nbret=ibid)
        call chpver('F', chamno, 'NOEU', 'DEPL_R', ier)
        call getvid(' ', 'CHAM_NO', nbval=0, nbret=ncham)
        nbm = -ncham
    end if
!
! RECUPERATION DE LA FREQUENCE MINIMALE ET MAX DE LA PLAGE
! DE FREQUENCE ETUDIEE
!
    call getvr8(' ', 'FREQ_INIT', scal=finit, nbret=nfinit)
    call getvr8(' ', 'FREQ_FIN', scal=ffin, nbret=nfin)
    if ((ffin-finit) .lt. r8prem()) then
        call utmess('F', 'MODELISA6_97')
    end if
!
    if (nfinit .lt. 0) then
        if (ncham .ne. 0) then
            call utmess('F', 'MODELISA6_98')
        end if
        valr = fmin
        call utmess('I', 'MODELISA9_15', sr=valr)
        finit = fmin
    end if
    if (nfin .lt. 0) then
        if (ncham .ne. 0) then
            call utmess('F', 'MODELISA6_99')
        end if
        valr = fmax
        call utmess('I', 'MODELISA9_16', sr=valr)
        ffin = fmax
    end if
!
! DISCRETISATION FREQUENTIELLE
    call getvis(' ', 'NB_POIN', scal=nbpoin, nbret=npoin)
!
! PAS FREQUENTIEL
    df = (ffin-finit)/(nbpoin-1)
    if (df .lt. r8prem()) then
        call utmess('F', 'MODELISA7_1')
    end if
!
! CALCUL DE L'ACCEPTANCE
!
    call getvid(' ', 'SPEC_TURB', scal=spectr, nbret=ibid)
    call jeveuo(spectr//'           .VARE', 'L', vr=vare)
    call jeveuo(spectr//'           .VATE', 'L', vk16=vate)

! RECUPERATION DES CONSTANTES DU SPECTRES DU
! MODELE 5 : CONSTANT PUIS NUL POUR FR > 10
!
    val_spec = vate(1)

    if (vate(1) .eq. 'SPEC_CORR_CONV_1') then
        uflui = vare(3)
        rho = vare(4)
        fcoupu = vare(5)
        kste = vare(6)
        dhyd = vare(7)
! LONGUEURS DE CORRELATION
        long1 = vare(1)
        long2 = vare(2)
! VITESSE CONVECTIVE RADIALE (METHODE AU-YANG)
        uc = vare(8)*uflui
! VITESSE CONVECTIVE ORTHORADIALE (METHODE AU-YANG)
        ut = vare(9)*uflui

! CALCUL DE LA FREQUENCE DE COUPURE PRONE PAR LE MODELE
! ET COMPARAISON AVEC LA FREQUENCE DE COUPURE DONNEE PAR
! L UTILISATEUR
!
        fmodel = 10.d0*uflui/dhyd
        fcoup_red = 0.d0
        if (fcoupu .le. fmodel) then
            valr = fcoupu
            call utmess('I', 'MODELISA9_17', sr=valr)
            valr = fmodel
            call utmess('I', 'MODELISA9_18', sr=valr)
            call utmess('I', 'MODELISA9_19')
            fcoup = fcoupu
            fcoup_red = fcoup*dhyd/uflui
        else
            valr = fcoupu
            call utmess('I', 'MODELISA9_20', sr=valr)
            valr = fmodel
            call utmess('I', 'MODELISA9_21', sr=valr)
            call utmess('I', 'MODELISA9_22')
            fcoup = fmodel
            fcoup_red = 10.d0
        end if
!
! RECUPERATION DE LA METHOD DE LA FONCTION
! DE COHERENCE
!
        method = vate(11) (1:8)

    else if (vate(1) .eq. 'SPEC_CORR_CONV_2') then
        uflui = vare(1)
        fcoup = vare(3)
        method = vate(5) (1:8)
        fonct = vate(2)
        long1f = vate(3)
        long2f = vate(4)
        uc = vare(4)*uflui
        ut = vare(2)*uflui
    else if (vate(1) .eq. 'SPEC_CORR_CONV_3') then
        fonct = vate(2)
        goto 10
    end if
!
!
! RECUPERATION DES DIRECTIONS DU PLAN DE LA PLANCHE
    if (method(1:6) .eq. 'CORCOS') then
        call getvr8(' ', 'VECT_X', nbval=0, nbret=nvecx)
        nvecx = -nvecx
        if (nvecx .gt. 0) then
            AS_ALLOCATE(vr=vecx, size=3)
            call getvr8(' ', 'VECT_X', nbval=nvecx, vect=vecx, nbret=nbid)
        end if
        call getvr8(' ', 'VECT_Y', nbval=0, nbret=nvecy)
        nvecy = -nvecy
        if (nvecy .gt. 0) then
            AS_ALLOCATE(vr=vecy, size=3)
            call getvr8(' ', 'VECT_Y', nbval=nvecy, vect=vecy, nbret=nbid)
        end if
        if (nvecx .lt. 0 .or. nvecy .lt. 0) then
            call utmess('F', 'MODELISA7_2')
        end if
!
! VECTEUR Z LOCAL = VECT-X VECTORIEL VECT-Y
        AS_ALLOCATE(vr=vecz, size=3)
        vecz(1) = vecx(1+1)*vecy(1+2)-vecy(1+1)*vecx(1+2)
        vecz(1+1) = vecx(1+2)*vecy(1)-vecy(1+2)*vecx(1)
        vecz(1+2) = vecx(1)*vecy(1+1)-vecy(1)*vecx(1+1)
        do in = 1, 3
            dir(1, in) = vecx(in)
            dir(2, in) = vecy(in)
            dir(3, in) = vecz(in)
        end do
    else if (method(1:7) .eq. 'AU_YANG') then
        yang = .true.
        call getvr8(' ', 'VECT_X', nbval=0, nbret=nvecx)
        nvecx = -nvecx
        if (nvecx .gt. 0) then
            call getvr8(' ', 'VECT_X', nbval=nvecx, vect=dir(1, 1), nbret=nbid)
        end if
        call getvr8(' ', 'ORIG_AXE', nbval=0, nbret=nveco)
        nveco = -nveco
        if (nveco .gt. 0) then
            call getvr8(' ', 'ORIG_AXE', nbval=nveco, vect=dir(1, 2), nbret=nbid)
        end if
        if (nvecx .lt. 0 .or. nveco .lt. 0) then
            call utmess('F', 'MODELISA7_3')
        end if
    end if
!
! VALEURS NON DEPENDANTES DE LA FREQUENCE
!
10  continue
    if (vate(1) .eq. 'SPEC_CORR_CONV_3') then
        call accep2(base(1:8), nbm, pg, phi, sphi)
    else
        call accep1(base(1:8), ligrmo, nbm, dir, yang)
    end if
!
!
! CAS SPEC_CORR_CONV_1 ET 2
    mxval = nbm*(nbm+1)/2
    chnumi = nomres//'.NUMI'
    call wkvect(chnumi, 'G V I', mxval, lnumi)
    chnumj = nomres//'.NUMJ'
    call wkvect(chnumj, 'G V I', mxval, lnumj)
    chvale = nomres//'.VALE'
    call jecrec(chvale, 'G V R', 'NU', 'DISPERSE', 'VARIABLE', &
                mxval)
    chfreq = nomres//'.DISC'
    call wkvect(chfreq, 'G V R', nbpoin, lfreq)
!
    do iff = 0, nbpoin-1
        f = finit+iff*df
        zr(lfreq+iff) = f
    end do
!
!  POUR LE CAS SPEC_CORR_CONV_3
    if (vate(1) .eq. 'SPEC_CORR_CONV_3') then
! TABLE CONTENANT LES FONCTIONS DE FORME
        is = vate(2)
        do iff = 0, nbpoin-1
            f = finit+iff*df
            zr(lfreq+iff) = f
            call evalis(is, pg, phi, sphi, f, &
                        iff, nomres)
        end do
    else
        ij = 0
        do im2 = 1, nbm
!
            do im1 = im2, nbm
                ij = ij+1
!
                zi(lnumi-1+ij) = im1
                zi(lnumj-1+ij) = im2
!
                call jecroc(jexnum(chvale, ij))
                if (im1 .eq. im2) then
                    nbabs = nbpoin
                else
                    nbabs = 2*nbpoin
                end if
!
                call jeecra(jexnum(chvale, ij), 'LONMAX', nbabs)
                call jeecra(jexnum(chvale, ij), 'LONUTI', nbabs)
                call jeveuo(jexnum(chvale, ij), 'E', lvale)
!
! BOUCLE SUR LES FREQUENCES ET REMPLISSAGE DU .VALE
! IE VALEURS DES INTERSPECTRS
!
                ier = 0
                do iff = 0, nbpoin-1
                    f = finit+iff*df
                    if (f .gt. fcoup) then
                        prs = 0.d0
                    else if (vate(1) .eq. 'SPEC_CORR_CONV_2') then
                        puls = deuxpi*f
                        call fointe('F', fonct, 1, ['PULS'], [puls], &
                                    prs, ier)

! APPELS AUX LONGUEURS DE CORRELATION
                        call fointe('F', long1f, 1, ['PULS'], [puls], &
                                    long1, ier)
                        call fointe('F', long2f, 1, ['PULS'], [puls], &
                                    long2, ier)

! APPEL A LA FONCTION DE CALCUL D ACCEPTANCE
                        call accept(f, nbm, method, im2, im1, &
                                    jc, dir, uc, ut, &
                                    long1, long2, val_spec)
                    else
                        prs = dspprs(kste, uflui, dhyd, rho, f, fcoup_red)
                        call accept(f, nbm, method, im2, im1, &
                                    jc, dir, uc, ut, &
                                    long1, long2, val_spec)

                    end if
                    if (im1 .eq. im2) then
                        zr(lvale+iff) = prs*jc
                    else
                        zr(lvale+2*iff) = prs*jc
                        zr(lvale+2*iff+1) = 0.d0
                    end if
                end do
            end do
        end do
    end if
!
    AS_DEALLOCATE(vr=vecx)
    AS_DEALLOCATE(vr=vecy)
    AS_DEALLOCATE(vr=vecz)
!
    if (vate(1) .eq. 'SPEC_CORR_CONV_3') then
    else
        call jedetc('V', '&&329', 1)
        call jedetc('V', '&&V.M', 1)
        call jedetc('V', '&&GROTAB.TAB', 1)
    end if
!
    call jedema()
end subroutine
