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
!
subroutine manopg(model, ligrez, optioz, paramz, mnogaz)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/alchml.h"
#include "asterfort/assert.h"
#include "asterfort/cescre.h"
#include "asterfort/cesexi.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/elrfno.h"
#include "asterfort/elraga.h"
#include "asterfort/elrfvf.h"
#include "asterfort/indk32.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/cormgi.h"
#include "asterfort/initel.h"
#include "asterfort/jedup1.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/manopx.h"
#include "asterfort/nbelem.h"
#include "asterfort/nbptca.h"
#include "asterfort/typele.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
!
    character(len=8), intent(in) :: model
    character(len=*) :: ligrez, mnogaz, optioz, paramz
! ------------------------------------------------------------------
! BUT: CREER LE CHAM_ELEM_S MNOGAZ QUI CONTIENDRA LA MATRICE
!      DE PASSAGE NOEUDS -> POINTS DE GAUSS POUR LES ELEMENTS
!      DU LIGREL ET POUR LA FAMILLE (OPTIOZ/PARAMZ)
! ------------------------------------------------------------------
!     ARGUMENTS:
! LIGREZ  IN/JXIN  K19 : LIGREL
! OPTIOZ,PARAMZ  IN  K* : OPTION ET PARAMETRE PERMETTANT DE DETERMINER
!                         LA FAMILLE DE PG UTILISEE.
! MNOGAZ  IN/JXOUT K19 : CHAM_ELEM_S (VARI_R) DE TYPE 'ELEM'
! ------------------------------------------------------------------
! REMARQUES :
!  MNOGAZ(IMA) EST UN VECTEUR V DE REELS DE DIMENSION 2 + NBNO*NBPG
!    V(1) : NBNO
!    V(2) : NBPG
!    V(2+NBNO*(IPG-1)+INO) : MATRICE DE PASSAGE (IPG,INO)
!
!  ATTENTION :
!     1) LES MAILLES TARDIVES SONT IGNOREES.
!     2) POUR ECONOMISER LE VOLUME DE MNOGAZ, ON UTILISE LE FAIT QUE
!        LES MATRICES DE PASSAGE DES ELEMENTS D'UN MEME GREL SONT
!        IDENTIQUES CAR ELLES NE DEPENDENT QUE DE L'ELREFA.
!
!        EXCEPTION : LA FAMILLE XFEM DES ELEMENTS XFEM EST UNE FAMILLE
!        SPECIALE DONT LA POSITION DES POINTS DEPEND DU DECOUPAGE DE
!        CHAQUE ELEMENT. POUR CETTE FAMILLE, ON NE PEUT PAS FAIRE
!        L'ECONOMIE DANS MNOGAZ.
!
!        ON UTILISE LA CONVENTION :
!          SI MNOGAZ(IMA,1) > 0 : LA MAILLE IMA EST LA 1ERE D'UN GREL
!              SA MATRICE EST STOCKEE ET ELLE SERT DE REFERENCE POUR
!              LES AUTRES
!          SI MNOGAZ(IMA,1) < 0 : LA MAILLE IMA N'EST PAS LA 1ERE
!              D'UN GREL. MNOGAZ(IMA,1)= -IMAREF
!              IMAREF EST LA MAILLE DE REFERENCE POUR IMA
!
!-----------------------------------------------------------------------
!
    integer(kind=8) :: nbpgmx, nbflmx
    parameter(nbpgmx=27, nbflmx=20)
    integer(kind=8) :: nbma, ima, jcesd, jcesl, iad, jnbpg
    integer(kind=8) :: ilcnx1, k, jfpgl
    integer(kind=8) :: nec, kfpg, ndim, nno, npg, kp, ino
    integer(kind=8) ::  nbgrel, nel, nute, imolo, jmolo, jecono
    integer(kind=8) :: igr, iel, lont1
    integer(kind=8) :: jnbno, jdime, iret, ncpmax, nbfam, kfam, nbpgt, iad0, iad1
    integer(kind=8) :: nblfpg, nuflpg, nufgpg, jliel, jliel1
    integer(kind=8) :: jcesgl, jcesgv, jcesgd, nbpt, nbsp
    character(len=3) :: exixfm
    character(len=8) :: ma, nomgd, famil, elrefe, param
    character(len=8) :: lielrf(nbflmx), lifapg(nbflmx)
    character(len=16) :: pheno, option, nomte, nofpg, valk(2)
    character(len=19) :: mnoga, ligrel, celmod, ligre1, chsgeo
    character(len=24) :: obnbpg, obnbno, obdime, kecono
    character(len=32) :: noflpg
    real(kind=8) :: xpg(3*nbpgmx), ff(MT_NNOMAX), poipg(MT_NNOMAX)
    aster_logical :: econom
    integer(kind=8), pointer :: nolocfpg(:) => null()
    integer(kind=8), pointer :: celd(:) => null()
    character(len=32), pointer :: pnlocfpg(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    integer(kind=8), pointer :: mailref(:) => null()
!     ------------------------------------------------------------------
    call jemarq()
!
    mnoga = mnogaz
    ligrel = ligrez
    option = optioz
    param = paramz
!
    obnbpg = '&&MANOPG.NBPG'
    obnbno = '&&MANOPG.NBNO'
    obdime = '&&MANOPG.DIME'
!
    call jelira(ligrel//'.LIEL', 'NMAXOC', nbgrel)
    call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=ma)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call dismoi('PHENOMENE', ligrel, 'LIGREL', repk=pheno)
    call jeveuo(jexatr(ma//'.CONNEX', 'LONCUM'), 'L', ilcnx1)
!
    call jeveuo('&CATA.TE.PNLOCFPG', 'L', vk32=pnlocfpg)
    call jelira('&CATA.TE.NOLOCFPG', 'LONMAX', nblfpg)
    call jeveuo('&CATA.TE.NOLOCFPG', 'L', vi=nolocfpg)
!
!
!   0.1 CALCUL DE '&&MANOPG.ECONO' ET '&&MANOPG.CHSGEO'
!   ------------------------------------------------------------------
!     '.ECONO': V(IGR) = 1 => LE GREL IGR EST STOCKE 'ECONOMIQUE'
    kecono = '&&MANOPG.ECONO'
    chsgeo = '&&MANOPG.CHSGEO'
    call manopx(model, ligrel, option, param, chsgeo, exixfm, &
                kecono)
    call jeveuo(kecono, 'L', jecono)
    if (exixfm .eq. 'OUI') then
        call jeveuo(chsgeo//'.CESD', 'L', jcesgd)
        call jeveuo(chsgeo//'.CESV', 'L', jcesgv)
        call jeveuo(chsgeo//'.CESL', 'L', jcesgl)
    end if
!
!
!   0.2 ON FABRIQUE UN "FAUX" LIGREL (LIGRE1) N'AYANT QUE LE NOMBRE
!      NECESSAIRE D'ELEMENTS PAR GREL POUR DIMINUER LA TAILLE DE MNOGA
!   ------------------------------------------------------------------
!
!   '&&MANOPG.MAILREF': OBJET DONNANT POUR CHAQUE MAILLE SA MAILLE
!                       DE REFERENCE
    AS_ALLOCATE(vi=mailref, size=nbma)
!
    ligre1 = '&&MANOPG.LIGRE1'
    ASSERT(nbgrel .gt. 0)
    call jecrec(ligre1//'.LIEL', 'V V I', 'NU', 'CONTIG', 'VARIABLE', &
                nbgrel)
    lont1 = 0
    do igr = 1, nbgrel
        econom = zi(jecono-1+igr) .eq. 1
        nel = nbelem(ligrel, igr)
        if (econom .and. nel .gt. 1) then
            lont1 = lont1+2
        else
            lont1 = lont1+nel+1
        end if
    end do
!
    call jeecra(ligre1//'.LIEL', 'LONT', lont1, ' ')
    do igr = 1, nbgrel
        econom = zi(jecono-1+igr) .eq. 1
        nel = nbelem(ligrel, igr)
        if (econom .and. nel .gt. 1) then
            call jecroc(jexnum(ligre1//'.LIEL', igr))
            call jeecra(jexnum(ligre1//'.LIEL', igr), 'LONMAX', 2)
            call jeveuo(jexnum(ligre1//'.LIEL', igr), 'E', jliel1)
            call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', jliel)
            zi(jliel1-1+1) = zi(jliel-1+1)
            zi(jliel1-1+2) = zi(jliel-1+nel+1)
            do iel = 1, nel
                ima = zi(jliel-1+iel)
                if (ima .lt. 0) cycle
                if (iel .eq. 1) then
                    mailref(ima) = +zi(jliel-1+1)
                else
                    mailref(ima) = -zi(jliel-1+1)
                end if
            end do
        else
            call jecroc(jexnum(ligre1//'.LIEL', igr))
            call jeecra(jexnum(ligre1//'.LIEL', igr), 'LONMAX', nel+1)
            call jeveuo(jexnum(ligre1//'.LIEL', igr), 'E', jliel1)
            call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', jliel)
            zi(jliel1-1+nel+1) = zi(jliel-1+nel+1)
            do iel = 1, nel
                zi(jliel1-1+iel) = zi(jliel-1+iel)
                ima = zi(jliel-1+iel)
                if (ima .lt. 0) cycle
                mailref(ima) = +zi(jliel-1+iel)
            end do
        end if
    end do
    call jedup1(ligrel//'.LGRF', 'V', ligre1//'.LGRF')
    call jedup1(ligrel//'.NBNO', 'V', ligre1//'.NBNO')
    call jedup1(ligrel//'.NEMA', 'V', ligre1//'.NEMA')
    call cormgi('V', ligre1)
    call initel(ligre1)

!
!
!
!   1. ON RECUPERE LES NOMBRES DE POINTS DE GAUSS ET DE NOEUDS :
!   ------------------------------------------------------------
    call nbptca(ligre1, option, param, obnbpg, obnbno)
    call jeveuo(obnbpg, 'L', jnbpg)
    call jeveuo(obnbno, 'L', jnbno)
!
!
!   2. ALLOCATION DU CHAM_ELEM_S MNOGA :
!   ---------------------------------------------------------------
    call wkvect(obdime, 'V V I', nbma, jdime)
    ncpmax = 0
    do ima = 1, nbma
        zi(jdime-1+ima) = zi(jnbpg-1+ima)*zi(jnbno-1+ima)+2
        ncpmax = max(ncpmax, zi(jdime-1+ima))
    end do
    call cescre('V', mnoga, 'ELEM', ma, 'VARI_R', &
                -ncpmax, ' ', [-1], [-1], zi(jdime))
!
!
!   3. ALLOCATION D'UN CHAMP MODELE POUR DETERMINER LES FAMILLES
!      DE POINTS DE GAUSS UTILISEES.
!   ---------------------------------------------------------------
    celmod = '&&MANOPG.CELMOD'
    call alchml(ligre1, option, param, 'V', celmod, &
                iret, ' ')
    if (iret .ne. 0) then
        valk(1) = param
        valk(2) = option
        call utmess('F', 'CALCULEL7_7', nk=2, valk=valk)
    end if
    call jeveuo(celmod//'.CELD', 'L', vi=celd)
!
!
!   4. REMPLISSAGE DE MNOGA :
!   ---------------------------------------------------------------
    call jeveuo(mnoga//'.CESD', 'L', jcesd)
    call jeveuo(mnoga//'.CESL', 'E', jcesl)
    call jeveuo(mnoga//'.CESV', 'E', vr=cesv)
!
!
    do igr = 1, nbgrel
        econom = zi(jecono-1+igr) .eq. 1
        call jeveuo(jexnum(ligrel//'.LIEL', igr), 'L', jliel)
        nel = nbelem(ligrel, igr)
        nute = typele(ligrel, igr)
        call jenuno(jexnum('&CATA.TE.NOMTE', nute), nomte)
        imolo = celd(celd(4+igr)+2)
        if (imolo .eq. 0) goto 1
!
        call jeveuo(jexnum(ligre1//'.LIEL', igr), 'L', jliel1)
!
!
!       4.1 DETERMINATION DE LA LISTE DES FAMILLES DE PG :
!       -----------------------------------------------------------
!       => NBFAM, LIELRF, LIFAPG
        call jeveuo(jexnum('&CATA.TE.MODELOC', imolo), 'L', jmolo)
        call dismoi('NOM_GD', celmod, 'CHAMP', repk=nomgd)
        call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
        kfpg = zi(jmolo-1+4+nec+1)
!
!       -- FAMILLE "LISTE"
        if (kfpg .lt. 0) then
!          FAMILLE "LISTE" :
            call jelira(jexnum('&CATA.TE.FPG_LISTE', -kfpg), 'LONMAX', nbfam)
            nbfam = nbfam-1
            ASSERT(nbfam .le. nbflmx)
            call jeveuo(jexnum('&CATA.TE.FPG_LISTE', -kfpg), 'L', jfpgl)
            elrefe = zk8(jfpgl-1+nbfam+1)
            do k = 1, nbfam
                lielrf(k) = elrefe
                noflpg = nomte//elrefe//zk8(jfpgl-1+k)
                nuflpg = indk32(pnlocfpg, noflpg, 1, nblfpg)
                nufgpg = nolocfpg(nuflpg)
                call jenuno(jexnum('&CATA.TM.NOFPG', nufgpg), nofpg)
                ASSERT(elrefe .eq. nofpg(1:8))
                lifapg(k) = nofpg(9:16)
            end do
!
!       -- FAMILLE "ORDINAIRE"
        else
            nbfam = 1
            call jenuno(jexnum('&CATA.TM.NOFPG', kfpg), nofpg)
            lielrf(1) = nofpg(1:8)
            lifapg(1) = nofpg(9:16)
        end if
!
!
!       4.2 BOUCLE SUR LA/LES FAMILLE(S) :
!       ------------------------------------------
        nbpgt = 0
        do kfam = 1, nbfam
            elrefe = lielrf(kfam)
            famil = lifapg(kfam)

            call elrfno(elrefe, ndim=ndim, nno=nno)
!
!         4.2.1 CALCUL DE NPG ET XPG (SI FAMILLE != XFEM)
!         ------------------------------------------------
            if (famil(1:4) .ne. 'XFEM') then
                call elraga(elrefe, famil, ndim, npg, xpg, &
                            poipg)
                ASSERT(npg .le. nbpgmx)
!
!         4.2.2 CALCUL DE NPG (SI FAMILLE == 'XFEM...')
!         ------------------------------------------------
            else
                ima = zi(jliel-1+1)
                npg = zi(jcesgd-1+5+4*(ima-1)+1)
            end if
!
!         4.2.3 ECRITURE DANS MANOPG :
!         ------------------------------------------
            do iel = 1, nel
                ima = zi(jliel-1+iel)
                if (ima .lt. 0) goto 3
!
                call cesexi('C', jcesd, jcesl, ima, 1, &
                            1, 1, iad0)
                iad0 = abs(iad0)
                ASSERT(iad0 .gt. 0)
!
!           -- SI CE N'EST PAS UNE MAILLE DE REFERENCE :
                if (mailref(ima) .lt. 0) then
                    zl(jcesl-1+iad0-1+1) = .true.
                    cesv(iad0-1+1) = mailref(ima)
                    goto 3
                end if
                ASSERT(mailref(ima) .eq. ima)
!
!
!           -- LES 2 PREMIERES CMPS : NNO ET NPG :
                if (kfam .eq. 1) then
                    zl(jcesl-1+iad0-1+1) = .true.
                    zl(jcesl-1+iad0-1+2) = .true.
                    cesv(iad0-1+1) = nno
                    cesv(iad0-1+2) = npg
                else
                    ASSERT(nint(cesv(iad0-1+1)) .eq. nno)
                    cesv(iad0-1+2) = cesv(iad0-1+2)+npg
                end if
!
                call cesexi('C', jcesd, jcesl, ima, 1, &
                            1, 2+nno*(nbpgt+npg), iad)
                ASSERT(iad .lt. 0)
                iad = iad0+2+nbpgt*nno
!
!           -- LES NNO*NPG AUTRES CMPS :
                do kp = 1, npg
                    if (famil(1:4) .eq. 'XFEM') then
                        ASSERT(.not. econom)
                        nbpt = zi(jcesgd-1+5+4*(ima-1)+1)
                        nbsp = zi(jcesgd-1+5+4*(ima-1)+2)
                        ASSERT(nbpt .eq. npg)
                        ASSERT(nbsp .eq. 1)
                        call cesexi('C', jcesgd, jcesgl, ima, kp, &
                                    1, 1, iad1)
                        ASSERT(iad1 .gt. 0)
                        call elrfvf(elrefe, zr(jcesgv-1+iad1), ff, nno)
                    else
                        call elrfvf(elrefe, xpg(ndim*(kp-1)+1), ff, nno)
                    end if
                    do ino = 1, nno
                        zl(jcesl-1+iad-1+nno*(kp-1)+ino) = .true.
                        cesv(iad-1+nno*(kp-1)+ino) = ff(ino)
                    end do
                end do
!
3               continue
            end do
            nbpgt = nbpgt+npg
        end do
1       continue
    end do
!
!
    call detrsd('CHAMP', celmod)
    call detrsd('CHAMP', chsgeo)

    call detrsd('LIGREL', ligre1)
    call jedetr(obnbpg)
    call jedetr(obnbno)
    call jedetr(obdime)
    AS_DEALLOCATE(vi=mailref)
    call jedetr(kecono)
!
    call jedema()
end subroutine
