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

subroutine dc_polycristal(nboccp, sdcomp)
!
! person_in_charge: jean-luc.flejou at edf.fr
!
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/eulnau.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
!
    integer(kind=8), intent(in) :: nboccp
    character(len=8), intent(in) :: sdcomp
!
! --------------------------------------------------------------------------------------------------
!
! DEFI_COMPOR
!
! POLYCRISTAL
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: mono, chaine
    character(len=16) :: loca, noms(6)
    character(len=24) :: nomvar(100)
    real(kind=8) :: fvol, orie(3), dl, da, euler(3), fvolt, mu_loca
    integer(kind=8) :: iocc, nloc, ndl, nda, itbint, nums(3)
    integer(kind=8) :: i, nmono, imi, ipk, ipi, ipr, iorie, irra
    integer(kind=8) :: ncpri, ncprk, ncprr, jcprk, jcprr, jcpri, nvit, lmk, ifvol, ipl
    integer(kind=8) :: imono, nbmono, nvloc, indvar
    integer(kind=8) :: nbsyst, nbsysm
    character(len=24), pointer :: cprk(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    call getvtx(' ', 'LOCALISATION', scal=loca, nbret=nloc)
    dl = 0.d0
    da = 0.d0
    nvloc = 1
    call getvr8(' ', 'MU_LOCA', scal=mu_loca, nbret=nloc)
    if (loca .eq. 'BETA') then
        call getvr8(' ', 'DL', scal=dl, nbret=ndl)
        call getvr8(' ', 'DA', scal=da, nbret=nda)
        nvloc = nvloc+2
    end if
!
!
!     organisation de CPRI :
!     1 : TYPE =2 pour POLYCRISTAL
!     2 : NBPHAS pour POLYCRISTAL
!     3 : NVITOT pour POLYCRISTAL
!     4 : NOMBRE DE MONOCRISTAUX différents
!     5 : NBFAMILLES DE SYS GLIS pour Phase 1
!     6 : Numero du MONO 1
!     7 : NVI du Mono 1
!     8 : NBFAMILLES DE SYS GLIS pour Phase 2
!     9 : Numero du MONO 2
!     10 : NVI du Mono 2
!      etc...
!     avant dernier : dimension de CPRK
!     nombre de paramètres de localisation
!
    ncpri = 4+3*nboccp+1+1+1
    call wkvect(sdcomp//'.CPRI', 'G V I', ncpri, ipi)
    zi(ipi) = 2
    zi(ipi+1) = nboccp
    call wkvect('&&OP0059.LISTEMONO', 'V V K8', nboccp, ipl)
!
    nbmono = 0
    ncprk = 0
!
    do iocc = 1, nboccp
        call getvid('POLYCRISTAL', 'MONOCRISTAL', iocc=iocc, scal=mono, nbret=nmono)
!        On ne stocke pas les doublons
        imono = indik8(zk8(ipl), mono, 1, nbmono)
        if (imono .eq. 0) then
            nbmono = nbmono+1
            zk8(ipl-1+nbmono) = mono
            zi(ipi-1+4+3*(iocc-1)+2) = nbmono
            call jelira(mono//'.CPRK', 'LONMAX', lmk)
            ncprk = ncprk+lmk+2
        else
            zi(ipi-1+4+3*(iocc-1)+2) = imono
        end if
    end do
    ncprk = ncprk+1
    if (nbmono .gt. 5) then
        call utmess('F', 'COMPOR2_16', si=itbint)
    else
        zi(ipi-1+4) = nbmono
    end if
    irra = 0
!
!     organisation de CPRK :
!     On ne stocke que les monocristaux DIFFERENTS
!     1   : Nom méthode localisation
!     2   : Nom Monocristal 1 + NBFAM + CPRK du monocristal 1
!     n+2 : Nom Monocristal 2 + NBFAM + CPRK du monocristal 2
!     ..: etc...
    call wkvect(sdcomp//'.CPRK', 'G V K24', ncprk, ipk)
    jcprk = 1
    itbint = 0
    do imono = 1, nbmono
        mono = zk8(ipl-1+imono)
        call jelira(mono//'.CPRK', 'LONMAX', lmk)
        call jeveuo(mono//'.CPRK', 'L', vk24=cprk)
        call jeveuo(mono//'.CPRI', 'L', imi)
!        RECOPIE DU VECTEUR K16 DU MONOCRISTAL DANS CELUI DU POLY
        zk24(ipk-1+jcprk+1) = mono
        write (zk24(ipk-1+jcprk+2), '(I24)') zi(imi-1+5)
        do i = 1, lmk
            zk24(ipk-1+jcprk+2+i) = cprk(i)
        end do
        jcprk = jcprk+lmk+2
        if (cprk(3) .eq. 'MONO_DD_CC_IRRA') irra = 1
        if (cprk(3) .eq. 'MONO_DD_CFC_IRRA') irra = 2
!
    end do
!
    ncprr = 4*nboccp+3
!
    call wkvect(sdcomp//'.CPRR', 'G V R', ncprr, ipr)
    jcprr = 0
    jcpri = 4
    nvit = 0
    fvolt = 0.d0
    nbsyst = 0
    nbsysm = 0
    do iocc = 1, nboccp
        imono = zi(ipi-1+4+3*(iocc-1)+2)
        mono = zk8(ipl-1+imono)
        call jeveuo(mono//'.CPRI', 'L', imi)
        zi(ipi-1+jcpri+1) = zi(imi-1+5)
        zi(ipi-1+jcpri+3) = zi(imi-1+7)
!
!        NVI DU MONOCRISTAL : 6+4*NS + 3 + (NS SI IRRA) +3 )
!       (EVP + NS(ALPHAS, GAMMAS, PS) + RHO_IRRA + 3)
!        NOMBRE DE VAR INT POLYCRISTAL SANS IRRA :
!         7   + 6*NG  +NG*(NS*(ALPHAS,GAMMAS,PS))+6*NG+1
!        EVP+P+EVPG*NG+NG*(NS*(ALPHAS,GAMMAS,PS))+SIG*NG+1)
!        =7   + 6*NG  +NG*(NMONO-9)+6*NG+1=NG*(NMONO+3)+8
!        NOMBRE DE VAR INT POLYCRISTAL AVEC IRRA :
!         7   + 6*NG  +NG*(NS*(ALPHAS,GAMMAS,PS))+12*NG+6*NG+1
!        =7   + 6*NG  +NG*(NMONO-9)+6*NG+1=NG*(NMONO+3)+8
!
        nbsysm = max(nbsysm, zi(imi-1+8))
        nbsyst = nbsyst+zi(imi-1+8)
!        ON ENLEVE LES TAUS ET LES 3 VARIABLES INTERNES P,SCLIV,INDIC
        nvit = nvit+zi(imi-1+7)+3-zi(imi-1+8)
        jcpri = jcpri+3
        call getvr8('POLYCRISTAL', 'FRAC_VOL', iocc=iocc, scal=fvol, nbret=ifvol)
        call getvr8('POLYCRISTAL', 'ANGL_REP', iocc=iocc, nbval=3, vect=orie, &
                    nbret=iorie)
        if (iorie .eq. 0) then
            call getvr8('POLYCRISTAL', 'ANGL_EULER', iocc=iocc, nbval=3, vect=euler, &
                        nbret=iorie)
            call eulnau(euler, orie)
        end if
        fvolt = fvolt+fvol
        zr(ipr-1+jcprr+1) = fvol
        zr(ipr-1+jcprr+2) = orie(1)
        zr(ipr-1+jcprr+3) = orie(2)
        zr(ipr-1+jcprr+4) = orie(3)
        jcprr = jcprr+4
    end do
!
    if (abs(fvolt-1.d0) .gt. 1.d-3) then
        call utmess('F', 'COMPOR2_8', sr=fvolt)
    end if
!
    zr(ipr-1+jcprr+1) = mu_loca
    zr(ipr-1+jcprr+2) = dl
    zr(ipr-1+jcprr+3) = da
!
!      NOMBRE DE VAR INT TOTAL + 8 (TENSEUR B OU EVP + NORME+INDIC)
    nvit = nvit+8
    zi(ipi-1+3) = nvit
    zi(ipi-1+ncpri-2) = jcprk
    zi(ipi-1+ncpri-1) = nvloc
!
    zk24(ipk) = loca
!
!     IMPRESSION DES VARIABLES INTERNES
    indvar = 0
    noms(1) = 'POLYCRISTAL'
    noms(2) = loca
    nums(1) = nboccp
    nums(2) = nbmono
    nums(3) = nvit
    call utmess('I', 'COMPOR4_53', sk=loca, ni=3, vali=nums)
!
    nomvar(1) = 'EPSPXX'
    nomvar(2) = 'EPSPYY'
    nomvar(3) = 'EPSPZZ'
    nomvar(4) = 'EPSPXY'
    nomvar(5) = 'EPSPXZ'
    nomvar(6) = 'EPSPYZ'
    nomvar(7) = 'EPSPEQ'
    do i = 1, 7
        call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=i)
    end do
    indvar = indvar+7
    call utmess('I', 'COMPOR4_54', si=indvar+1)
    nomvar(1) = 'EPSPXX(GRAIN_I)'
    nomvar(2) = 'EPSPYY(GRAIN_I)'
    nomvar(3) = 'EPSPZZ(GRAIN_I)'
    nomvar(4) = 'EPSPXY(GRAIN_I)'
    nomvar(5) = 'EPSPXZ(GRAIN_I)'
    nomvar(6) = 'EPSPYZ(GRAIN_I)'
    do i = 1, 6
        call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=indvar+i)
    end do
    indvar = indvar+6*nboccp
    call utmess('I', 'COMPOR4_56', si=indvar)
!
    call utmess('I', 'COMPOR4_55', si=indvar+1)
!
    do i = 1, nbsysm
        call codent(i, 'G', chaine)
        nomvar(3*i-2) = 'ALPHA'//chaine
        nomvar(3*i-1) = 'GAMMA'//chaine
        nomvar(3*i) = 'P'//chaine
    end do
    do i = 1, 3*nbsysm
        call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=indvar+i)
    end do
!
    indvar = indvar+3*nbsyst
!
    call utmess('I', 'COMPOR4_56', si=indvar)
!
    if (irra .eq. 1) then
        call utmess('I', 'COMPOR4_54', si=indvar+1)
        do i = 1, 12
            call codent(i, 'G', chaine)
            nomvar(i) = 'RHO_IRRA_'//chaine
        end do
        do i = 1, 12
            call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=indvar+i)
        end do
        indvar = indvar+12*nboccp
        call utmess('I', 'COMPOR4_56', si=indvar)
    end if
!
    if (irra .eq. 2) then
        call utmess('I', 'COMPOR4_54', si=indvar+1)
        do i = 1, 12
            call codent(i, 'G', chaine)
            nomvar(i) = 'RHO_LOOPS_'//chaine
        end do
        do i = 1, 12
            call codent(i, 'G', chaine)
            nomvar(12+i) = 'PHI_VOIDS_'//chaine
        end do
        do i = 1, 24
            call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=indvar+i)
        end do
        indvar = indvar+24*nboccp
        call utmess('I', 'COMPOR4_56', si=indvar)
    end if
!
    call utmess('I', 'COMPOR4_54', si=indvar+1)
!
    do i = 1, 6
        nomvar(1) = 'SIGMAXX(GRAIN_I)'
        nomvar(2) = 'SIGMAYY(GRAIN_I)'
        nomvar(3) = 'SIGMAZZ(GRAIN_I)'
        nomvar(4) = 'SIGMAXY(GRAIN_I)'
        nomvar(5) = 'SIGMAXZ(GRAIN_I)'
        nomvar(6) = 'SIGMAYZ(GRAIN_I)'
    end do
    do i = 1, 6
        call utmess('I', 'COMPOR4_20', sk=nomvar(i), si=indvar+i)
    end do
    indvar = indvar+6*nboccp
    call utmess('I', 'COMPOR4_56', si=indvar)
    indvar = indvar+1
    nomvar(1) = 'INDIPLAS'
    call utmess('I', 'COMPOR4_55', sk=nomvar(1), si=nvit)
    ASSERT(indvar .eq. nvit)
! FIN ------------------------------------------------------------------
    call jedema()
end subroutine
