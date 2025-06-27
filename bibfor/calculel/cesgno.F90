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
subroutine cesgno(ces1, celfpg, ces2)
!
    implicit none
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/cesexi.h"
#include "asterfort/dismoi.h"
#include "asterfort/elraca.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jni002.h"
#include "asterfort/nuelrf.h"
!
    character(len=24) :: celfpg
    character(len=19) :: ces2, ces1
! ------------------------------------------------------------------
! BUT: TRANSFORMER UN CHAM_ELEM_S/ELGA EN CHAM_ELEM_S/ELNO
! ------------------------------------------------------------------
!     ARGUMENTS:
! CES1  IN/JXIN  K19 : CHAM_ELEM_S A TRANSFORMER
!
! CELFPG IN/JXVAR  K24 :
!    NOM DE L'OBJET DECRIVANT LES FAMILLES DE P.G. DE CES1 (OU ' ')
!    CET OBJET N'EST UTILISE QUE SI ELGA -> ELNO
!    CET OBJET EST OBTENU PAR LA ROUTINE CELFPG.F
!  ATTENTION : CET OBJET EST EFFACE PENDANT L'OPERATION.
!
!
! CES2    IN/JXVAR K19 : CHAM_ELEM_S RESULTAT
!         REMARQUE : LE CHAM_ELEM_S EST DEJA ALLOUE.
!-----------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbpgmx = 27
    character(len=8) :: elrf, fapg1, fapg(MT_NBFAMX)
    integer(kind=8) :: nbfpg, nbpg(MT_NBFAMX), ndiml, nnol, nnosl
!
    integer(kind=8) :: ima, ncmp, icmp, ino, isp, nno
    integer(kind=8) :: nbma, iret
    integer(kind=8) :: npg, ipg, nujni, nbobj
    integer(kind=8) :: jces1d, jces1l, jces1v, jces1c, iad1, nbpt1, nbsp1
    integer(kind=8) :: jces2d, jces2l, jces2v, iad2, nbpt2, nbsp2
    integer(kind=8) :: jmat, jganol, ivfl, jdfd2l, jcoopl, ipoidl, npgl, lonfam
    integer(kind=8) :: ifam, decal, jvr, idfdel, nufpg, avance, jnofpg
    character(len=8) :: ma, nomgd
    character(len=3) :: tsca
    character(len=16) :: schema
    character(len=24) :: liobj(10)
    real(kind=8) :: vrpg(nbpgmx), vrno(MT_NNOMAX), sr
    complex(kind=8) :: vcpg(nbpgmx), vcno(MT_NNOMAX), sc
    character(len=8), pointer :: cesk(:) => null()
    integer(kind=8) :: ndim
    integer(kind=8) :: nnos, ipoids, jcoopg, ivf, idfde, jdfd2, jgano
!
!
!     ------------------------------------------------------------------
    call jemarq()
!
!
!
!     1. RECUPERATION DE :
!        MA     : NOM DU MAILLAGE
!        NOMGD  : NOM DE LA GRANDEUR
!        NCMP   : NOMBRE DE CMPS DE CES1
!        TSCA   : TYPE SCALAIRE DE LA GRANDEUR : R/C/I ...
!        NBMA   : NOMBRE DE MAILLES DU MAILLAGE
!        JCES1XX : ADRESSES CES1
!        JCES2XX : ADRESSES CES2
!        JNOFPG : ADRESSE DE CELFPG
!     --------------------------------------------------------------
    call exisd('CHAM_ELEM_S', ces1, iret)
    ASSERT(iret .gt. 0)
    call jeveuo(ces1//'.CESK', 'L', vk8=cesk)
    call jeveuo(ces1//'.CESC', 'L', jces1c)
    call jeveuo(ces1//'.CESD', 'L', jces1d)
    call jeveuo(ces1//'.CESV', 'L', jces1v)
    call jeveuo(ces1//'.CESL', 'L', jces1l)
    ma = cesk(1)
    nomgd = cesk(2)
    call dismoi('NB_MA_MAILLA', ma, 'MAILLAGE', repi=nbma)
    call jelira(ces1//'.CESC', 'LONMAX', ncmp)
    call dismoi('TYPE_SCA', nomgd, 'GRANDEUR', repk=tsca)
    ASSERT(tsca .eq. 'R' .or. tsca .eq. 'C')
!
    call jeveuo(ces2//'.CESD', 'L', jces2d)
    call jeveuo(ces2//'.CESV', 'E', jces2v)
    call jeveuo(ces2//'.CESL', 'E', jces2l)
!
    ASSERT(celfpg .ne. ' ')
    call jeveuo(celfpg, 'E', jnofpg)
!
!
!
!     3. REMPLISSAGE DES OBJETS .CESL ET .CESV :
!     ------------------------------------------
!
!     POUR DES RAISONS DE PERFORMANCE, ON TRAITE TOUTES LES MAILLES
!     CORRESPONDANT AU MEME SCHEMA DE POINTS DE GAUSS
!
!     BOUCLE TANT QUE CELFPG N'EST PAS TOTALEMENT EFFACE (AVANCE=1):
!     ---------------------------------------------------------------
10  continue
    schema = ' '
    avance = 0
!
    do ima = 1, nbma
        if (zk16(jnofpg-1+ima) .eq. ' ') goto 110
        if (schema .eq. ' ') schema = zk16(jnofpg-1+ima)
        if (zk16(jnofpg-1+ima) .ne. schema) goto 110
!
        avance = avance+1
        schema = zk16(jnofpg-1+ima)
        elrf = schema(1:8)
        fapg1 = schema(9:16)
        zk16(jnofpg-1+ima) = ' '
!
!
!
        if (schema(9:12) .eq. 'XFEM') then
!
!           3.1 : CALCUL DE LA MATRICE DE PASSAGE GA->NO
!                 (ON NE LE FAIT QUE POUR LA 1ERE MAILLE DU SCHEMA)
!           OUT : NPG,NNO,JMAT
!           --------------------------------------------------------
!
            if (avance .eq. 1) then
                ASSERT(schema(1:8) .eq. elrf)
!
! --------- Get list of integration schemes of geometric support
                call elraca(elrf, nbfpg_=nbfpg, fapg_=fapg, nbpg_=nbpg, ndim_=ndiml, &
                            nno_=nnol, nnos_=nnosl)
!
! --------- Get index for integration scheme
                nufpg = indik8(fapg1, fapg1, 1, nbfpg)
                ASSERT(nufpg .gt. 0)
!
                call jeveuo('&INEL.'//elrf//'.ELRA_R', 'L', jvr)
                decal = 0
                do ifam = 1, nufpg-1
                    npgl = nbpg(ifam)
!
                    lonfam = npgl
                    lonfam = lonfam+npgl*ndiml
                    lonfam = lonfam+npgl*nnol
                    lonfam = lonfam+npgl*nnol*ndiml
                    lonfam = lonfam+npgl*nnol*ndiml*ndiml
                    lonfam = lonfam+2+npgl*nnol
!
                    decal = decal+lonfam
                end do
!
                npgl = nbpg(nufpg)
!
                ipoidl = jvr+decal
                jcoopl = ipoidl+npgl
                ivfl = jcoopl+npgl*ndiml
                idfdel = ivfl+npgl*nnol
                jdfd2l = idfdel+npgl*nnol*ndiml
                jganol = jdfd2l+npgl*nnol*ndiml*ndiml
!
                ndim = ndiml
                nnos = nnosl
                nno = nnol
                npg = npgl
                ipoids = ipoidl
                jcoopg = jcoopl
                ivf = ivfl
                idfde = idfdel
                jdfd2 = jdfd2l
                jgano = jganol
!
                ASSERT(nno .le. MT_NNOMAX)
                ASSERT(npg .le. nbpgmx)
            end if
!
            nbpt1 = zi(jces1d-1+5+4*(ima-1)+1)
            nbsp1 = zi(jces1d-1+5+4*(ima-1)+2)
            nbpt2 = zi(jces2d-1+5+4*(ima-1)+1)
            nbsp2 = zi(jces2d-1+5+4*(ima-1)+2)
            ASSERT(nbsp1 .eq. nbsp2)
            ASSERT(nbpt2 .eq. nno)
!
!
            do icmp = 1, ncmp
                call cesexi('C', jces1d, jces1l, ima, 1, &
                            1, icmp, iad1)
                if (iad1 .le. 0) goto 300
                do isp = 1, nbsp1
                    do ino = 1, nno
                        call cesexi('C', jces2d, jces2l, ima, ino, &
                                    isp, icmp, iad2)
                        ASSERT(iad2 .lt. 0)
                        if (tsca .eq. 'R') then
                            zr(jces2v-1-iad2) = 0.d0
                        else
                            zc(jces2v-1-iad2) = 0.d0
                        end if
                        zl(jces2l-1-iad2) = .true.
                    end do
                end do
300             continue
            end do
            goto 110
!
!
        else
!
!
!           3.1 : CALCUL DE LA MATRICE DE PASSAGE GA->NO
!                 (ON NE LE FAIT QUE POUR LA 1ERE MAILLE DU SCHEMA)
!           OUT : NPG,NNO,JMAT
!           --------------------------------------------------------
            if (avance .eq. 1) then
!
! --------- Get list of integration schemes of geometric support
                call elraca(elrf, nbfpg_=nbfpg, fapg_=fapg, nbpg_=nbpg, ndim_=ndiml, &
                            nno_=nnol, nnos_=nnosl)
!
! --------- Get index for integration scheme
                nufpg = indik8(fapg, fapg1, 1, nbfpg)
                ASSERT(nufpg .gt. 0)
!
                call nuelrf(elrf, nujni)
                ASSERT(nujni .eq. 2)
                call jni002(elrf, 10, liobj, nbobj)
                call jeveuo('&INEL.'//elrf//'.ELRA_R', 'L', jvr)
!
                decal = 0
                do ifam = 1, nufpg-1
                    npgl = nbpg(ifam)
!
                    lonfam = npgl
                    lonfam = lonfam+npgl*ndiml
                    lonfam = lonfam+npgl*nnol
                    lonfam = lonfam+npgl*nnol*ndiml
                    lonfam = lonfam+npgl*nnol*ndiml*ndiml
                    lonfam = lonfam+2+npgl*nnol
!
                    decal = decal+lonfam
                end do
!
                npgl = nbpg(nufpg)
!
                ipoidl = jvr+decal
                jcoopl = ipoidl+npgl
                ivfl = jcoopl+npgl*ndiml
                idfdel = ivfl+npgl*nnol
                jdfd2l = idfdel+npgl*nnol*ndiml
                jganol = jdfd2l+npgl*nnol*ndiml*ndiml
                nno = nint(zr(jganol-1+1))
                npg = nint(zr(jganol-1+2))
                jmat = jganol+2
                ASSERT(nno .le. MT_NNOMAX)
                ASSERT(npg .le. nbpgmx)
            end if
!
!
!           3.2 : MULTIPLICATION PAR LA MATRICE
!           ------------------------------------
            nbpt1 = zi(jces1d-1+5+4*(ima-1)+1)
            nbsp1 = zi(jces1d-1+5+4*(ima-1)+2)
            nbpt2 = zi(jces2d-1+5+4*(ima-1)+1)
            nbsp2 = zi(jces2d-1+5+4*(ima-1)+2)
            ASSERT(nbsp1 .eq. nbsp2)
            ASSERT(nbpt1 .eq. npg)
            ASSERT(nbpt2 .eq. nno)
!
!
            do icmp = 1, ncmp
                call cesexi('C', jces1d, jces1l, ima, 1, &
                            1, icmp, iad1)
                if (iad1 .le. 0) goto 100
!
                do isp = 1, nbsp1
!
!               -- RECOPIE DANS VXPG :
                    do ipg = 1, npg
                        call cesexi('C', jces1d, jces1l, ima, ipg, &
                                    isp, icmp, iad1)
                        ASSERT(iad1 .gt. 0)
                        if (tsca .eq. 'R') then
                            vrpg(ipg) = zr(jces1v-1+iad1)
!
                        else
                            vcpg(ipg) = zc(jces1v-1+iad1)
                        end if
                    end do
!
!               -- MULTIPLICATION :
                    if (tsca .eq. 'R') then
                        do ino = 1, nno
                            sr = 0.d0
                            do ipg = 1, npg
                                sr = sr+zr(jmat-1+(ino-1)*npg+ipg)*vrpg(ipg)
                            end do
                            vrno(ino) = sr
                        end do
!
                    else
                        do ino = 1, nno
                            sc = dcmplx(0.d0, 0.d0)
                            do ipg = 1, npg
                                sc = sc+zr(jmat-1+(ino-1)*npg+ipg)*vcpg(ipg)
                            end do
                            vcno(ino) = sc
                        end do
                    end if
!
!               -- RECOPIE DE VXNO :
                    do ino = 1, nno
                        call cesexi('C', jces2d, jces2l, ima, ino, &
                                    isp, icmp, iad2)
                        ASSERT(iad2 .lt. 0)
                        if (tsca .eq. 'R') then
                            zr(jces2v-1-iad2) = vrno(ino)
                        else
                            zc(jces2v-1-iad2) = vcno(ino)
                        end if
                        zl(jces2l-1-iad2) = .true.
                    end do
                end do
100             continue
            end do
        end if
!
110     continue
    end do
    if (avance .gt. 0) goto 10
!
    call jedema()
end subroutine
