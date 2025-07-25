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
subroutine giecas(nfic, ndim, nbobj)
    implicit none
!
!     ARGUMENTS:
!     ----------
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/giecma.h"
#include "asterfort/giinco.h"
#include "asterfort/jacopo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
#include "asterfort/uttrii.h"
#include "asterfort/wkvect.h"
    integer(kind=8) :: nfic, ndim, nbobj
! ----------------------------------------------------------------------
!     BUT: ECRIRE LE FICHIER DE MAILLAGE ASTER A PARTIR DES OBJETS
!          CREES PAR GILIRE ( '&&GILIRE.....')
!
!     IN : NFIC : UNITE D'ECRITURE
!          NDIM : DIMENSION DU PROBLEME (2D OU 3D)
!          NBOBJ: NOMBRE D'OBJETS (AU SENS GIBI)
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: vali
!     VARIABLES LOCALES:
    character(len=7) :: k7bid, k7nom(7)
    character(len=8) :: tymail, nomobj, nomno, k8nom(7), nomobg
    aster_logical :: magoui, trouve, indir
!
    character(len=1) :: cbid
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ianema
    integer(kind=8) :: iaptin, icok
    integer(kind=8) :: icoma, iecrit, iecrma, ii, il, ima, imb
    integer(kind=8) :: ino, inu, inutri, iret, itot, itrnu, j
    integer(kind=8) :: jj, k, kk, l, nbelc, nbele, nbelim
    integer(kind=8) :: nbelt, nbfois, nbno, nbnono, nbnoto, nbobno, nbrest
    integer(kind=8) :: nbsoob, ncoo, nmelim, nono, numno
    character(len=8), pointer :: objet_nom(:) => null()
    real(kind=8), pointer :: coordo(:) => null()
    integer(kind=8), pointer :: objet_num(:) => null()
    character(len=8), pointer :: point_nom(:) => null()
    integer(kind=8), pointer :: descobj(:) => null()
    character(len=8), pointer :: vnomobj(:) => null()
    integer(kind=8), pointer :: soob(:) => null()
    integer(kind=8), pointer :: ssob(:) => null()
    integer(kind=8), pointer :: point_num(:) => null()
    integer(kind=8), pointer :: cumul_ele(:) => null()
!-----------------------------------------------------------------------
    data cbid/' '/
!
!
!     -- ON INITIALISE LA COLLECTION QUI CONTIENT LES CORRESPONDANCES
!     -- ENTRE LES NUMEROTATIONS LOCALES DES NOEUDS GIBI ET ASTER:
    call jemarq()
    call giinco()
!
    call jeveuo('&&GILIRE.COORDO   ', 'L', vr=coordo)
    call jelira('&&GILIRE.COORDO   ', 'LONMAX', ncoo)
!
    call jeexin('&&GILIRE.NOMOBJ', iret)
    if (iret .eq. 0) then
        call utmess('F', 'PREPOST_46')
    end if
    call jeveuo('&&GILIRE.NOMOBJ', 'L', vk8=vnomobj)
    call jeveuo('&&GILIRE.DESCOBJ', 'L', vi=descobj)
    call jeveuo('&&GILIRE.CUMUL_ELE', 'L', vi=cumul_ele)
!
    call jeveuo('&&GILIRE.OBJET_NOM', 'L', vk8=objet_nom)
    call jeveuo('&&GILIRE.OBJET_NUM', 'L', vi=objet_num)
!
    call jeveuo('&&GILIRE.NUMANEW', 'L', ianema)
!
!
!     -----------------------------------------------------------------
!     --ECRITURE DU TITRE:
!     -----------------------------------------------------------------
    write (nfic, *) 'TITRE'
    write (nfic, *) '%  GIBI FECIT'
    write (nfic, *) 'FINSF'
    write (nfic, *) '%'
!
!     -----------------------------------------------------------------
!     --ECRITURE DES NOEUDS:
!     -----------------------------------------------------------------
    if (ndim .eq. 3) then
        write (nfic, *) 'COOR_3D'
    else if (ndim .eq. 2) then
        write (nfic, *) 'COOR_2D'
    else if (ndim .eq. 1) then
        write (nfic, *) 'COOR_1D'
    else
        call utmess('F', 'PREPOST_53')
    end if
    indir = .false.
    call jeexin('&&GILIRE.INDIRECT', iret)
    if (iret .ne. 0) then
        indir = .true.
        call jelira('&&GILIRE.INDIRECT', 'LONMAX', nbnoto)
        call jeveuo('&&GILIRE.INDIRECT', 'L', iaptin)
        call wkvect('&&GILIRE.NOENOM', 'V V I', nbnoto, inutri)
        call jacopo(nbnoto, 'I', iaptin, inutri)
        call uttrii(zi(inutri), nbnoto)
        nbelim = (ncoo/ndim)-nbnoto
        if (nbelim .gt. 0) then
            vali = nbelim
            call utmess('I', 'PREPOST5_19', si=vali)
        end if
    else
        nbnoto = ncoo/ndim
    end if
!
    do ino = 1, nbnoto
        if (indir) then
            nono = zi(inutri-1+ino)
        else
            nono = ino
        end if
        call codent(nono, 'G', k7bid)
        write (nfic, 1001) 'N'//k7bid, (coordo(ndim*(nono-1)+j), j= &
                                        1, ndim)
    end do
!
    write (nfic, *) 'FINSF'
    write (nfic, *) '%'
!
!     -----------------------------------------------------------------
!     --ECRITURE DES MAILLES:
!     -----------------------------------------------------------------
!
    call jelira('&&GILIRE.OBJET_NOM', 'LONMAX', nbobno)
    call wkvect('&&GILIRE.OBJTRI_NUM', 'V V I', nbobj, itrnu)
    call wkvect('&&GILIRE.ECRIGRM', 'V V L', nbobj, iecrit)
!
    call jeveuo('&&GILIRE.NUMANEW', 'L', ianema)
    call jelira('&&GILIRE.NUMANEW', 'LONUTI', itot)
!
!   CALCUL DU NB TOT D'ELEMENTS
    call wkvect('&&GILIRE.ECRMAIL', 'V V L', itot, iecrma)
!
!
    do il = 1, nbobj
        zl(iecrit+il-1) = .false.
    end do
    imb = 1
    do ima = 1, nbobno
        ii = objet_num(ima)
        if (.not. (zl(iecrit+ii-1))) then
            zi(itrnu+imb-1) = ii
            imb = imb+1
            zl(iecrit+ii-1) = .true.
        end if
    end do
!
    do i = 1, nbobno
        ii = objet_num(i)
        nbsoob = descobj(4*(ii-1)+1)
        nomobj = vnomobj(2*(ii-1)+1)
        if (nbsoob .ne. 0) then
            call jeveuo('&&GILIRE'//nomobj//'.SOUSOB', 'L', vi=soob)
            do kk = 1, nbsoob
                jj = soob(kk)
                if (.not. (zl(iecrit+jj-1))) then
                    zi(itrnu+imb-1) = jj
                    zl(iecrit+jj-1) = .true.
                    imb = imb+1
                end if
            end do
        end if
    end do
!
! ON SUPPRIME UN IMB CAR ON EN COMPTE UN DE PLUS DANS LA FIN DE BOUCLE
!
    imb = imb-1
!
! ON TRIE LA TABLE
!
    if (imb .gt. 1) then
        call uttrii(zi(itrnu), imb)
    end if
!
    icoma = 0
    nbelt = 0
    nbelc = 0
!
    do i = 1, nbobj
        trouve = .false.
        do jj = 1, imb
            ii = zi(itrnu+jj-1)
            if (i .eq. ii) then
                trouve = .true.
                goto 13
            end if
        end do
13      continue
!
        nbno = descobj(4*(i-1)+3)
        nbele = descobj(4*(i-1)+4)
        nomobj = vnomobj(2*(i-1)+1)
        tymail = vnomobj(2*(i-1)+2)
        nbelt = nbelt+nbele
        if (trouve) nbelc = nbelc+nbele
!
!        -- SI L'OBJET EST 1 OBJET SIMPLE , ON ECRIT SES MAILLES:
        if (nbele .gt. 0) then
            call giecma(nfic, trouve, nbele, nomobj, tymail, &
                        nbno, zl(iecrma), icoma)
        end if
    end do
    if (nbelc .gt. 9999999) then
        vali = nbelc
        call utmess('F', 'PREPOST6_2', si=vali)
    end if
    nmelim = nbelt-nbelc
    if (nmelim .gt. 0) then
        vali = nmelim
        call utmess('I', 'PREPOST5_20', si=vali)
    end if
!
!     -----------------------------------------------------------------
!     --ECRITURE DES GROUP_NO:
!     -----------------------------------------------------------------
!
    call jeexin('&&GILIRE.POINT_NOM', iret)
    if (iret .gt. 0) then
        call jeveuo('&&GILIRE.POINT_NOM', 'L', vk8=point_nom)
        call jeveuo('&&GILIRE.POINT_NUM', 'L', vi=point_num)
        call jelira('&&GILIRE.POINT_NOM', 'LONMAX', nbnono)
    else
        nbnono = 0
    end if
!
    do i = 1, nbnono
        nomno = point_nom(i)
        if (nomno(1:1) .eq. '#') goto 3
        numno = point_num(i)
        call codent(numno, 'G', k7bid)
        write (nfic, *) 'GROUP_NO'
        write (nfic, 1002) nomno, 'N'//k7bid
        write (nfic, *) 'FINSF'
        write (nfic, *) '%'
3       continue
    end do
!
!     -----------------------------------------------------------------
!     --ECRITURE DES GROUP_MA:
!     -----------------------------------------------------------------
!
    call jelira('&&GILIRE.OBJET_NOM', 'LONMAX', nbobno)
    do ii = 1, nbobj
        trouve = .false.
        do inu = 1, nbobno
            if (objet_num(inu) .eq. ii) then
                trouve = .true.
                nomobg = objet_nom(inu)
                if (nomobg(1:1) .eq. '#') goto 21
                write (nfic, *) 'GROUP_MA'
                write (nfic, *) '  ', nomobg
                nbsoob = descobj(4*(ii-1)+1)
                if (nbsoob .eq. 0) then
!
!           -- ON FAIT COMME SI L'OBJET SE CONTENAIT LUI-MEME:
                    nbsoob = 1
                    magoui = .true.
                else
                    magoui = .false.
                    nomobj = vnomobj(2*(ii-1)+1)
                    call jeveuo('&&GILIRE'//nomobj//'.SOUSOB', 'L', vi=ssob)
                end if
                do j = 1, nbsoob
!
!        -- L'OBJET EST 1 OBJET COMPOSE, ON ECRIT SES MAILLES:
                    if (magoui) then
                        jj = ii
                    else
                        jj = ssob(j)
                    end if
                    nomobj = vnomobj(2*(jj-1)+1)
                    nbno = descobj(4*(jj-1)+3)
                    nbele = descobj(4*(jj-1)+4)
                    nbfois = nbele/7
                    nbrest = nbele-7*nbfois
                    icok = cumul_ele(jj)
!
                    do k = 1, nbfois
                        do kk = 1, 7
                            icok = icok+1
                            call codent(zi(ianema-1+icok), 'G', k7nom(kk))
                            k8nom(kk) = 'M'//k7nom(kk)
                        end do
                        write (nfic, 1003) (k8nom(l), l=1, 7)
                    end do
!
                    do kk = 1, nbrest
                        icok = icok+1
                        call codent(zi(ianema-1+icok), 'G', k7nom(kk))
                        k8nom(kk) = 'M'//k7nom(kk)
                    end do
                    write (nfic, 1003) (k8nom(l), l=1, nbrest)
!
                end do
                write (nfic, *) 'FINSF'
                write (nfic, *) '%'
            end if
21          continue
        end do
!
    end do
!
!     -- ON ECRIT LE "FIN" FINAL ET ON REMBOBINE LE FICHIER:
!     ------------------------------------------------------
    write (nfic, *) 'FIN'
    rewind (nfic)
!
    call jedetc('V', '&&GILIRE', 1)
1001 format(1x, a8, 1x, 1pd21.14, 1x, 1pd21.14, 1x, 1pd21.14)
1002 format(1x, a8, 1x, a8)
1003 format(7(1x, a8))
!
    call jedema()
end subroutine
