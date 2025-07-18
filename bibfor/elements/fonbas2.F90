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

subroutine fonbas2(noma, basnof, typm, fonoeu, coorfond, nbnoff, absfon, &
                   basloc, abscur, lnno, ltno)
!
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnscre.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/wkvect.h"

    integer(kind=8)           :: nbnoff
    character(len=8)  :: noma, typm
    character(len=19) :: basnof, basloc, lnno, ltno
    character(len=24) :: fonoeu, absfon, abscur, coorfond
!
! FONCTION REALISEE:
!
!     CALCUL :
!        - DE LA BASE LOCALE EN CHAQUE NOEUD DU MAILLAGE
!        - DES LEVEL-SETS    EN CHAQUE NOEUD DU MAILLAGE
!
!
!     ENTREES:
!        NOMA   : NOM DU MAILLAGE
!        BASNOF : BASE AUX NOEUDS DU FOND DE FISSURE
!        FONOEU : NOM DES NOEUDS DU FOND DE FISSURE
!        NBNOFF : NOMBRE DE NOEUDS AU FOND DE FISSURE
!        TYPM   : TYPE DE FOND DE FISSURE : LIN OU QUAD
!        ABSFON : ABSCISSE CURVILGNE DE CHAQUE NOEUD DU FOND DE FISSURE
!        COORFOND : COORDONNE DES NOEUDS DU FOND
!     SORTIES:
!        BASLOC : BASE LOCALE EN CHAQUE NOEUD DU MAILLAGE
!        ABSCUR : ABSCISSE CURVILIGNE DU PROJETE SUR LE FOND DE FISSURE
!                 DE CHAQUE NOEUD DU MAILLAGE
!        LTNO   : LEVEL-SETS TANGENTS EN CHAQUE NOEUD DU MAILLAGE
!        LNNO   : LEVLE-SETS NORMAUX EN CHAQUE NOEUD DU MAILLAGE
!-----------------------------------------------------------------------
!
    integer(kind=8) :: ibid, indica, indicb, ina, inb, ino, jnoe
    integer(kind=8) :: iseg, jbas, jabsf, jabscur, coorfd
    integer(kind=8) :: jgsl, jlnsl, jltsl
    integer(kind=8) :: k, nbno, ndim, nseg
    real(kind=8) :: d, dmin, eps, norm2, s, sn, xln, xlt
    real(kind=8) :: xa, ya, za, xb, yb, zb, xm, ym, zm
    real(kind=8) :: xab, yab, zab, xam, yam, zam, xnm, ynm, znm
    real(kind=8) :: n(3), nm(3), vdira(3), vnora(3), vdirb(3), vnorb(3)
    real(kind=8) :: vdirn(3), vnorn(3)
    character(len=8) :: licmp(9)
    character(len=16) :: casfon
    character(len=19) :: cnsbas, cnsln, cnslt
    real(kind=8), pointer :: gsv(:) => null()
    real(kind=8), pointer :: lnsv(:) => null()
    real(kind=8), pointer :: ltsv(:) => null()
    real(kind=8), pointer :: vale(:) => null()
!
    data licmp/'X1', 'X2', 'X3',&
     &             'X4', 'X5', 'X6',&
     &             'X7', 'X8', 'X9'/
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
!     ------------------------------------------------------------------
!     INITIALISATIONS
!     ------------------------------------------------------------------
!
!     RECUPERATION DES INFORMATIONS RELATIVES AU MAILLAGE
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbno)
    call jeveuo(noma//'.COORDO    .VALE', 'L', vr=vale)
!
!     NSEG : NOMBRE DE "SEGMENTS" DU FOND A TRAITER
    if (ndim .eq. 2) then
        nseg = 1
        casfon = ' '
    else if (ndim .eq. 3) then
        casfon = 'LINEAIRE'
        nseg = nbnoff-1
!       CAS QUADRATIQUE
        if (typm .eq. 'SEG3') then
            casfon = 'QUADRATIQUE'
            nseg = (nbnoff-1)/2
        end if
    end if
!
!     INITIALISATION DES CHAMPS SIMPLES DES LEVEL-SETS
    cnslt = '&&FONBAS2.CNSLT'
    cnsln = '&&FONBAS2.CNSLN'
    call cnscre(noma, 'NEUT_R', 1, 'X1', 'V', &
                cnslt)
    call cnscre(noma, 'NEUT_R', 1, 'X1', 'V', &
                cnsln)
!
    call jeveuo(cnslt//'.CNSV', 'E', vr=ltsv)
    call jeveuo(cnslt//'.CNSL', 'E', jltsl)
    call jeveuo(cnsln//'.CNSV', 'E', vr=lnsv)
    call jeveuo(cnsln//'.CNSL', 'E', jlnsl)
!
!     INITIALISATION DU CHAMP SIMPLE DE LA BASE LOCALE
    cnsbas = '&&FONBAS2.CNSBAS'
    call cnscre(noma, 'NEUT_R', ndim*3, licmp, 'V', &
                cnsbas)
    call jeveuo(cnsbas//'.CNSV', 'E', vr=gsv)
    call jeveuo(cnsbas//'.CNSL', 'E', jgsl)
!
    call jeveuo(basnof, 'L', jbas)
!
!   RECUPERATION DES NOMS DES NOEUDS DU FOND DE FISSURE
    call jeveuo(fonoeu, 'L', jnoe)
!
!   RECUPERATION DE L'ABSCISSE CURVILIGNE DES NOEUDS DU FOND DE FISSURE
    call jeveuo(absfon, 'L', jabsf)
!
!   CREATION DU VECTEUR CONTENANT L'ABSCISSE CURVILIGNE DU PROJETE
!   SUR LE FOND DE FISSURE DE CHAQUE NOEUD DU MAILLAGE
    call wkvect(abscur, 'G V R', nbno, jabscur)
!
!   COORDONNE DES NOEUDS DU FOND
    call jeveuo(coorfond, 'L', coorfd)
!     ------------------------------------------------------------------
!                BOUCLE SUR LES NOEUDS DU MAILLAGE
!     ------------------------------------------------------------------
!
    eps = 1.d-12
    do ino = 1, nbno
!
!       COORD DU NOEUD M DU MAILLAGE
        xm = vale(1+(ino-1)*3+1-1)
        ym = vale(1+(ino-1)*3+2-1)
        zm = vale(1+(ino-1)*3+3-1)
!
!       CAS 2D : LE PROJETE EST TRIVIAL !
        if (ndim .eq. 2) then
!
!         COORD PT N
            n(1) = zr(coorfd)
            n(2) = zr(coorfd+1)
            n(3) = 0.d0
!
!         VECTEUR NM
            nm(1) = xm-n(1)
            nm(2) = ym-n(2)
!
!         STOCKAGE DES VECTEURS DE LA BASE
            do k = 1, ndim
!               STOCKAGE COORDONNEES DU PROJETE
                gsv(3*ndim*(ino-1)+k) = n(k)
                zl(jgsl-1+3*ndim*(ino-1)+k) = .true.
!               STOCKAGE VECTEUR DIRECTION
                gsv(3*ndim*(ino-1)+k+2) = zr(jbas-1+k+ndim)
                zl(jgsl-1+3*ndim*(ino-1)+k+2) = .true.
!               STOCKAGE VECTEUR NORMAL
                gsv(3*ndim*(ino-1)+k+4) = zr(jbas-1+k)
                zl(jgsl-1+3*ndim*(ino-1)+k+4) = .true.
            end do
!
!         STOCKAGE DES LEVEL-SETS
            lnsv((ino-1)+1) = nm(1)*zr(jbas-1+1)+nm(2)*zr(jbas-1+ &
                                                          2)
            ltsv((ino-1)+1) = nm(1)*zr(jbas-1+3)+nm(2)*zr(jbas-1+ &
                                                          4)
            zl(jlnsl-1+(ino-1)+1) = .true.
            zl(jltsl-1+(ino-1)+1) = .true.
!
!       CAS 3D : RECHERCHE DU PROJETE PUIS STOCKAGE DES VECTEURS
        else if (ndim .eq. 3) then
!
!         RECHERCHE DU PROJETE DE INO SUR LE FOND DE FISSURE
!         --------------------------------------------------
            dmin = r8maem()
!
!         BOUCLE SUR LES "SEGMENTS" DU FOND DE FISSURE
            do iseg = 1, nseg
!
                if (casfon .eq. 'LINEAIRE') then
                    ina = iseg
                    inb = iseg+1
                else if (casfon .eq. 'QUADRATIQUE') then
                    ina = 2*iseg-1
                    inb = 2*iseg+1
                end if
!
!           COORD DES POINTS A ET B, EXTREMITES DU SEGMENT ISEG
                xa = zr(coorfd-1+3*(ina-1)+1)
                ya = zr(coorfd-1+3*(ina-1)+2)
                za = zr(coorfd-1+3*(ina-1)+3)
                xb = zr(coorfd-1+3*(inb-1)+1)
                yb = zr(coorfd-1+3*(inb-1)+2)
                zb = zr(coorfd-1+3*(inb-1)+3)
!
!           VECTEUR AB ET AM
                xab = xb-xa
                yab = yb-ya
                zab = zb-za
                xam = xm-xa
                yam = ym-ya
                zam = zm-za
!
!           PARAM S (PRODUIT SCALAIRE...)
                s = xab*xam+yab*yam+zab*zam
                norm2 = xab*xab+yab*yab+zab*zab
                s = s/norm2
!
!           SI N EN DEHORS DU SEGMENT AB
                if ((s-1) .ge. eps) s = 1.d0
                if (s .le. eps) s = 0.d0
!
!           COORD DU PROJETE DE M SUR ISEG: N
                xnm = xm-(s*xab+xa)
                ynm = ym-(s*yab+ya)
                znm = zm-(s*zab+za)
!
!           DISTANCE MN
                d = sqrt(xnm*xnm+ynm*ynm+znm*znm)
!
                if (d .lt. (dmin*(1-abs(r8prem())*1.d04))) then
                    dmin = d
                    sn = s
                    indica = ina
                    indicb = inb
!
                    n(1) = s*xab+xa
                    n(2) = s*yab+ya
                    n(3) = s*zab+za
                end if
!
            end do
!           ABSCISSE CURVILIGNE DU NOEUD N SUR LE FRONT DE FISSURE
            zr(jabscur-1+ino) = (1-sn)*zr(jabsf-1+indica)+sn*zr(jabsf-1+indicb)
!
!         CALCUL DES VECTEURS DE LA BASE LOCALE AU POINT PROJETE
!         ------------------------------------------------------
!
            nm(1) = xm-n(1)
            nm(2) = ym-n(2)
            nm(3) = zm-n(3)
!
            do k = 1, ndim
!
                vnora(k) = zr(jbas-1+6*(indica-1)+k)
                vdira(k) = zr(jbas-1+6*(indica-1)+k+ndim)
                vnorb(k) = zr(jbas-1+6*(indicb-1)+k)
                vdirb(k) = zr(jbas-1+6*(indicb-1)+k+ndim)
                vnorn(k) = sn*vnorb(k)+(1-sn)*vnora(k)
                vdirn(k) = sn*vdirb(k)+(1-sn)*vdira(k)
!
!               STOCKAGE COORDONNEES DU PROJETE
                gsv(3*ndim*(ino-1)+k) = n(k)
                zl(jgsl-1+3*ndim*(ino-1)+k) = .true.
!               STOCKAGE VECTEUR DIRECTION
                gsv(3*ndim*(ino-1)+k+3) = vdirn(k)
                zl(jgsl-1+3*ndim*(ino-1)+k+3) = .true.
!               STOCKAGE VECTEUR NORMAL
                gsv(3*ndim*(ino-1)+k+6) = vnorn(k)
                zl(jgsl-1+3*ndim*(ino-1)+k+6) = .true.
!
            end do
!
!         STOCKAGE DES LEVEL-SETS
            xln = nm(1)*vnorn(1)+nm(2)*vnorn(2)+nm(3)*vnorn(3)
            xlt = nm(1)*vdirn(1)+nm(2)*vdirn(2)+nm(3)*vdirn(3)
            lnsv((ino-1)+1) = xln
            ltsv((ino-1)+1) = xlt
!
            zl(jlnsl-1+(ino-1)+1) = .true.
            zl(jltsl-1+(ino-1)+1) = .true.
!
!       CAS NI 2D NI 3D
        else
!
            ASSERT(.false.)
!
        end if
!
    end do
!
!
! --- CREATION DES CHAM_NO
!
!     ENREGISTREMENT DE .LTNO, .LNNO ET .BASLOC DANS LA SD FOND_FISS
    call cnscno(cnslt, ltno(1:13)//'.NUMEQ', 'NON', 'G', ltno, &
                'F', ibid)
    call cnscno(cnsln, ltno(1:13)//'.NUMEQ', 'NON', 'G', lnno, &
                'F', ibid)
    call cnscno(cnsbas, basloc(1:13)//'.NUMEQ', 'NON', 'G', basloc, &
                'F', ibid)
!
!
!     MENAGE
    call detrsd('CHAM_NO_S', cnsln)
    call detrsd('CHAM_NO_S', cnslt)
    call detrsd('CHAM_NO_S', cnsbas)
!
    call jedema()
end subroutine
