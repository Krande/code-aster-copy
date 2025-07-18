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
subroutine gcour2(resu, noma, nomno, coorn, nbnoeu, &
                  trav1, trav2, trav3, fonoeu, chfond, &
                  basfon, nomfiss, connex, stok4, liss, &
                  nbre, milieu, ndimte, norfon)
    implicit none
!
! FONCTION REALISEE:
!
! 1.  POUR CHAQUE NOEUD DU FOND DE FISSURE GAMM0 ON RECUPERE
!     LE TRIPLET ( MODULE(THETA), RINF, RSUP )
!
! 2.  PUIS ON  CALCULE LA DIRECTION DES CHAMPS THETA
!     APPEL A GDIREC
!
! 3.  ENSUITE ON CALCULE LES CHAMPS THETA SUR TOUS LES NOEUDS DU
!     MAILLAGE
!
!     ------------------------------------------------------------------
! ENTREE:
!        RESU   : NOM DU CONCEPT RESULTAT
!        NOMA   : NOM DU CONCEPT MAILLAGE
!        NOMNO  : NOM DE L'OBJET CONTENANT LES NOEUDS DU MAILLAGE
!        COORN  : NOM DE L'OBJET CONTENANT LES COORDONNEES DU MAILLAGE
!        NBNOEU : NOMBRE DE NOEUDS DE GAMM0
!        FONOEU : NOMS DES NOEUDS DU FOND DE FISSURE
!        CHFOND : NOM DE L'OBJET CONTENANT LES COORDONNEES DES NOEUDS DE GAMM0
!        NOMFISS: NOM DU CONCEPT FOND_FISS
!        TRAV1  : RINF
!        TRAV2  : RSUP
!        LISS   : TYPE DE LISSAGE
!        NBRE   : DEGRE DES POLYNOMES DE LEGENDRE
!                     SINON 0
!        CONNEX: .TRUE.  : FOND DE FISSURE FERME
!                .FALSE. : FOND DE FISSURE DEBOUCHANT
! SORTIE:
!        STOK4  : DIRECTION DU CHAMP THETA
!                 LISTE DE CHAMPS_NO THETA
!        TRAV3 : MODULE(THETA)
!        MILIEU: .TRUE.  : ELEMENT QUADRATIQUE
!                .FALSE. : ELEMENT LINEAIRE
!     ------------------------------------------------------------------
!
!
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getres.h"
#include "asterc/r8maem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/gdinor.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    character(len=24) :: trav1, trav2, trav3, fonoeu
    character(len=24) :: obj3, norm, numgam, chamno, chfond, basfon
    character(len=24) :: stok4, coorn, nomno, indicg
    character(len=24) :: liss, norfon
    character(len=16) :: k16b, nomcmd
    character(len=8) :: nomfiss, resu, noma, k8b
    character(len=6) :: kiord
    character(len=1), parameter :: base = 'G'
!
    integer(kind=8) :: nbnoeu, iadrt1, iadrt2, iadrt3, itheta, ifon
    integer(kind=8) :: in2, iadrco, jmin, ielinf, iadnum, jvect
    integer(kind=8) :: iadrno, num, indic, iadrtt, nbre, nbptfd
    integer(kind=8) :: iret, numa, ndimte, iebas, iftyp, nec
!
    real(kind=8) :: xi1, yi1, zi1, xj1, yj1, zj1
    real(kind=8) :: xij, yij, zij, eps, d, tei, tej
    real(kind=8) :: xm, ym, zm, xim, yim, zim, s, dmin, smin, xn, yn, zn
    real(kind=8) :: rii, rsi, alpha, valx, valy, valz, norm2
!
    aster_logical :: milieu, connex
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, idesc, inorfon
    integer(kind=8) :: ienorm, irefe, j, jresu, k, nbel
!-----------------------------------------------------------------------
    call jemarq()
!
    call getres(k8b, k16b, nomcmd)
!
    eps = 1.d-06
    call jeveuo(trav1, 'L', iadrt1)
    call jeveuo(trav2, 'L', iadrt2)
    call jeveuo(trav3, 'E', iadrt3)
    call jeveuo(fonoeu, 'L', iadrno)
    call jeveuo(coorn, 'L', iadrco)
    call jeveuo(chfond, 'L', ifon)
!
! VERIFICATION SI MAILLAGE QUADRATIQUE OU NON
!
    call jeveuo(nomfiss//'.FOND.TYPE', 'L', iftyp)
    milieu = .false.
    if (zk8(iftyp) .eq. 'SEG3' .or. zk8(iftyp) .eq. 'NOE3') then
        milieu = .true.
    end if
!
! VERIFICATION PRESENCE NB_POINT_FOND
!
    call getvis('THETA', 'NB_POINT_FOND', iocc=1, nbval=1, nbret=nbptfd)
!
! SI PRESENCE NB_POINT_FOND ALORS NOEUDS MILIEU ABSENTS
!
    if (nbptfd .ne. 0) then
        milieu = .false.
    end if
!
! INTERDICTION D AVOIR NB_POINT_FOND AVEC LISSAGES
! LAGRANGE_NO_NO, MIXTE OU LEGENDRE
!
    if (nbptfd .ne. 0) then
        if ((liss .eq. 'LEGENDRE') .or. (liss .eq. 'MIXTE') .or. &
            (liss .eq. 'LAGRANGE_NO_NO')) then
            call utmess('F', 'RUPTURE1_73')
        end if
    end if
!
! RECUPERATION  DES NUMEROS DE NOEUDS DE GAMM0
!
    numgam = '&&COURON.NUMGAMM0'
    call wkvect(numgam, 'V V I', nbnoeu, iadnum)
    do j = 1, nbnoeu
        zi(iadnum+j-1) = char8_to_int(zk8(iadrno+j-1))
    end do
!
!  SI LEVRE_INF EST DEFINIE DANS LE CONCEPT FOND
!
    obj3 = nomfiss//'.LEVREINF.MAIL'
    call jeexin(obj3, ielinf)
!
!  SI NORMALE EST DEFINIE DANS LE CONCEPT FOND
!
    norm = nomfiss//'.NORMALE        '
    call jeexin(norm, ienorm)
!
    stok4 = '&&COURON.DIREC'
    call wkvect(stok4, 'V V R', 3*nbnoeu, in2)
    call jeexin(basfon, iebas)
!
!   DETERMINATION DE LA DIRECTION DE PROPAGATION
    if (iebas .ne. 0) then
!       CAS D'UNE FISSURE : LA DIRECTION EST A PRENDRE DANS BASEFOND
        call jeveuo(basfon, 'L', jvect)
        do i = 1, nbnoeu
            zr(in2+(i-1)*3+1-1) = zr(jvect-1+6*(i-1)+4)
            zr(in2+(i-1)*3+2-1) = zr(jvect-1+6*(i-1)+5)
            zr(in2+(i-1)*3+3-1) = zr(jvect-1+6*(i-1)+6)
        end do
    else if (ienorm .ne. 0) then
!       CAS D'UNE ENTAILLE : BASEFOND N'EXISTE PAS
        call gdinor(norm, nbnoeu, iadnum, coorn, in2)
    else
        ASSERT(.FALSE.)
    end if
!
!
    norfon = '&&NORM.STOCK'
    call wkvect(norfon, 'V V R', 3*nbnoeu, inorfon)
!
!   stockage des directions des normales au fond de fissure
    call jeexin(nomfiss//'.BASEFOND', iebas)
!
    if (iebas .ne. 0) then
!       * cas general : la base du fond de fissure est definie et on
!                       copie la normale
        call jeveuo(basfon, 'L', jvect)
        do i = 1, nbnoeu
            zr(inorfon+(i-1)*3+1-1) = zr(jvect-1+(i-1)*6+4)
            zr(inorfon+(i-1)*3+2-1) = zr(jvect-1+(i-1)*6+5)
            zr(inorfon+(i-1)*3+3-1) = zr(jvect-1+(i-1)*6+6)
        end do
    else
!   * cas particulier : la base locale n'est pas définie,
!       e.g. si l'utilisateur donne le champ de normale dans
!       DEFI_FOND_FISS
!       on copie la direction du champ theta et on aura donc
!       theta . n = 1, pour un champ theta norme
        do i = 1, nbnoeu
            zr(inorfon+(i-1)*3+1-1) = zr(in2+(i-1)*3+1-1)
            zr(inorfon+(i-1)*3+2-1) = zr(in2+(i-1)*3+2-1)
            zr(inorfon+(i-1)*3+3-1) = zr(in2+(i-1)*3+3-1)
        end do
    end if
!
! ALLOCATION D UN OBJET INDICATEUR DU CHAMP THETA SUR GAMMO
!
    call dismoi('NB_NO_MAILLA', noma, 'MAILLAGE', repi=nbel)
!
    indicg = '&&COURON.INDIC        '
!
    call wkvect(indicg, 'V V I', nbel, indic)
!
! ALLOCATION DES OBJETS POUR STOCKER LE CHAMP_NO THETA ET LA DIRECTION
! TYPE CHAM_NO ( DEPL_R) AVEC PROFIL NOEUD CONSTANT (3 DDL)
!
    if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO') .or. (liss .eq. 'MIXTE')) then
        ndimte = nbnoeu
    else
        ndimte = nbre+1
    end if
!
    call wkvect(resu, base//' V K24', ndimte+1, jresu)
!
! CREATION DES NDIMTE+1 CHAMPS_NO ET VALEUR SUR GAMMA0
!
    do k = 1, ndimte+1
        call codent(k, 'D0', kiord)
        chamno = resu(1:8)//'_CHAM'//kiord//'     '
        zk24(jresu+k-1) = chamno
        call jeexin(chamno(1:19)//'.DESC', iret)
        ASSERT(iret .ge. 0 .and. iret .le. 100)
        if (iret .eq. 0) then
            call jedetr(chamno(1:19)//'.DESC')
            call jedetr(chamno(1:19)//'.REFE')
            call jedetr(chamno(1:19)//'.VALE')
        end if
!  .DESC
        chamno(20:24) = '.DESC'
        call dismoi('NB_EC', 'DEPL_R', 'GRANDEUR', repi=nec)
        call wkvect(chamno, base//' V I', 2+nec, idesc)
!
        call jeecra(chamno, 'DOCU', cval='CHNO')
        call jenonu(jexnom('&CATA.GD.NOMGD', 'DEPL_R'), numa)
        zi(idesc+1-1) = numa
        zi(idesc+2-1) = -3
        zi(idesc+3-1) = 14
!  .REFE
        chamno(20:24) = '.REFE'
        call wkvect(chamno, base//' V K24', 4, irefe)
        zk24(irefe+1-1) = noma//'                '
!  .VALE
        chamno(20:24) = '.VALE'
        call wkvect(chamno, base//' V R', 3*nbel, itheta)
!
        if (k .ne. (ndimte+1)) then
            if ((liss .eq. 'LAGRANGE') .or. (liss .eq. 'LAGRANGE_NO_NO') .or. &
                (liss .eq. 'MIXTE')) then
                if (nbptfd .eq. 0) then
                    do i = 1, nbnoeu
                        num = zi(iadnum+i-1)
                        zr(itheta+(num-1)*3+1-1) = 0.d0
                        zr(itheta+(num-1)*3+2-1) = 0.d0
                        zr(itheta+(num-1)*3+3-1) = 0.d0
                        zi(indic+num-1) = 1
                    end do
                    num = zi(iadnum+k-1)
                    iadrtt = iadrt3+(k-1)*nbnoeu+k-1
                    zr(iadrtt) = 1.d0
                    zr(itheta+(num-1)*3+1-1) = zr(iadrtt)*zr(in2+(k-1)*3+1-1)
                    zr(itheta+(num-1)*3+2-1) = zr(iadrtt)*zr(in2+(k-1)*3+2-1)
                    zr(itheta+(num-1)*3+3-1) = zr(iadrtt)*zr(in2+(k-1)*3+3-1)
                    if (connex .and. (k .eq. 1)) then
                        num = zi(iadnum+ndimte-1)
                        iadrtt = iadrt3+(k-1)*nbnoeu+ndimte-1
                        zr(iadrtt) = 1.d0
                        zr(itheta+(num-1)*3+1-1) = zr(iadrtt)*zr(in2+(ndimte-1)*3+1-1)
                        zr(itheta+(num-1)*3+2-1) = zr(iadrtt)*zr(in2+(ndimte-1)*3+2-1)
                        zr(itheta+(num-1)*3+3-1) = zr(iadrtt)*zr(in2+(ndimte-1)*3+3-1)
                    end if
                    if (connex .and. (k .eq. ndimte)) then
                        num = zi(iadnum+1-1)
                        iadrtt = iadrt3+(k-1)*nbnoeu+1-1
                        zr(iadrtt) = 1.d0
                        zr(itheta+(num-1)*3+1-1) = zr(iadrtt)*zr(in2+(1-1)*3+1-1)
                        zr(itheta+(num-1)*3+2-1) = zr(iadrtt)*zr(in2+(1-1)*3+2-1)
                        zr(itheta+(num-1)*3+3-1) = zr(iadrtt)*zr(in2+(1-1)*3+3-1)
                    end if
                else
                    zr(iadrt3-1+(k-1)*nbnoeu+k) = 1.d0
                    if (connex .and. (k .eq. 1)) then
                        iadrtt = iadrt3+(k-1)*nbnoeu+ndimte-1
                        zr(iadrtt) = 1.d0
                    end if
                    if (connex .and. (k .eq. ndimte)) then
                        iadrtt = iadrt3+(k-1)*nbnoeu+1-1
                        zr(iadrtt) = 1.d0
                    end if
                end if
            else
                if (nbptfd .eq. 0) then
                    do i = 1, nbnoeu
                        num = zi(iadnum+i-1)
                        iadrtt = iadrt3+(k-1)*nbnoeu+i-1
                        zr(itheta+(num-1)*3+1-1) = zr(iadrtt)*zr(in2+(i-1)*3+1-1)
                        zr(itheta+(num-1)*3+2-1) = zr(iadrtt)*zr(in2+(i-1)*3+2-1)
                        zr(itheta+(num-1)*3+3-1) = zr(iadrtt)*zr(in2+(i-1)*3+3-1)
                        zi(indic+num-1) = 1
                    end do
                else
                    do i = 1, nbnoeu
                        zr(iadrt3-1+(k-1)*nbnoeu+i) = 1.d0
                    end do
                end if
            end if
        else
!     STOCKAGE DE LA DIRECTION DU CHAMPS THETA SUR LE FOND DE FISSURE
            if (nbptfd .eq. 0) then
                do i = 1, nbnoeu
                    num = zi(iadnum+i-1)
                    zr(itheta+(num-1)*3+1-1) = zr(in2+(i-1)*3+1-1)
                    zr(itheta+(num-1)*3+2-1) = zr(in2+(i-1)*3+2-1)
                    zr(itheta+(num-1)*3+3-1) = zr(in2+(i-1)*3+3-1)
                end do
            end if
        end if
    end do
!
!         BOUCLE SUR LES NOEUDS M COURANTS DU MAILLAGE SANS GAMMO
!         POUR CALCULER PROJ(M)=N
!
    do i = 1, nbel
        if ((zi(indic+i-1) .ne. 1) .or. (nbptfd .ne. 0)) then
            zr(itheta+(i-1)*3+1-1) = 0.d0
            zr(itheta+(i-1)*3+2-1) = 0.d0
            zr(itheta+(i-1)*3+3-1) = 0.d0
            xm = zr(iadrco+(i-1)*3+1-1)
            ym = zr(iadrco+(i-1)*3+2-1)
            zm = zr(iadrco+(i-1)*3+3-1)
            dmin = r8maem()
            jmin = 0
            smin = 0.d0
            do j = 1, nbnoeu-1
                xi1 = zr(ifon-1+4*(j-1)+1)
                yi1 = zr(ifon-1+4*(j-1)+2)
                zi1 = zr(ifon-1+4*(j-1)+3)
                xj1 = zr(ifon-1+4*(j-1+1)+1)
                yj1 = zr(ifon-1+4*(j-1+1)+2)
                zj1 = zr(ifon-1+4*(j-1+1)+3)
                xij = xj1-xi1
                yij = yj1-yi1
                zij = zj1-zi1
                xim = xm-xi1
                yim = ym-yi1
                zim = zm-zi1
                s = xij*xim+yij*yim+zij*zim
                norm2 = xij*xij+yij*yij+zij*zij
                s = s/norm2
                if ((s-1) .ge. eps) then
                    s = 1.d0
                end if
                if (s .le. eps) then
                    s = 0.d0
                end if
                xn = s*xij+xi1
                yn = s*yij+yi1
                zn = s*zij+zi1
                d = sqrt((xn-xm)*(xn-xm)+(yn-ym)*(yn-ym)+(zn-zm)*(zn-zm))
                if (d .lt. (dmin*(1-abs(r8prem())*1.d04))) then
                    dmin = d
                    jmin = j
                    smin = s
                end if
            end do
            rii = (1-smin)*zr(iadrt1+jmin-1)+smin*zr(iadrt1+jmin+1-1)
            rsi = (1-smin)*zr(iadrt2+jmin-1)+smin*zr(iadrt2+jmin+1-1)
            alpha = (dmin-rii)/(rsi-rii)
            do k = 1, ndimte+1
                call codent(k, 'D0', kiord)
                chamno = resu(1:8)//'_CHAM'//kiord//'     '
                chamno(20:24) = '.VALE'
                call jeveuo(chamno, 'E', itheta)
                if (k .ne. (ndimte+1)) then
                    iadrtt = iadrt3+(k-1)*nbnoeu+jmin-1
                    tei = zr(iadrtt)
                    tej = zr(iadrtt+1)
                    valx = (1-smin)*zr(in2+(jmin-1)*3+1-1)*tei
                    valx = valx+smin*zr(in2+(jmin+1-1)*3+1-1)*tej
                    valy = (1-smin)*zr(in2+(jmin-1)*3+2-1)*tei
                    valy = valy+smin*zr(in2+(jmin+1-1)*3+2-1)*tej
                    valz = (1-smin)*zr(in2+(jmin-1)*3+3-1)*tei
                    valz = valz+smin*zr(in2+(jmin+1-1)*3+3-1)*tej
!
                    if ((abs(alpha) .le. eps) .or. (alpha .lt. 0)) then
                        zr(itheta+(i-1)*3+1-1) = valx
                        zr(itheta+(i-1)*3+2-1) = valy
                        zr(itheta+(i-1)*3+3-1) = valz
                    else if ((abs(alpha-1) .le. eps) .or. ((alpha-1) .gt. 0)) &
                        then
                        zr(itheta+(i-1)*3+1-1) = 0.d0
                        zr(itheta+(i-1)*3+2-1) = 0.d0
                        zr(itheta+(i-1)*3+3-1) = 0.d0
                    else
                        zr(itheta+(i-1)*3+1-1) = (1-alpha)*valx
                        zr(itheta+(i-1)*3+2-1) = (1-alpha)*valy
                        zr(itheta+(i-1)*3+3-1) = (1-alpha)*valz
                    end if
                else
                    zr(itheta+(i-1)*3+1-1) = 0.d0
                    zr(itheta+(i-1)*3+2-1) = 0.d0
                    zr(itheta+(i-1)*3+3-1) = 0.d0
                end if
            end do
        end if
    end do
!
! DESTRUCTION D'OBJETS DE TRAVAIL
!
    call jedetr(indicg)
    call jedetr(numgam)
!
    call jedema()
!
!
end subroutine
