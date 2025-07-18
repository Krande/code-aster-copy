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

subroutine fonnor2(resu, noma, cnxinv, typm, basnof)
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dis2no.h"
#include "asterfort/dismoi.h"
#include "asterfort/fonext.h"
#include "asterfort/fonno1.h"
#include "asterfort/fonno2.h"
#include "asterfort/fonno3.h"
#include "asterfort/fonno4.h"
#include "asterfort/fonno5.h"
#include "asterfort/fonno52.h"
#include "asterfort/fonno62.h"
#include "asterfort/fonno7.h"
#include "asterfort/fonno8.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/char8_to_int.h"
!
    character(len=8) :: resu, noma, typm
    character(len=19) :: cnxinv, basnof
! FONCTION REALISEE:
!
!     CALCUL EN 2D ET 3D DE LA BASE LOCALE POUR DEFI_FOND_FISS
!
!
!     ENTREES:
!        RESU       : NOM DU CONCEPT RESULTAT DE L'OPERATEUR
!        NOMA       : NOM DU MAILLAGE
!        CNXINV     : CONNECTIVITE INVERSE
!        TYPM       : TYPE DE FOND : LIN OU QUAD
!        BASNOF     : BASE LOCALE EN FOND DE FISSURE
!-----------------------------------------------------------------------
!
    integer(kind=8) :: j, jnoe1, jbasno, jbasse, jtail, k
    integer(kind=8) :: jborl, jdirol, jnvdir, jnor
    integer(kind=8) :: i, ina, inb, iseg, iret, nbnose, nbnoff, inc, inor
    integer(kind=8) :: na, nb, nret, ndim, nbnoel, nseg, nbmax, nbmac, inoext
    integer(kind=8) :: indic(4), noe(4, 4), indr(2), tablev(2), inoseg, nblev
    integer(kind=8), pointer :: connex(:) => null()
    real(kind=8) :: vdir(2, 3), vnor(2, 3), norme, vecdir(3), hmax, hmaxpr
    real(kind=8) :: vect(3), sens, dist, disttemp
    real(kind=8), pointer :: geom(:) => null()
    character(len=6) :: tyfond
    character(len=8) :: noeua, nompro
    character(len=8), pointer :: mail(:) => null()
    character(len=16) :: casfon
    character(len=19) :: basseg, macofo
    parameter(nompro='FONNOR2')
!     -----------------------------------------------------------------
!
    call jemarq()
!
!     ------------------------------------------------------------------
!     INITIALISATIONS
!     ------------------------------------------------------------------
!
    indr = 0
    nbnoel = 0
    dist = 0
    disttemp = 0
    noe = 0.d0
!
!     RECUPERATION DES INFORMATIONS RELATIVES AU MAILLAGE
!
!
!     RECUPERATION DU CONCEPT DU MAILLAGE
    call getvid(' ', 'MAILLAGE', scal=noma, nbret=nret)
!
!     RECUPERATION DU NOMBRE DE NOEUDS DU MAILLAGE
    call dismoi('DIM_GEOM', noma, 'MAILLAGE', repi=ndim)
!
!     RECUPERATION DU TYPE DU FOND DE FISSURE OUVERT OU FERME
    call getvtx('FOND_FISS', 'TYPE_FOND', iocc=1, scal=tyfond, nbret=iret)
!
!     RECUPERATION DES NOEUDS DU FOND DE FISSURE
!
    call jeexin(resu//'.FOND.NOEU', iret)
    if (iret .ne. 0) then
!       RECUPERATION DE L'ADRESSE DES NOEUDS DE FOND DE FISSURE
        call jeveuo(resu//'.FOND.NOEU', 'L', jnoe1)
!       RECUPERATION DU NOMBRE DE NOEUD
        call jelira(resu//'.FOND.NOEU', 'LONUTI', nbnoff)
    else
        ASSERT(.FALSE.)
    end if
!
!
!     VERIFICATION DE LA PRESENCE DE NORMALE
    call jeexin(resu//'.NORMALE', inor)
!
!     INITIALISATION DU SENS DU VECTEUR NORMAL
    sens = 1.d0
!
!     RECUPERATION DU TYPE DE MAILLE EN FOND DE FISSURE EN 3D
    if (ndim .eq. 3) then
!       SUIVANT LE CAS QUADRATIQUE/LINEAIRE EN 3D DEUX MAILLES SONT
!       CONNECTEES SI ELLES ONT AU MOINS NBMAX NOEUDS EN COMMUN
!       NBNOSE : NOMBRE DE NOEUDS PAR "SEGMENT" DE FOND DE FISSURE
        if (typm .eq. 'SEG2') then
            casfon = 'LINEAIRE'
            nbnose = 2
            nbmax = 3
        else if (typm .eq. 'SEG3') then
            casfon = 'QUADRATIQUE'
            nbnose = 3
            nbmax = 6
        else
            nbnose = 2
            nbmax = 3
!
        end if
    else
        nbnose = 2
        nbmax = 3
    end if
    if (ndim .eq. 2) casfon = '2D'

!
!     RECUPERATION DU VECTEUR .NORMALE S'IL EXISTE ET REMPLISSAGE DE
!     VNOR. S'IL N'EXISTE PAS, VNOR SERA REMPLI DANS FONNO5
    if (inor .ne. 0) then
        call jeveuo(resu//'.NORMALE', 'L', jnor)
        do i = 1, 3
!     DEUX VECTEUR NOMALE IDENTIQUES (UN PAR LEVRE)
            vnor(1, i) = zr(jnor-1+i)
            vnor(2, i) = zr(jnor-1+i)
        end do
    end if
!
!
!     ALLOCATION DU VECTEUR DES BASES LOCALES PAR NOEUD DU FOND  :
!           - VECTEUR DIRECTION DE PROPA
!           - VECTEUR NORMAL (A LA SURFACE)
    call wkvect(basnof, 'G V R', 2*ndim*nbnoff, jbasno)
!
!     ALLOCATION DU VECTEUR DES TAILLES DE MAILLES MAX PAR NOEUD DU FOND
    call wkvect(resu//'.FOND.TAILLE_R', 'G V R', nbnoff, jtail)
!
!
!     NSEG : NOMBRE DE "SEGMENTS" DU FOND A TRAITER
    if (ndim .eq. 2) then
        ASSERT(nbnoff .eq. 1)
        nseg = 1
    else if (ndim .eq. 3) then
        if (.not. (nbnoff .gt. 1)) then
            call utmess('F', 'RUPTURE0_92', sk=noma)
        end if
        if (casfon .eq. 'LINEAIRE') nseg = nbnoff-1
        if (casfon .eq. 'QUADRATIQUE') nseg = (nbnoff-1)/2
    end if
!
!     VECTEUR TEMPORAIRE DES BASES LOCALES PAR SEGMENT DU FOND
    basseg = '&&'//nompro//'.BASSEG'
    call wkvect(basseg, 'V V R', 2*ndim*nseg, jbasse)
!
!     VECTEUR PERMETTANT DE SAVOIR SI LE VECTEUR DE DIRECTION DE
!     PROPAGATION (VDIR) A ETE RECALCULE OU NON AUX POINTS
!     EXTREMITES DE FONFIS
    call wkvect('&&FONNOR2.LBORD', 'V V L', nbnoff, jborl)
!
!     VECTEUR CONTENANT LES VDIR INITIAUX (CAD SANS MODIFICATION
!     DES VECTEURS AUX POINTS EXTREMITES DE FONFIS)
    call wkvect('&&FONNOR2.VDIROL', 'V V R', 3*nbnoff, jdirol)
!
!     VECTEUR CONTENANT 0 OU 1 AUX POINTS EXTREMITES DE FONFIS:
!     0: LE PRODUIT SCALAIRE ENTRE LA NORMALE A LA FACE DE BORD ET
!        LE VDIR INITIAL ESI INFERIEUR A 0
!     1: LE PRODUIT SCALAIRE EST SUPERIEUR OU EGAL A 0
    call wkvect('&&FONNOR2.NVDIR', 'V V I', nbnoff, jnvdir)
!
    do i = 1, nbnoff
        zl(jborl-1+i) = .false.
    end do
!
!     ------------------------------------------------------------------
!     BOUCLE SUR LES "SEGMENTS" DU FOND DE FISSURE
!     ------------------------------------------------------------------
!
    do iseg = 1, nseg
!
!       INDICES DES NOEUDS DU SEGMENT :
!       NOEUDS SOMMETS (INA ET INB), NOEUD MILIEU (INC)
        if (casfon .eq. '2D') then
            ina = iseg
        else if (casfon .eq. 'LINEAIRE') then
            ina = iseg
            inb = iseg+1
        else if (casfon .eq. 'QUADRATIQUE') then
            ina = 2*iseg-1
            inb = 2*iseg+1
            inc = 2*iseg
        end if
!
!       NUMEROS (ABSOLUS) DU PREMIER NOEUDS SOMMETS DU SEGMENT : NA
        noeua = zk8(jnoe1-1+ina)
        na = char8_to_int(noeua)
        if (ndim .eq. 3) then
!       EN 3D : NB EST LE NUMERO (ABSOLU) DU DEUXIEME NOEUD SOMMETS DU SEGMENT
            nb = char8_to_int(zk8(jnoe1-1+inb))
            if (iseg .eq. 1) then
                inoext = na
                inoseg = nb
            end if
            if (iseg .eq. nseg) then
                inoext = nb
                inoseg = na
            end if
        else if (ndim .eq. 2) then
!       EN 2D : NB EST LE NUMERO (ABSOLU) D'UN NOEUD DE L'EXTREMITE DE LA
!       LEVRE SUPERIEURE QUE L'ON CHERCHE COMME ETANT UN NOEUD "ELOIGNE"
!       DU FOND

!       ON RECUPERE LES LEVRES DE LA FISSURE ET LES COORDONNEES DES NOEUDS
            call jeveuo(resu//'.LEVRESUP.MAIL', 'L', vk8=mail)
            call jelira(resu//'.LEVRESUP.MAIL', 'LONUTI', nblev)
            call jeveuo(noma//'.COORDO    .VALE', 'L', vr=geom)
            do i = 1, nblev
!               POUR CHAQUE SEGMENT, ON RECUPERE LA CONNEXITE
                iret = char8_to_int(mail(i))
                call jeveuo(jexnum(noma//'.CONNEX', iret), 'L', vi=connex)

!               ON TESTE UNIQUEMENT UN NOEUD PAR SEGMENT (SUFFISANT)
                disttemp = dis2no(geom, na, connex(1))
                if (disttemp .ge. dist) then
                    nb = connex(1)
                    dist = disttemp
                end if

            end do

        else
            ASSERT(.FALSE.)
        end if
!
!
!       1) RECUP DES NUMEROS DES MAILLES CONNECTEES AU SEGMENT DU FOND
!          -> REMPLISSAGE DE MACOFO
!       --------------------------------------------------------------
!
!       VECTEUR DES MAILLES CONNECTEES AU SEGMENT DU FOND
        macofo = '&&'//'NOMPRO'//'.MACOFOND'
        call fonno1(noma, cnxinv, ndim, na, nb, &
                    nbmac, macofo)
!
!
!       2) PARMI LES MAILLES CONNECTEES AU SEGMENT DU FOND, FILTRAGE DES
!          MAILLES CONNECTEES A 1 LEVRE (CAD AYANT UNE FACE LIBRE)
!          -> REMPLISSAGE DE TABLEV
!       ----------------------------------------------------------------
!
        call fonno2(macofo, noma, nbmac, nbnoff, nbnose, &
                    nbmax, noeua, tablev)
!
!
!
!       3) RECUP DES FACES / ARETES CONNECTEES AU FOND
!          POUR CHACUNE DES 2 MAILLES
!          -> REMPLISSAGE DE NOE
!       ----------------------------------------------------
!
        if (inor .eq. 0) then
!           CETTE OPERATION N'EST PAS FAITE SI ON N'A PAS LA NORMALE
            call fonno3(noma, tablev, ndim, na, nb, &
                        noe)
        end if
!
!
!       4) FILTRE DES FACES LIBRES
!          -> REMPLISSAGE DE INDIC
!       ----------------------------------------------------
!
        call fonno4(ndim, macofo, noma, nbmac, tablev, &
                    noe, nbnoff, indic)
        call jedetr(macofo)
!
!       5) CALCUL DES VECTEURS DE LA BASE LOCALE :
!          -> REMPLISSAGE DE VDIR ET VNOR
!            VNOR : VECTEUR NORMAL A LA SURFACE DE LA FISSURE
!            VDIR : VECTEUR DANS LA DIRECTION DE PROPAGATION
!        RQ : CHACUN CONTIENT EN FAIT 2 VECTEURS (UN PAR LEVRE)
!       --------------------------------------------------------
!
!
!       ON TEST L'EXISTENCE DU MOT CLE NORMALE CORRESPONDANT
!       AU CAS LEVRES DECOLLEES
!
        if (inor .eq. 0) then
            call fonno5(noma, indic, noe, na, nb, &
                        ndim, nbnoel, indr, vnor, vdir)
        else
            call fonno52(noma, na, nb, ndim, vnor, vdir)
        end if
!
!
!       6) DETERMINATION DU VRAI VECTEUR ET BASE PAR SEGMENT
!          -> REMPLISSAGE DE BASSEG
!       ----------------------------------------------------
!
        if (iseg .eq. 1) then
            call fonno8(resu, noma, tablev, vect)
        end if
!

        call fonno62(resu, noma, ndim, &
                     iseg, noe, indr, nbnoel, &
                     vnor, vdir, basseg, vect, sens)

!
!
!
!       7) EN 3D : BASE LOCALE : PASSAGE SEGMENTS -> NOEUDS
!          -> REMPLISSAGE DE BASNOF
!       ---------------------------------------------------
!
        if (ndim .eq. 3) then
!
!         MOYENNE POUR LES NOEUDS SOMMETS INA ET INB
!         DIRECT POUR LE NOEUD MILIEU INC
            do j = 1, 6
!
                zr(jbasno-1+6*(ina-1)+j) = (zr(jbasno-1+6*(ina-1)+j) &
                                            +zr(jbasse-1+6*(iseg-1)+j))/2.d0
!
                zr(jbasno-1+6*(inb-1)+j) = (zr(jbasno-1+6*(inb-1)+j) &
                                            +zr(jbasse-1+6*(iseg-1)+j))/2.d0
!
                if (casfon .eq. 'QUADRATIQUE') zr(jbasno-1+6*(inc-1)+j) = zr( &
                                                                          jbasse-1+6*(iseg-1)+j)
!
            end do
!
!         NORMALISATIONS
!         ZR(JBASNO-1+6*(INC-1)+J) DEJA NORMALISE
            call normev(zr(jbasno-1+6*(ina-1)+1), norme)
            call normev(zr(jbasno-1+6*(ina-1)+4), norme)
            call normev(zr(jbasno-1+6*(inb-1)+1), norme)
            call normev(zr(jbasno-1+6*(inb-1)+4), norme)
!
!         CORRECTION DU VECTEUR DE PROPAGATION AUX EXTREMITES DU FOND
!         SE TROUVANT SUR UN BORD DE LA STRUCTURE
            if (((iseg .eq. 1) .or. (iseg .eq. nseg)) .and. tyfond .ne. 'FERME') then
                call fonext(noma, cnxinv, jbasno, inoext, inoseg, &
                            nbnoff, jborl, jdirol, jnvdir, iseg)
            end if
!
        else if (ndim .eq. 2) then
!
            do j = 1, 4
                zr(jbasno-1+4*(ina-1)+j) = zr(jbasse-1+4*(iseg-1)+j)
            end do
!
        end if
!
!     DETERMINATION DE LA TAILLE DE MAILLE POUR LE NOEUD NA
!     (ET LE NOEUD MILIEU INC EN QUADRATIQUE)
        do k = 1, ndim
            vecdir(k) = zr(jbasno-1+2*ndim*(ina-1)+k+ndim)
        end do
!
        call fonno7(noma, cnxinv, ndim, na, vecdir, &
                    hmax)
!
        if ((casfon .eq. 'QUADRATIQUE') .and. (iseg .eq. 1)) then
            zr(jtail) = hmax
        else if ((casfon .eq. 'QUADRATIQUE') .and. (iseg .ne. 1)) then
!          TAILLE POUR LE NOEUD SOMMET NA COURANT
            zr(jtail-1+2*(iseg-1)+1) = hmax
!          RECUP TAILLE NOEUD SOMMET NA PRECEDENT
            hmaxpr = zr(jtail-1+2*(iseg-2)+1)
!          TAILLE POUR LE NOEUD MILIEU DU SEGMENT PRECEDENT
            zr(jtail-1+2*(iseg-2)+2) = (hmax+hmaxpr)/2
        else
            zr(jtail-1+iseg) = hmax
        end if
!
    end do
!
!     DANS LE CAS D'UN FOND FERME: CORRECTION DE BASEFOND
!     AUX NOEUDS EXTREMITES
    if ((tyfond .eq. 'FERME') .and. (ndim .eq. 3)) then
        do j = 1, 6
            zr(jbasno-1+6*(1-1)+j) = (zr(jbasse-1+6*(1-1)+j)+zr( &
                                      jbasse-1+6*(nseg-1)+j))/2.d0
            zr(jbasno-1+6*(nbnoff-1)+j) = zr(jbasno-1+6*(1-1)+j)
        end do
        call normev(zr(jbasno-1+6*(1-1)+1), norme)
        call normev(zr(jbasno-1+6*(1-1)+4), norme)
        call normev(zr(jbasno-1+6*(nbnoff-1)+1), norme)
        call normev(zr(jbasno-1+6*(nbnoff-1)+4), norme)
    end if
!
!     DETERMINATION DE LA TAILLE DE MAILLE POUR LE DERNIER NOEUD
    if (ndim .eq. 3) then
!
        do k = 1, ndim
            vecdir(k) = zr(jbasno-1+2*ndim*(inb-1)+k+ndim)
        end do
!
        call fonno7(noma, cnxinv, ndim, nb, vecdir, &
                    hmax)
        zr(jtail-1+nbnoff) = hmax
!
!         TAILLE DE MAILLE POUR LE DERNIER NOEUD MILIEU
        if (casfon .eq. 'QUADRATIQUE') then
            hmaxpr = zr(jtail-1+nbnoff-2)
            zr(jtail-1+nbnoff-1) = (hmax+hmaxpr)/2
        end if
!
!        DANS LE CAS D'UN FOND FERME: CORRECTION DE LA TAILLE DE
!        MAILLE AU PREMIER NOEUD
        if (tyfond .eq. 'FERME') zr(jtail-1+1) = zr(jtail-1+nbnoff)
!
    end if
!
!     MENAGE
    call jedetr(basseg)
    call jedetr('&&FONNOR2.LBORD')
    call jedetr('&&FONNOR2.VDIROL')
    call jedetr('&&FONNOR2.NVDIR')
!
    call jedema()
end subroutine
