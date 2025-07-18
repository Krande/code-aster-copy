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

subroutine xptfon(noma, ndim, nmafon, cnslt, cnsln, &
                  cnxinv, jmafon, nxptff, jfon, nfon, &
                  jbas, jtail, fiss, goinop, listpt, &
                  orient, typdis, nbmai, operation_opt)
!
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/cengra.h"
#include "asterfort/cnocns.h"
#include "asterfort/confac.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/intfac.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/normev.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "asterfort/xextre.h"
#include "asterfort/xfabor.h"
#include "asterfort/xfocoh.h"
#include "asterfort/xnorme.h"
#include "asterfort/xtailm.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/as_allocate.h"
#include "asterfort/int_to_char8.h"
!
    integer(kind=8) :: nmafon, jmafon, jfon, nfon, jbas, jtail, nxptff
    character(len=8) :: noma, fiss
    character(len=16) :: typdis
    character(len=16), intent(in), optional :: operation_opt
    character(len=19) :: cnslt, cnsln, cnxinv, listpt
    aster_logical :: orient, goinop
!     ------------------------------------------------------------------
! person_in_charge: samuel.geniaut at edf.fr
!
!       RECHERCHE DES POINTS DU FOND DE FISSURE DANS LE CADRE DE XFEM
!
!  ENTREES :
!     NOMA         :    NOM DE L'OBJET MAILLAGE
!     NMAFON       :    NOMBRE DE MAILLES DE LA ZONE FOND DE FISSURE
!     JMAFON       :    MAILLES DE LA ZONE FOND DE FISSURE
!     NXPTFF       :    NOMBRE MAXIMUM DE POINTS DU FOND DE FISSURE
!     CNSLT,CNSLN  :    LEVEL-SETS
!     CNXINV       :    CONNECTIVITE INVERSE
!     FISS         :    SD FISS_XFEM (POUR RECUP DES GRADIENTS)
!     GOINOP       :    .TRUE. SI OPOO10 AVEC UPWIND-SIMPLEXE/GRILLE/3D
!                       .FALSE. SINON
!
!  SORTIES :
!     JFON         :   ADRESSE DES POINTS DU FOND DE FISSURE
!     JBAS         :   ADRESSE DES DIRECTIONS DE PROPAGATION
!     NFON         :   NOMBRE DE POINTS DU FOND DE FISSURE
!     LISTPT       :   LISTE DES INDICES DES POINTS DU FOND DESORDONNES
!     ORIENT       :   .TRUE.  SI LE FOND PEUT ETRE ORIENTE
!                      .FALSE. SI LE FOND NE PEUT PAS ETRE ORIENTE
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: ipt, ima, i, j, ibid, ndim, ino, k, ifq, ifm, niv
    integer(kind=8) :: nmaabs, nbf, nbnoma, nuno, nunoa, nunob, nunoc, nunod
    integer(kind=8) :: fa(6, 8), ibid3(12, 3), vecind(5)
    integer(kind=8) :: jconx1, jconx2, jcoor
    integer(kind=8) :: jglsn, jglst, igeom, itypma
    integer(kind=8) :: indipt, jborl, jdirol, jnvdir, jlistp
    integer(kind=8) :: nbfacb, iptbor(2), nbptma, ndime, indptf(3), codret
    integer(kind=8) :: nunopa, nunopb, nunopc, nunopd
    integer(kind=8) :: snuno, pnuno, inuno, snunop, pnunop, inunop
    integer(kind=8) :: nbmai
    real(kind=8) :: m(3), p(3), gln(3), glt(3), coorg(3), vectn(12)
    real(kind=8) :: normi, r3bid(3)
    character(len=8) :: typma, nommai, alias
    character(len=19) :: grlt, chgrt, grln, chgrn
    character(len=16) :: operation
    aster_logical :: fabord, indic
    real(kind=8), pointer :: lsn(:) => null()
    real(kind=8), pointer :: lst(:) => null()
    aster_logical, pointer :: ptbord(:) => null()
    integer(kind=8), pointer :: typmail(:) => null()
    real(kind=8), pointer :: gn(:) => null()
    real(kind=8), pointer :: gt(:) => null()
    real(kind=8), pointer :: lnsv(:) => null()
    real(kind=8), pointer :: ltsv(:) => null()
! ----------------------------------------------------------------------
    call jemarq()
    call infdbg('XFEM', ifm, niv)
!
    if (present(operation_opt)) then
        operation = operation_opt
    else
        operation = 'RIEN'
    end if
!
    call jeveuo(noma//'.COORDO    .VALE', 'L', jcoor)
    call jeveuo(noma//'.CONNEX', 'L', jconx1)
    call jeveuo(jexatr(noma//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call jeveuo(noma//'.TYPMAIL', 'L', vi=typmail)
!
    call jeveuo(cnslt//'.CNSV', 'L', vr=ltsv)
    call jeveuo(cnsln//'.CNSV', 'L', vr=lnsv)
!
    if (operation .ne. 'FRONT_COHE_BRUT') then
!
!       GRADIENT LST
        if (goinop) then
            grlt = fiss//'.GRI.GRLTNO'
        else
            grlt = fiss//'.GRLTNO'
        end if
        chgrt = '&&XPTFON.GRLN'
        call cnocns(grlt, 'V', chgrt)
        call jeveuo(chgrt//'.CNSV', 'L', vr=gt)
!
!       GRADIENT LSN
        if (goinop) then
            grln = fiss//'.GRI.GRLNNO'
        else
            grln = fiss//'.GRLNNO'
        end if
        chgrn = '&&XPTFON.GRLT'
        call cnocns(grln, 'V', chgrn)
        call jeveuo(chgrn//'.CNSV', 'L', vr=gn)
!
!       Initialisation de fissure cohesive
        if (typdis .eq. 'COHESIF' .and. operation .eq. 'RIEN') then
            call xfocoh(jbas, jconx1, jconx2, jcoor, jfon, &
                        cnsln, chgrn, chgrt, noma, listpt, ndim, nfon, &
                        nxptff, orient, nbmai)
        end if
!
!       on refait le test pour eviter un Warning à la compilation
        if (typdis .eq. 'COHESIF' .and. operation .eq. 'RIEN') goto 999
!
        AS_ALLOCATE(vl=ptbord, size=nxptff)
!
!       VECTEUR PERMETTANT DE SAVOIR SI LE VECTEUR DE DIRECTION DE
!       PROPAGATION (VDIR) A ETE RECALCULE OU NON AUX POINTS
!       EXTREMITES DE FONFIS
        call wkvect('&&XPTFON.LBORD', 'V V L', nxptff, jborl)
!
!         VECTEUR CONTENANT LES VDIR INITIAUX (CAD SANS MODIFICATION
!         DES VECTEURS AUX POINTS EXTREMITES DE FONFIS)
        call wkvect('&&XPTFON.VDIROL', 'V V R', 3*nxptff, jdirol)
!
!         VECTEUR CONTENANT 0 OU 1 AUX POINTS EXTREMITES DE FONFIS:
!         0: LE PRODUIT SCALAIRE ENTRE LA NORMALE A LA FACE DE BORD ET
!            LE VDIR INITIAL ESI INFERIEUR A 0
!         1: LE PRODUIT SCALAIRE EST SUPERIEUR OU EGAL A 0
        call wkvect('&&XPTFON.NVDIR', 'V V I', nxptff, jnvdir)
!
        do i = 1, nxptff
            ptbord(i) = .false.
            zl(jborl-1+i) = .false.
        end do
    end if
!
!     COMPTEUR : NOMBRE DE POINTS DE FONFIS TROUVES
    ipt = 0
!
    orient = .true.
!     LISTE DES INDICES DES POINTS DU FOND DESORDONNES
    call wkvect(listpt, 'V V I', 2*nmafon, jlistp)
!
!     BOUCLE SUR LES MAILLES DE MAFOND
    do ima = 1, nmafon
!
        nmaabs = zi(jmafon-1+ima)
        itypma = typmail(nmaabs)
        call jenuno(jexnum('&CATA.TM.NOMTM', itypma), typma)
        call dismoi('DIM_TOPO', typma, 'TYPE_MAILLE', repi=ndime)
!
!       ON SE RECREE UN ENVIRONNEMENT COMME DANS UN TE
!       POUR LSN, LST, GRLST, GRLST ET IGEOM
!       AFIN DE POUVOIR APPELER INTFAC
        nbnoma = zi(jconx2+nmaabs)-zi(jconx2+nmaabs-1)
        AS_ALLOCATE(vr=lsn, size=nbnoma)
        AS_ALLOCATE(vr=lst, size=nbnoma)
        if (operation .ne. 'FRONT_COHE_BRUT') then
            call wkvect('&&XPTFON.GRLSN', 'V V R', nbnoma*ndim, jglsn)
            call wkvect('&&XPTFON.GRLST', 'V V R', nbnoma*ndim, jglst)
        end if
        call wkvect('&&XPTFON.IGEOM', 'V V R', nbnoma*ndim, igeom)
        do ino = 1, nbnoma
            nuno = zi(jconx1-1+zi(jconx2+nmaabs-1)+ino-1)
            lsn(ino) = lnsv(nuno)
            lst(ino) = ltsv(nuno)
            do j = 1, ndim
                if (operation .ne. 'FRONT_COHE_BRUT') then
                    zr(jglsn-1+ndim*(ino-1)+j) = gn(ndim*(nuno-1)+j)
                    zr(jglst-1+ndim*(ino-1)+j) = gt(ndim*(nuno-1)+j)
                end if
                zr(igeom-1+ndim*(ino-1)+j) = zr(jcoor-1+3*(nuno-1)+j)
            end do
        end do
!
!       CALCUL DU CENTRE DE GRAVITE DE LA MAILLE
        call cengra(noma, nmaabs, coorg)
!
        call confac(typma, ibid3, ibid, fa, nbf)
!
!       INITIALISATIONS
        nbfacb = 0
        iptbor(1) = 0
        iptbor(2) = 0
        nbptma = 0
        do i = 1, 3
            indptf(i) = 0
        end do
!
!       BOUCLE SUR LES FACES
        do ifq = 1, nbf
!         TYPE DE FACE
            if (fa(ifq, 4) .eq. 0) then
                alias = 'TR3'
            else
                alias = 'QU4'
            end if
!
            indipt = 0
!         RECHERCHE DES INTERSECTION ENTRE LE FOND DE FISSURE ET LA FACE
            if (operation .ne. 'FRONT_COHE_BRUT') then
                call intfac(noma, nmaabs, ifq, fa, nbnoma, &
                            lst, lsn, ndim, 'OUI', jglsn, &
                            jglst, igeom, m, indptf, gln, &
                            glt, codret)
            else
                call intfac(noma, nmaabs, ifq, fa, nbnoma, &
                            lst, lsn, ndim, 'NON', ibid, &
                            ibid, igeom, m, indptf, r3bid, &
                            r3bid, codret)
            end if
!
            if (codret .eq. 0) goto 200
!
!         VERIFICATION SI CE POINT A DEJA ETE TROUVE
            nunoa = zi(jconx1-1+zi(jconx2+nmaabs-1)+fa(ifq, 1)-1)
            nunob = zi(jconx1-1+zi(jconx2+nmaabs-1)+fa(ifq, 2)-1)
            nunoc = zi(jconx1-1+zi(jconx2+nmaabs-1)+fa(ifq, 3)-1)
            if (alias .eq. 'QU4') then
                nunod = zi(jconx1-1+zi(jconx2+nmaabs-1)+fa(ifq, 4)-1)
            end if
!
            do j = 1, ipt
                indic = .false.
!
!           VERIFICATION POUR UN POINT INTERIEUR (cf Doc R7.02.12)
                if ((indptf(1) .eq. 3) .and. (int(zr(jfon-1+11*(j-1)+9)) .eq. 3)) then
                    nunopa = int(zr(jfon-1+11*(j-1)+5))
                    nunopb = int(zr(jfon-1+11*(j-1)+6))
                    nunopc = int(zr(jfon-1+11*(j-1)+7))
                    snuno = nunoa+nunob+nunoc
                    pnuno = nunoa*nunob*nunoc
                    inuno = nunoa*nunob+nunoa*nunoc+nunob*nunoc
                    snunop = nunopa+nunopb+nunopc
                    pnunop = nunopa*nunopb*nunopc
                    inunop = nunopa*nunopb+nunopa*nunopc+nunopb*nunopc
                    if (alias .eq. 'QU4') then
                        nunopd = int(zr(jfon-1+11*(j-1)+8))
                        snuno = nunoa+nunob+nunoc+nunod
                        pnuno = nunoa*nunob*nunoc*nunod
                        inuno = nunoa*nunob+nunoa*nunoc+nunob*nunoc+ &
                                nunoa*nunod+nunob*nunod+nunoc*nunod
                        snunop = nunopa+nunopb+nunopc+nunopd
                        pnunop = nunopa*nunopb*nunopc*nunopd
                        inunop = nunopa*nunopb+nunopa*nunopc+nunopb* &
                                 nunopc+nunopa*nunopd+nunopb*nunopd+nunopc* &
                                 nunopd
                    end if
                    if ((snuno .eq. snunop) .and. (pnuno .eq. pnunop) .and. (inuno .eq. &
                                                                             inunop)) then
                        indic = .true.
                    end if
!           VERIFICATION POUR UN POINT ARETE (cf Doc R7.02.12)
                else if ((indptf(1) .eq. 2) .and. (int(zr(jfon-1+11*(j-1) &
                                                          +9)) .eq. 2)) then
                    snuno = indptf(2)+indptf(3)
                    pnuno = (indptf(2))*(indptf(3))
                    nunopa = int(zr(jfon-1+11*(j-1)+10))
                    nunopb = int(zr(jfon-1+11*(j-1)+11))
                    snunop = nunopa+nunopb
                    pnunop = nunopa*nunopb
                    if ((snuno .eq. snunop) .and. (pnuno .eq. pnunop)) then
                        indic = .true.
                    end if
!           VERIFICATION POUR UN POINT SOMMET (cf Doc R7.02.12)
                else if ((indptf(1) .eq. 1) .and. (int(zr(jfon-1+11*(j-1) &
                                                          +9)) .eq. 1)) then
                    nunopa = int(zr(jfon-1+11*(j-1)+10))
                    if (indptf(2) .eq. nunopa) then
                        indic = .true.
                    end if
                end if
!
                if (indic) then
                    if (operation .ne. 'FRONT_COHE_BRUT') then
!               CALCUL DE LA TAILLE MAX DE LA MAILLE IMA ET MISE A
!               JOUR DU VECTEUR DES TAILLES DE MAILLES
                        call xtailm(ndim, glt, nmaabs, typma, jcoor, &
                                    jconx1, jconx2, j, jtail)
                    end if
                    if (ndim .ne. 3) goto 200
                    indipt = j
                    goto 300
                end if
            end do
!
!         CE POINT N'A PAS DEJA ETE TROUVE, ON LE GARDE
            ipt = ipt+1
!         AUGMENTER NXPTFF
            ASSERT(ipt .le. nxptff)
!
!         STOCKAGE DES COORDONNEES DU POINT M,
!         DE LA BASE LOCALE (GRADIENT DE LSN ET LST) ET
!         DES NUMEROS DES SOMMETS DE LA FACE CONTENANT M
            do k = 1, ndim
                zr(jfon-1+11*(ipt-1)+k) = m(k)
                if (operation .ne. 'FRONT_COHE_BRUT') then
                    zr(jbas-1+2*ndim*(ipt-1)+k) = gln(k)
                    zr(jbas-1+2*ndim*(ipt-1)+k+ndim) = glt(k)
                end if
            end do
!           l'ecriture des numeros des sommets de la face contenant M dans zr(jfon) est
!           indispensable pour le remplissage du vecteur .NOFACPTFON dans xffcr et xlmail
            do k = 1, 3
                zr(jfon-1+11*(ipt-1)+4+k) = zi(jconx1-1+zi(jconx2+nmaabs-1)+fa(ifq, k)-1)
                zr(jfon-1+11*(ipt-1)+8+k) = indptf(k)
            end do
            if (alias .eq. 'QU4') then
                zr(jfon-1+11*(ipt-1)+8) = zi(jconx1-1+zi(jconx2+nmaabs-1)+fa(ifq, 4)-1)
            else
                zr(jfon-1+11*(ipt-1)+8) = 0
            end if
!
            indipt = ipt
!
!         CALCUL DE LA TAILLE MAX DE LA MAILLE IMA
            if (operation .ne. 'FRONT_COHE_BRUT') then
                call xtailm(ndim, glt, nmaabs, typma, jcoor, &
                            jconx1, jconx2, ipt, jtail)
            end if
!
            if (ndim .ne. 3) goto 200
!
300         continue
!
!           ON VERIFIE SI LA FACE COURANTE EST UNE FACE DE BORD
!           CELA N'A DE SENS QU'EN 3D
            if (operation .ne. 'FRONT_COHE_BRUT') then
                call xfabor(noma, cnxinv, nunoa, nunob, nunoc, &
                            fabord)
!
!               SI LA FACE EST UNE FACE DE BORD ON PREND SA NORMALE
                if (fabord) then
                    if (.not. ptbord(indipt)) then
                        ptbord(indipt) = .true.
                    end if
                    if (ndim .eq. 3) then
                        call xnorme(indipt, iptbor, vectn, nbfacb, nunoa, &
                                    nunob, nunoc, jcoor, coorg)
                    end if
                end if
            end if
!
!         STOCKAGE DES INDICES DES POINTS DU FOND APPARTENANT A
!         LA MAILLE
            if (ndime .eq. 3) then
                if (nbptma .eq. 0) then
                    zi(jlistp-1+2*(ima-1)+1) = indipt
                    nbptma = nbptma+1
                    vecind(nbptma) = indipt
!
                elseif ((nbptma .eq. 1) .and. (indipt .ne. zi(jlistp-1+2*( &
                                                              ima-1)+1))) then
                    zi(jlistp-1+2*(ima-1)+2) = indipt
                    nbptma = nbptma+1
                    vecind(nbptma) = indipt
!
                else if (nbptma .gt. 1) then
                    do j = 1, nbptma
                        if (indipt .eq. vecind(j)) goto 400
                    end do
!             ORIENTATION DU FOND PAS POSSIBLE A REALISER
                    orient = .false.
!
                    nbptma = nbptma+1
                    vecind(nbptma) = indipt
!
!             IMPRESSION DES POINTS DU FOND APPARTENANT A LA MAILLE
                    if (ifq .eq. nbf) then
                        nommai = int_to_char8(nmaabs)
                        call utmess('I', 'XFEM_51', sk=nommai, si=nbptma)
!
                        write (ifm, 797)
797                     format(7x, 'X', 13x, 'Y', 13x, 'Z')
                        do j = 1, nbptma
                            p(1) = zr(jfon-1+11*(vecind(j)-1)+1)
                            p(2) = zr(jfon-1+11*(vecind(j)-1)+2)
                            p(3) = zr(jfon-1+11*(vecind(j)-1)+3)
                            write (ifm, 798) (p(k), k=1, 3)
798                         format(2x, 3(e12.5, 2x))
                        end do
                    end if
                end if
            end if
400         continue
!       FIN BOUCLE SUR LES FACES DE LA MAILLE
200         continue
        end do
!
!       CALCUL DES VECTEURS DE PROPAGATION AUX EXTREMITES
        if ((ndim .eq. 3) .and. (nbfacb .ne. 0) .and. &
            operation .ne. 'FRONT_COHE_BRUT') then
            call xextre(iptbor, vectn, nbfacb, jbas, jborl, &
                        jdirol, jnvdir)
        end if
!
!       DESTRUCTION DES VECTEURS LOCAUX A LA MAILLE
        AS_DEALLOCATE(vr=lsn)
        AS_DEALLOCATE(vr=lst)
        call jedetr('&&XPTFON.GRLSN')
        call jedetr('&&XPTFON.GRLST')
        call jedetr('&&XPTFON.IGEOM')
!
!     FIN BOUCLE SUR LES MAILLES
    end do
!
!     NORMALISATION DES NOUVEAUX VECTEURS DE DIRECTION DE
!     PROPAGATION
    if (operation .ne. 'FRONT_COHE_BRUT' .and. ndim .eq. 3) then
        do k = 1, nxptff
            if (ptbord(k)) then
                call normev(zr(jbas-1+6*(k-1)+4), normi)
            end if
        end do
!
    end if
!
    nfon = ipt
!
!     IMPRESSION DE TOUS LES POINTS DE FOND DE FISSURE
!     DANS LE CAS OU IL N'EST PAS POSSIBLE DE LES ORDONNER
    if (.not. orient) then
        call utmess('A', 'XFEM_52')
        call utmess('I', 'XFEM_35')
        write (ifm, 897)
897     format(7x, 'X', 13x, 'Y', 13x, 'Z')
        do i = 1, nfon
            p(1) = zr(jfon-1+11*(i-1)+1)
            p(2) = zr(jfon-1+11*(i-1)+2)
            p(3) = zr(jfon-1+11*(i-1)+3)
            write (ifm, 898) (p(k), k=1, 3)
898         format(2x, 3(e12.5, 2x))
        end do
    end if
!
!   IMPRESSION DU FRONT BRUT EN MODE DETECTION
    if (operation .eq. 'FRONT_COHE_BRUT' .and. niv .gt. 1) then
    do i = 1, nfon
        p(1) = zr(jfon-1+11*(i-1)+1)
        p(2) = zr(jfon-1+11*(i-1)+2)
        p(3) = zr(jfon-1+11*(i-1)+3)
        write (6, 899) (p(k), k=1, 3)
899     format(2x, 3(e12.5, 2x))
    end do
    end if
!
!
    if (operation .ne. 'FRONT_COHE_BRUT') then
        call jedetr('&&XPTFON.LBORD')
        call jedetr('&&XPTFON.VDIROL')
        call jedetr('&&XPTFON.NVDIR')
    end if
    AS_DEALLOCATE(vl=ptbord)
999 continue
    if (operation .ne. 'FRONT_COHE_BRUT') then
        call jedetr(chgrn)
        call jedetr(chgrt)
    end if
    call jedema()
end subroutine
