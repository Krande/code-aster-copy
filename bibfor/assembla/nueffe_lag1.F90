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
subroutine nueffe_lag1(nbLigr, listLigr, base, numeDofZ, renumZ, &
                       modelZ, modeLocZ, idenRelaZ)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/creprn.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jecrec.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nddl.h"
#include "asterfort/nudeeq_lag1.h"
#include "asterfort/nulili.h"
#include "asterfort/nuno1.h"
#include "asterfort/nunueq.h"
#include "asterfort/renuno.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nbLigr
    character(len=24), pointer :: listLigr(:)
    character(len=2), intent(in) :: base
    character(len=*), intent(in) :: numeDofZ, renumZ, modelZ
    character(len=*), intent(in) :: modeLocZ, idenRelaZ
!
! --------------------------------------------------------------------------------------------------
!
! Factor
!
! Numbering - Create NUME_EQUA objects (with simple Lagrange)
!
! --------------------------------------------------------------------------------------------------
!
! In  nbLigr         : number of LIGREL in list
! Ptr listLigr       : pointer to list of LIGREL
! In  numeDof        : name of numbering object (NUME_DDL)
! In  base           : JEVEUX base to create objects
!                      base(1:1) => NUME_DDL objects
!                      base(2:2) => NUME_EQUA objects
! In  renum          : method for renumbering equations (SANS/RCMK)
! In  modeLoc        : local mode for GRANDEUR numbering
! In  idenRela       : name of object for identity relations between dof
!
! --------------------------------------------------------------------------------------------------
!
! Attention : ne fait pas jemarq/jedema car nulili
!             recopie des adresses jeveux dans .ADNE et .ADLI
!             Ces objets seront utilises pendant la creation de la sd "stockage" (promor.F90)
!
! --------------------------------------------------------------------------------------------------
!
! Cette routine cree les objets suivants :
!  nume(1:14)//     .ADLI
!                   .ADNE
!              .NUME.DEEQ
!              .NUME.DELG
!              .NUME.LILI
!              .NUME.NEQU
!              .NUME.NUEQ
!              .NUME.PRNO
!              .NUME.REFN
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: modeLoc, idenRela
    character(len=8) :: gran_name, kbid
    character(len=3) :: kret
    integer(kind=8) :: n, igds, nec, nlili
    character(len=8) :: nomcmp
    character(len=8) :: mesh
    character(len=14) :: numeDof
    character(len=16) :: nomte
    character(len=24) :: nnli, num21, nuno, nomli
    character(len=24) :: num2, dsclag, exi1, newn, oldn, derli
    character(len=19) :: nume_equa
    character(len=24) :: nequ, refn
    character(len=24) :: lili, prno, nueq, deeq, delg
    integer(kind=8) :: nb_node_mesh, ilim, itypel, nb_dof, jdeeq, jdelg, nb_equa
    integer(kind=8) :: nb_iden_rela, nb_iden_dof, nb_iden_term
    integer(kind=8) :: i, iad, ianueq, icddlb
    integer(kind=8) :: iconx1, iconx2, iddlag, iderli
    integer(kind=8) :: idnocm, idprn1, idprn2, idref
    integer(kind=8) :: iec, iel, iexi1, ifm, igr, ilag, n0, n1, n2, nn, n22
    integer(kind=8) :: ili, inewn, ino, inum21
    integer(kind=8) :: inuno2, ioldn, iprnm, ire, iret
    integer(kind=8) :: j, jprno, k, l
    integer(kind=8) :: nbcmp, nbn, nb_node_subs
    integer(kind=8) :: nb_node, nbnonu, nbnore, nddl1, nddlb
    integer(kind=8) :: nel, niv, nlag, nno, nblag
    integer(kind=8) :: numa, nunoel, nume_late_node, nunoel_save
    integer(kind=8) :: vali(5)
    aster_logical :: lparallel_mesh
    integer(kind=8), pointer :: v_nnli(:) => null()
    integer(kind=8), pointer :: adli(:) => null()
    integer(kind=8), pointer :: bid(:) => null()
    integer(kind=8), pointer :: adne(:) => null()
    integer(kind=8), pointer :: qrns(:) => null()
    integer(kind=8), pointer :: p_nequ(:) => null()
    integer(kind=8), pointer :: idenRelaInfo(:) => null()
    integer(kind=8), pointer :: lagr_mult(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
!     NBNOM  : NOMBRE DE NOEUDS DU MAILLAGE
!     DERLI  : NOM DE L'OBJET NU.DERLI CREE SUR 'V'
!              DERLI(N3)= MAX DES N1 TELS QUE IL EXISTE UNE MAILLE SUP
!              DE TYPE SEG3 MS TELLE QUE N1 1ER, N3 3EME NOEUD DE MS
!     DSCLAG : NOM DE L'OBJET NU.DSCLAG CREE SUR 'V'
!              DIM=3*NBRE LE LAGR.
!              SI ILAG LAGRANGE DE BLOCAGE
!              DSCLAG(3*(ILAG-1)+1)= +NUM DU NOEUD PH. BLOQUE
!              DSCLAG(3*(ILAG-1)+2)= -NUM DU DDL DU NOEUD PH. BLOQUE
!              DSCLAG(3*(ILAG-1)+3)= +1 (LAGR).
!              SI ILAG LAGRANGE DE LIAISON
!              DSCLAG(3*(ILAG-1)+1)= 0
!              DSCLAG(3*(ILAG-1)+2)= 0
!              DSCLAG(3*(ILAG-1)+3)= +1 (LAGR).
!-----------------------------------------------------------------------
!     FONCTIONS LOCALES D'ACCES AUX DIFFERENTS CHAMPS DES
!     S.D. MANIPULEES DANS LE SOUS PROGRAMME
!-----------------------------------------------------------------------

!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS LIEL DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL
!     ZZLIEL(ILI,IGREL,J) =
!      SI LA JIEME MAILLE DU LIEL IGREL DU LIGREL ILI EST:
!          -UNE MAILLE DU MAILLAGE : SON NUMERO DANS LE MAILLAGE
!          -UNE MAILLE TARDIVE : -POINTEUR DANS LE CHAMP .NEMA

#define zzliel(ili,igrel,j) zi(adli(1+3*(ili-1)+1)-1+ zi(adli(1+3*(ili-1)+2)+igrel-1)+j-1)

!---- NBRE DE GROUPES D'ELEMENTS (DE LIEL) DU LIGREL ILI

#define zzngel(ili) adli(1+3* (ili-1))

!---- NBRE DE NOEUDS DE LA MAILLE TARDIVE IEL ( .NEMA(IEL))
!     DU LIGREL ILI REPERTOIRE .LILI
!     (DIM DU VECTEUR D'ENTIERS .LILI(ILI).NEMA(IEL) )

#define zznsup(ili,iel) zi(adne(1+3* (ili-1)+2)+iel) - zi(adne(1+3*(ili-1)+2)+iel-1 ) - 1

!---- NBRE D ELEMENTS DU LIEL IGREL DU LIGREL ILI DU REPERTOIRE .LILI
!     (DIM DU VECTEUR D'ENTIERS .LILI(ILI).LIEL(IGREL) )

#define zznelg(ili,igrel) zi(adli(1+3*(ili-1)+2)+igrel) - zi(adli(1+3*(ili-1)+2)+igrel-1) - 1

!---- NBRE D ELEMENTS SUPPLEMENTAIRE (.NEMA) DU LIGREL ILI DE .LILI

#define zznels(ili) adne(1+3* (ili-1))

!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS NEMA DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL
!     ZZNEMA(ILI,IEL,J) =  1.LE. J .GE. ZZNELS(ILI)
!      SI LE J IEME NOEUD DE LA MAILE TARDIVE IEL DU LIGREL ILI EST:
!          -UN NOEUD DU MAILLAGE : SON NUMERO DANS LE MAILLAGE
!          -UN NOEUD TARDIF : -SON NUMERO DANS LA NUMEROTATION LOCALE
!                              AU LIGREL ILI
!     ZZNEMA(ILI,IEL,ZZNELS(ILI)+1)=NUMERO DU TYPE_MAILLE DE LA MAILLE
!                                   IEL DU LIGREL ILI

#define zznema(ili,iel,j) zi( adne(1+3* (ili-1)+1)-1+ zi(adne(1+3* (ili-1)+2)+iel-1 )+j-1 )

!---- FONCTION D ACCES AUX ELEMENTS DES CHAMPS PRNO DES S.D. LIGREL
!     REPERTORIEES DANS LE CHAMP LILI DE NUME_DDL ET A LEURS ADRESSES
!     ZZPRNO(ILI,NUNOEL,1) = NUMERO DE L'EQUATION ASSOCIEES AU 1ER DDL
!                            DU NOEUD NUNOEL DANS LA NUMEROTATION LOCALE
!                            AU LIGREL ILI DE .LILI
!     ZZPRNO(ILI,NUNOEL,2) = NOMBRE DE DDL PORTES PAR LE NOEUD NUNOEL
!     ZZPRNO(ILI,NUNOEL,2+1) = 1ER CODE
!     ZZPRNO(ILI,NUNOEL,2+NEC) = NEC IEME CODE

#define izzprn(ili,nunoel,l) (idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
#define zzprno(ili,nunoel,l) zi( idprn1-1+zi(idprn2+ili-1)+ (nunoel-1)* (nec+2)+l-1)
!
! --------------------------------------------------------------------------------------------------
!
!    call jemarq() FORBIDDEN !
    call infniv(ifm, niv)
!
    numeDof = numeDofZ
    modeLoc = modeLocZ
    idenRela = idenRelaZ

! - Identity relations between dof
    nb_iden_rela = 0
    nb_iden_dof = 0
    nb_iden_term = 0
    if (idenRela .ne. ' ') then
        call jeveuo(idenRela(1:19)//'.INFO', 'L', vi=idenRelaInfo)
        nb_iden_rela = idenRelaInfo(1)
        nb_iden_term = idenRelaInfo(2)
        nb_iden_dof = idenRelaInfo(3)
    end if

! - SI LE CONCEPT : NU EXISTE DEJA, ON LE DETRUIT COMPLETEMENT :
    call detrsd('NUME_DDL', numeDof)

! - NOMS DES PRINCIPAUX OBJETS JEVEUX
    nume_equa = numeDof//'.NUME'
    lili = nume_equa(1:19)//'.LILI'
    prno = nume_equa(1:19)//'.PRNO'
    nueq = nume_equa(1:19)//'.NUEQ'
    deeq = nume_equa(1:19)//'.DEEQ'
    delg = nume_equa(1:19)//'.DELG'
    nequ = nume_equa(1:19)//'.NEQU'
    refn = nume_equa(1:19)//'.REFN'

! --- Objets temporaires
    nnli = numeDof//'.NNLI'
    nuno = numeDof//'.NUNO'
    exi1 = numeDof//'.EXI1'
    newn = numeDof//'.NEWN'
    oldn = numeDof//'.OLDN'
    derli = numeDof//'.DERLI'
    num21 = numeDof//'.NUM21'
    num2 = numeDof//'.NUM2'
    dsclag = numeDof//'.DESCLAG'

! - Create LILI objects
    call nulili(nbLigr, listLigr, lili, base(2:2), gran_name, &
                igds, mesh, nec, nlili, modeLocZ_=modeLoc)
    call jeveuo(numeDof//'     .ADLI', 'E', vi=adli)
    call jeveuo(numeDof//'     .ADNE', 'E', vi=adne)

! - Access to mesh objects
    call jeexin(mesh(1:8)//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(mesh(1:8)//'.CONNEX', 'L', iconx1)
        call jeveuo(jexatr(mesh(1:8)//'.CONNEX', 'LONCUM'), 'L', iconx2)
    end if
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=kret)
    lparallel_mesh = (kret .eq. 'OUI')
    call dismoi('NB_NO_MAILLA', mesh, 'MAILLAGE', repi=nb_node_mesh)
    call dismoi('NB_NL_MAILLA', mesh, 'MAILLAGE', repi=nb_node_subs)
    nb_node = nb_node_mesh+nb_node_subs

! --- LILI(1)='&MAILLA'
    ilim = 1

! --- ALLOCATION DE L'OBJET NU.NNLI NOMBRE DE NOEUDS DECLARES DANS
! --- LE LIGREL ILI DE LILI :
!     ---------------------
    call wkvect(nnli, 'V V I', nlili, vi=v_nnli)
    v_nnli(1) = nb_node
    call jecrec(nuno, 'V V I ', 'NU', 'CONTIG', 'VARIABLE', nlili)

! --- ALLOCATION DE PRNO :
!     -------------------------------------------------
    call jecrec(prno, base(2:2)//' V I', 'NU', 'CONTIG', 'VARIABLE', nlili)

! --- CALCUL DE N, CALCUL DES NNLI ET DU POINTEUR DE LONGUEUR DE
! --- PRNO :
! --- NBNOM NOMBRE DE NOEUDS TOTAL DU MAILLAGE :
!     ------------------------------------------

    call jeecra(jexnum(nuno, 1), 'LONMAX', nb_node)
    call jeecra(jexnum(prno, 1), 'LONMAX', nb_node*(nec+2))

! --- N CONTIENDRA LE NOMBRE TOTAL (MAX) DE NOEUDS DE NUME_DDL
! --- TOUS LES NOEUDS DU MAILLAGE + TOUS LES NOEUDS SUPL. DES LIGRELS :
!     ---------------------------------------------------------------
    n = nb_node
    nblag = 0
    do ili = 2, nlili
        call jenuno(jexnum(lili, ili), nomli)
        call jeexin(nomli(1:19)//'.NBNO', iret)
        if (iret .ne. 0) then

! ---    ACCES AU NOMBRE DE NOEUDS SUP DU LIGREL DE NOM NOMLI :
!        ----------------------------------------------------
            call jeveuo(nomli(1:19)//'.NBNO', 'L', iad)
            nbn = zi(iad)
        else
            nbn = 0
        end if

! ---    AFFECTATION DU CHAMP .NNLI DE NU :
!        --------------------------------
        v_nnli(ili) = nbn
!
        call jeecra(jexnum(nuno, ili), 'LONMAX', nbn)
        call jeecra(jexnum(prno, ili), 'LONMAX', nbn*(nec+2))
        n = n+nbn
        nblag = nblag+nbn
    end do

    call wkvect(derli, 'V V I', n+1, iderli)

    call jeveuo(prno, 'E', idprn1)
    call jeveuo(jexatr(prno, 'LONCUM'), 'L', idprn2)

! --- ALLOCATION DE LA COLLECTION NUNO NUMEROTEE DE VECTEUR DE LONGUEUR
! --- NNLI ET DE NBRE D'ELMT NLILI SUR LA BASE VOLATILE (.NUNO ET .NNLI
! --- SONT SUPPRIMES DE LA S.D. NUME_EQUA) :
!     ------------------------------------
    call jeveuo(jexatr(nuno, 'LONCUM'), 'L', inuno2)
    nlag = zi(inuno2+nlili)-zi(inuno2+1)
    ASSERT(nlag == nblag)
    ASSERT(n == (nb_node+nlag))
!
! --- RENUMEROTATION , CREATION DES OBJETS NU.EXI1, NU.NEWN ET NU.OLDN :
!     ----------------------------------------------------------------

    call renuno(numeDof, renumZ)
    call jeveuo(exi1, 'L', iexi1)
    call jeveuo(newn, 'L', inewn)
    call jeveuo(oldn, 'L', ioldn)

    if (nlag > 0) then
        call wkvect(dsclag, ' V V I', 2*nlag, iddlag)
!
        do ili = 2, nlili
!
            call jenuno(jexnum(lili, ili), nomli)
            call jeexin(nomli(1:19)//'.MULT', iret)
            if (lparallel_mesh .and. iret .ne. 0) then
                call jeveuo(nomli(1:19)//'.MULT', 'L', vi=lagr_mult)
            end if
!
            do iel = 1, zznels(ili)
                nn = zznsup(ili, iel)
                if (nn .eq. 2) then
                    n1 = zznema(ili, iel, 1)
                    n2 = zznema(ili, iel, 2)
                    if ((n1 .gt. 0) .and. (n2 .lt. 0)) then

! ---    TRANSFORMATION DE N2 , NUMERO DU PREMIER LAGRANGE DANS LA
! ---    NUMEROTATION LOCALE AU LIGREL EN SON NUMERO DANS LA
! ---    NUMEROTATION GLOBALE :
!        --------------------
                        n22 = -n2
                        n2 = -n2
                        n2 = zi(inuno2+ili-1)+n2-1
                        ilag = n2-nb_node

! ---    RECUPERATION DU NOEUD PHYSIQUE DE NUMERO LE PLUS GRAND
! ---    LIE AU LAGRANGE PAR LE TABLEAU DERLI, CETTE
! ---    VALEUR N'EST DIFFERENTE DE 0 QUE S'IL S'AGIT D'UNE
! ---    RELATION LINEAIRE :
!        -----------------
                        n0 = zi(iderli+n2)
!
                        if (n0 .gt. 0) then
                            zi(iddlag+2*(ilag-1)) = 0
                            zi(iddlag+2*(ilag-1)+1) = 0
                        else
                            zi(iddlag+2*(ilag-1)) = n1
                            zi(iddlag+2*(ilag-1)+1) = -1
                            zi(iderli+n2) = n1
                        end if
!
                        if (lparallel_mesh .and. iret .ne. 0) then
                            if (lagr_mult(n22) .gt. 1) then
                                zi(iddlag+2*(ilag-1)) = 0
                                zi(iddlag+2*(ilag-1)+1) = 0
                            end if
                        end if
                    end if
                end if
            end do
        end do
    end if
!
! ---  NBNORE EST LE NOMBRE DE NOEUDS DU MAILLAGE PARTICIPANTS A LA
! ---  NUMEROTATION:
!      ------------
!
    call jelira(oldn, 'LONUTI', nbnore)
    ASSERT(nbnore <= nb_node)
!
!
! ---  ALLOCATION DE NUM21 ET NUM2
! ---  NUM2(I) = NUMERO DANS LA NOUVELLE NUMEROTATION NUNO DU NOEUD I
! ---  DANS LA NUMEROTATION GLOBALE DE 1ER NIVEAU
! ---  SI NUM2(I)=J ALORS NUM21(J)=I :
!      -----------------------------
!
    call wkvect(num21, ' V V I', n+1, inum21)
!
    nbnonu = 0
! ---  BOUCLE SUR LES NOEUDS PHYSIQUES :
!      -------------------------------
    do ire = 1, nbnore
        nbnonu = nbnonu+1
        zi(inum21+nbnonu) = zi(ioldn-1+ire)
    end do
!
! --- Commencer la numérotation des noeuds tardifs après tous les noeuds physiques
!
    nume_late_node = nb_node
!
! ---  BOUCLE SUR LES NOEUDS DE LAGRANGE :
!      ---------------------------------
    do ire = 1, nlag
        nbnonu = nbnonu+1
        nume_late_node = nume_late_node+1
        zi(inum21+nbnonu) = nume_late_node
    end do
!
    if (nbnonu .ne. (nbnore+nlag)) then
        call utmess('F', 'ASSEMBLA_28')
    end if
!
! --- CALCUL DES DESCRIPTEURS DES PRNO
!     ================================

! --- DETERMINATION DES .QRNM ET DES .QRNS POUR CHAQUE LIGREL :
!     -------------------------------------------------------
    do ili = 2, nlili
        call jenuno(jexnum(lili, ili), nomli)
        call creprn(nomli, modeLoc, 'V', nomli(1:19)//'.QRNM', nomli(1:19)//'.QRNS')
    end do

! --- 1ERE ETAPE : NOEUDS DU MAILLAGE (PHYSIQUES ET LAGRANGES)
! --- SI NUNOEL NOEUD DU MAILLAGE
! --- PRNO(1,NUNOEL,L+2)= "SIGMA"(.QRNM(MA(NUNOEL))(L))
!     -------------------------------------------------
    do ili = 2, nlili
        call jenuno(jexnum(lili, ili), nomli)
        call jeexin(nomli(1:19)//'.QRNM', iret)
        if (iret .eq. 0) goto 150
        call jeveuo(nomli(1:19)//'.QRNM', 'L', iprnm)

        do i = 1, nbnore
            nunoel = zi(ioldn-1+i)
            do l = 1, nec
                iec = zi(iprnm-1+nec*(nunoel-1)+l)
                zi(izzprn(1, nunoel, l+2)) = ior(zzprno(1, nunoel, l+2), iec)
            end do
        end do
150     continue
    end do
!
! --- 2EME ETAPE : NOEUDS SUPPLEMENTAIRES DES LIGRELS:
! --- SI NUNOEL NOEUD TARDIF DU LIGREL ILI NOMLI = LILI(ILI)
! --- PRNO(ILI,NUNOEL,L+2)= NOMLI.QRNS(NUNOEL)(L) :
!     -------------------------------------------
    do ili = 2, nlili
        call jenuno(jexnum(lili, ili), nomli)
        call jeveuo(nomli(1:19)//'.QRNM', 'L', iprnm)
        call jeveuo(nomli(1:19)//'.NBNO', 'L', vi=bid)
        if (bid(1) .gt. 0) call jeveuo(nomli(1:19)//'.QRNS', 'L', vi=qrns)

        do igr = 1, zzngel(ili)
            nel = zznelg(ili, igr)

            if (nel .ge. 0) then
                itypel = zzliel(ili, igr, nel+1)
                call jenuno(jexnum('&CATA.TE.NOMTE', itypel), nomte)
                icddlb = 0
                nddlb = 0
            end if

            do j = 1, nel
                numa = zzliel(ili, igr, j)
                if (numa .lt. 0) then
                    numa = -numa
                    do k = 1, zznsup(ili, numa)
                        nunoel = zznema(ili, numa, k)

                        if (nunoel .gt. 0) then
                            do l = 1, nec
                                iec = zi(iprnm+nec*(nunoel-1)+l-1)
                                zi(izzprn(ilim, nunoel, l+2)) = ior( &
                                                                zzprno(ilim, nunoel, l+2), iec)
                            end do
                            nunoel_save = nunoel
                        else
                            nunoel = -nunoel

!                 -- CALCUL DU NUMERO DE LA CMP ASSO
                            if (icddlb .eq. 0) then
                                ASSERT(gran_name .eq. nomte(3:8))
                                nomcmp = nomte(10:16)

!                   "GLUTE" POUR TEMP_MIL, TEMP_INF, TEMP_SUP :
                                if (nomcmp .eq. 'TEMP_MI') then
                                    nomcmp = 'TEMP_MIL'
                                else if (nomcmp .eq. 'TEMP_IN') then
                                    nomcmp = 'TEMP_INF'
                                else if (nomcmp .eq. 'TEMP_SU') then
                                    nomcmp = 'TEMP_SUP'
                                end if

                                call jeveuo(jexnom('&CATA.GD.NOMCMP', gran_name), 'L', idnocm)
                                call jelira(jexnom('&CATA.GD.NOMCMP', gran_name), 'LONMAX', nbcmp, &
                                            kbid)
                                nddlb = indik8(zk8(idnocm), nomcmp, 1, nbcmp)
                                ASSERT(nddlb .ne. 0)
                                icddlb = 1
                            end if

                            do l = 1, nec
                                iec = qrns(1+nec*(nunoel-1)+l-1)
                                zi(izzprn(ili, nunoel, l+2)) = ior(zzprno(ili, nunoel, l+2), iec)
                            end do

                            ilag = zi(inuno2+ili-1)+nunoel-1
                            ilag = ilag-nb_node
                            zi(iddlag+2*(ilag-1)+1) = -zi(iddlag+2*(ilag-1)+1)*nddlb
                        end if
                    end do
                end if
            end do
        end do
    end do
!
! --- CALCUL DES ADRESSES DANS LES PRNO
!     =================================
    iad = 1
    do i = 1, n
        call nuno1(i, ili, nunoel, n, inum21, inuno2, nlili)
        if (ili .gt. 0) then
            nddl1 = zzprno(ili, nunoel, 2)

            if (nddl1 .eq. 0) then
                nddl1 = nddl(ili, nunoel, nec, idprn1, idprn2)
!
                zi(izzprn(ili, nunoel, 2)) = nddl1
            end if
            zi(izzprn(ili, nunoel, 1)) = iad
            iad = iad+nddl1
        end if
    end do
!
! - Total number of dof (physical and NOT physical)
!
    nb_dof = iad-1
!
! - Number of dof in identity relations
!
    nb_iden_dof = nb_iden_term-nb_iden_rela
!
! - Create NEQU object
!
    call jedetr(nequ)
    call wkvect(nequ, base(2:2)//' V I', 2, vi=p_nequ)
!
! - Number of dof for computation (number of equations in system)
!
    nb_equa = nb_dof-nb_iden_dof
    p_nequ(1) = nb_equa
!
! - Total number of dof (with identity relations)
!
    p_nequ(2) = nb_dof
!
    if (niv .ge. 1 .and. .not. lparallel_mesh .and. modeLoc .eq. ' ') then

! ---   CALCUL DE NNO : NOMBRE DE NOEUDS DU MAILLAGE PORTEURS DE DDLS :
!       ----------------------------------------------------------------
        nno = 0
        call jeveuo(jexnum(prno, 1), 'L', jprno)
        do ino = 1, nb_node
            if (zi(jprno-1+(ino-1)*(2+nec)+2) .gt. 0) nno = nno+1
        end do
        vali(1) = nb_dof
        vali(2) = nb_dof-nlag
        vali(3) = nno
        vali(4) = nlag
        vali(5) = nlag
        call utmess('I', 'FACTOR_1', ni=5, vali=vali)
    end if

    call wkvect(refn, base(2:2)//' V K24', 5, idref)
    zk24(idref) = mesh
    zk24(idref+1) = gran_name
    zk24(idref+2) = modelZ
!
! - Create NUEQ object
!
    call jedetr(nueq)
    call wkvect(nueq, base(2:2)//' V I', nb_dof, ianueq)
!
! - Set NUEQ object
!
    call nunueq(mesh, nume_equa, nb_dof, igds, idenRela)
!
! - Create DEEQ object
!
    call jedetr(deeq)
    call wkvect(deeq, base(2:2)//' V I', 2*nb_equa, jdeeq)
!
! - Create DELG object
!
    call jedetr(delg)
    call wkvect(delg, base(2:2)//' V I', nb_equa, jdelg)
!
! - Set DEEQ and DELG objects with non-physical nodes
!
    call nudeeq_lag1(mesh, nb_node_mesh, nb_node_subs, numeDof, nb_equa, &
                     igds, iddlag)

! --- DESTRUCTION DES .QRNM ET DES .QRNS DE CHAQUE LIGREL :
!     ---------------------------------------------------
    do ili = 1, nlili
        call jenuno(jexnum(lili, ili), nomli)
        call jedetr(nomli(1:19)//'.QRNM')
        call jedetr(nomli(1:19)//'.QRNS')
    end do

    call jedetr(exi1)
    call jedetr(num2)
    call jedetr(num21)
    call jedetr(nuno)
    call jedetr(nnli)
    call jedetr(dsclag)

    call jedetr(newn)
    call jedetr(oldn)

!   call jedema() FORBIDDEN !
!
end subroutine
