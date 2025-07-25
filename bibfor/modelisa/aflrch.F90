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
subroutine aflrch(lisrez, chargz, type_liai, elim, detr_lisrez, l_preallocz)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/assert.h"
#include "asterfort/cragch.h"
#include "asterfort/craglc.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/impre2.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/nocart.h"
#include "asterfort/noligr.h"
#include "asterfort/ordlrl.h"
#include "asterfort/utmess.h"
#include "asterfort/nbec.h"
#include "asterfort/jexnum.h"
#include "asterfort/jexatr.h"
#include "asterfort/getvtx.h"
#include "asterfort/agdual.h"
#include "asterfort/char8_to_int.h"
!
    character(len=*), intent(in) :: lisrez
    character(len=*), intent(in) :: chargz
    character(len=*), intent(in) :: type_liai
    character(len=*), intent(in), optional :: elim
    aster_logical, intent(in), optional :: detr_lisrez
    aster_logical, intent(in), optional :: l_preallocz
!
! -------------------------------------------------------
!  affectation de l'objet de type  liste_rela et de nom
!  lisrel a l'objet de type charge et de nom charge
! -------------------------------------------------------
!     si la sd_liste_rela n'existe pas, on ne fait rien.
!     la charge doit exister au prealable
!
!     la sd_liste_rela est detruite a la fin de la routine
! -------------------------------------------------------
!  lisrel        - in    - k24  - : nom de la sd liste_rela
!                - jxvar -      -
! -------------------------------------------------------
!  charge        - in    - k8   - : nom de la sd charge
!                - jxvar -      -
! -----------------------------------------------------------
!  type_liai     - in    - k8   - : 'LIN'/'NLIN'
!    'NLIN' : les relations sont lineaires si TYPE_CHARGE='SUIV'
! -----------------------------------------------------------
!  (f) elim  : /'OUI' : on veut eliminer les doublons
!              /'NON' : on ne veut pas eliminer les doublons
! -----------------------------------------------------------
!  (f) detr_lisrez : .true. : on detruit la sd liste de relation
! -----------------------------------------------------------
    character(len=19) :: lisrel
    character(len=8) :: charge
    integer(kind=8) :: nmocl
    parameter(nmocl=300)
    complex(kind=8) :: betac
    character(len=3) :: klag2
    character(len=4) :: typval, typcoe
    character(len=24) :: valk(2)
    character(len=7) :: typcha
    character(len=19) :: betaf
    character(len=8) :: mod, nomgd, nomnoe, kelim
    character(len=8) :: noma, cmp, nomcmp(nmocl), ctype1, ctype2
    character(len=9) :: nomte
    character(len=19) :: ca1, ca2
    character(len=19) :: ligrmo, ligrch
    integer(kind=8) :: ntypel(nmocl)
    real(kind=8) :: beta
    integer(kind=8) :: i, icmp, iddl, idecal, ifm, igrel, gd1, gd2, nbRet
    integer(kind=8) :: in, indsur, inema, inema0, ino, inom, ipntrl, irela
    integer(kind=8) :: iret, j, jprnm, jrlbe, jrlco, ier
    integer(kind=8) :: jrlcof, jrldd, jrlno, idnoeu, jrlpo
    integer(kind=8) ::    jvale1, jvale2, jvalv1, jvalv2, kddl, nec1, nec2
    integer(kind=8) ::    ncmpmx1, ncmpmx2, jnocmp1, jnocmp2, jnoma1, jnoma2
    integer(kind=8) ::    jnoli1, jnoli2, jdesc1, jdesc2, jncmp1, jncmp2, nedit
    integer(kind=8) ::    jlima01, jlima02, jlimac1, jlimac2, lontav1, lontav2
    integer(kind=8) :: nbcmp, nec, nbnema, nbrela, nbteli, nbterm, nddla
    integer(kind=8) :: jliel0, jlielc, jnema0, jnemac, nbrela2, nbterm2
    character(len=3), parameter :: rapide = 'OUI'
    aster_logical :: detr_lisrel, l_lag1, l_prealloc, lcolle

    integer(kind=8) :: niv, numel, nunewm, iexi, jlgns
    character(len=8), pointer :: lgrf(:) => null()
    integer(kind=8), pointer :: rlnr(:) => null()
    character(len=8), pointer :: rltc(:) => null()
    integer(kind=8), pointer :: rlnt(:) => null()
    integer(kind=8), pointer :: rlsu(:) => null()
    integer(kind=8), pointer :: nbno(:) => null()
    character(len=8), pointer :: rltv(:) => null()
! --------- FIN  DECLARATIONS  VARIABLES LOCALES --------
    call jemarq()
    lisrel = lisrez
    call jeexin(lisrel//'.RLCO', iexi)
    if (iexi .eq. 0) goto 999
!
    charge = chargz
    call infniv(ifm, niv)
!
!
    call dismoi('NOM_MODELE', charge, 'CHARGE', repk=mod)
    ligrmo = mod(1:8)//'.MODELE'
    call jeveuo(ligrmo//'.LGRF', 'L', vk8=lgrf)
    noma = lgrf(1)
    lcolle = .false.
    call jeexin(noma//'.NOMNOE', ier)
    if (ier .ne. 0) then
        lcolle = .true.
    end if
    call dismoi('TYPE_CHARGE', charge, 'CHARGE', repk=typcha)
!
    l_lag1 = .false.
    call getvtx(' ', 'DOUBLE_LAGRANGE', scal=klag2, nbret=nbRet)
    if (nbRet .ne. 0) then
        if (klag2 .eq. 'NON') then
            l_lag1 = .true.
        end if
    end if

    if (typcha(1:4) .eq. 'MECA') then
        ligrch = charge//'.CHME.LIGRE'
        nomgd = 'DEPL_R'
        if (.not. l_lag1) then
            nomte = 'D_DEPL_R_'
        else
            nomte = 'D_DEPL_RS'
        end if
    else if (typcha(1:4) .eq. 'THER') then
        ligrch = charge//'.CHTH.LIGRE'
        nomgd = 'TEMP_R'
        if (.not. l_lag1) then
            nomte = 'D_TEMP_R_'
        else
            nomte = 'D_TEMP_RS'
        end if
    else if (typcha(1:4) .eq. 'ACOU') then
        ligrch = charge//'.CHAC.LIGRE'
        nomgd = 'PRES_C'
        if (.not. l_lag1) then
            nomte = 'D_PRES_C_'
        else
            ASSERT(.false.)
        end if
    end if

    if (type_liai .ne. 'LIN') then
        ASSERT(nomgd .eq. 'DEPL_R')
    end if

!   -- Rearrangement des relations de la liste de relations
!   -- par ordre de numero de noeud croissant
!   -- et suppression des relations redondantes en
!   -- appliquant le principe de surcharge

    kelim = 'OUI'
    if (present(elim)) then
        ASSERT(elim .eq. 'OUI' .or. elim .eq. 'NON')
        kelim = elim
    end if
    if (kelim .eq. 'OUI') then
        call ordlrl(charge, lisrel, nomgd)
    end if
!
    detr_lisrel = .true.
    if (present(detr_lisrez)) then
        detr_lisrel = detr_lisrez
    end if
!
    l_prealloc = .false.
    if (present(l_preallocz)) then
        l_prealloc = l_preallocz
    end if
!
    if (ligrch(12:13) .eq. 'TH') then
        ca1 = charge//'.CHTH.CMULT'
        ca2 = charge//'.CHTH.CIMPO'
    else if (ligrch(12:13) .eq. 'ME') then
        ca1 = charge//'.CHME.CMULT'
        ca2 = charge//'.CHME.CIMPO'
    else if (ligrch(12:13) .eq. 'AC') then
        ca1 = charge//'.CHAC.CMULT'
        ca2 = charge//'.CHAC.CIMPO'
    else
        ASSERT(.false.)
    end if
!
    call jeveuo(lisrel//'.RLTC', 'L', vk8=rltc)
    typcoe = rltc(1) (1:4)
!
! --- TYPE DES VALEURS AU SECOND MEMBRE DES RELATIONS
    call jeveuo(lisrel//'.RLTV', 'L', vk8=rltv)
    typval = rltv(1) (1:4)
!
! --- NOMBRE DE RELATIONS DE LA LISTE DE RELATIONS
    call jeveuo(lisrel//'.RLNR', 'L', vi=rlnr)
    nbrela = rlnr(1)
!
! --- NOMBRE TOTAL DE TERMES IMPLIQUES DANS LES RELATIONS
! --- DE LA LISTE DE RELATIONS (SERT AU REDIMENSIONNEMENT
! --- DU LIGREL DE CHARGE ET DES CARTES .CMULT ET .CIMPO
! --- DE LA CHARGE)
    call jeveuo(lisrel//'.RLPO', 'L', jrlpo)
    nbteli = zi(jrlpo+nbrela-1)
! --- SI ON N'A PAS AU PREALABLE ALLOUE LE LIGREL DE CHARGE
!     ET LES CARTES CMULT ET CIMPO C'EST LE MOMENT DE LE FAIRE
    if (.not. l_prealloc) then
!
! --- VERIFICATION DE L'ADEQUATION DE LA TAILLE DU LIGREL DE
! --- CHARGE A SON AFFECTATION PAR LES MAILLES TARDIVES DUES
! --- AUX RELATIONS LINEAIRES
! --- SI LE LIGREL DE CHARGE N'EXISTE PAS, ON LE CREE
        call craglc(nbteli, ligrch)
!
! --- VERIFICATION DE L'ADEQUATION DE LA TAILLE DES CARTES
! --- .CMULT ET .CIMPO DE LA CHARGE A LEUR AFFECTATION
! --- PAR LES MAILLES TARDIVES DUES AUX RELATIONS LINEAIRES
! --- SI LES CARTES .CMULT ET .CIMPO N'EXISTENT PAS, ON
! --- LES CREE
        call cragch(nbteli, typcoe, typval, ligrch)
    end if
!
!
!
    call jeveuo(jexnom('&CATA.GD.NOMCMP', nomgd), 'L', inom)
    call jelira(jexnom('&CATA.GD.NOMCMP', nomgd), 'LONMAX', nbcmp)
    nddla = nbcmp-1
    ASSERT(nddla .le. nmocl)
    ASSERT(nbcmp .le. nmocl)
    do i = 1, nbcmp
        nomcmp(i) = zk8(inom-1+i)
        call jenonu(jexnom('&CATA.TE.NOMTE', nomte//nomcmp(i) (1:7)), ntypel(i))
    end do
    call dismoi('NB_EC', nomgd, 'GRANDEUR', repi=nec)
!
!
    ASSERT(nec .le. 10)
    call jeveuo(ligrmo//'.PRNM', 'L', jprnm)

!   -- les cartes ca1 et ca2 doivent obligatoirement avoir ete
!   -- crees au prealable
!   -----------------------------------------------------------
    call jeexin(ca1//'.DESC', iret)
    ASSERT(iret .gt. 0)

    call jeveuo(ca1//'.DESC', 'E', jdesc1)
    call jeveuo(ca1//'.NOMA', 'E', jnoma1)
    call jeveuo(ca1//'.NOLI', 'E', jnoli1)
    call jeveuo(ca1//'.VALE', 'E', jvale1)
    call jeveuo(ca1//'.NCMP', 'E', jncmp1)
    call jeveuo(ca1//'.VALV', 'E', jvalv1)
    call jelira(ca1//'.VALV', 'TYPELONG', cval=ctype1)
    call jeveuo(ca1//'.LIMA', 'E', jlima01)
    call jeveuo(jexatr(ca1//'.LIMA', 'LONCUM'), 'E', jlimac1)
    call jelira(ca1//'.LIMA', 'LONT', lontav1)
    gd1 = zi(jdesc1-1+1)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd1), 'L', jnocmp1)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd1), 'LONMAX', ncmpmx1)
    nec1 = nbec(gd1)

    call jeveuo(ca2//'.DESC', 'E', jdesc2)
    call jeveuo(ca2//'.NOMA', 'E', jnoma2)
    call jeveuo(ca2//'.NOLI', 'E', jnoli2)
    call jeveuo(ca2//'.VALE', 'E', jvale2)
    call jeveuo(ca2//'.NCMP', 'E', jncmp2)
    call jeveuo(ca2//'.VALV', 'E', jvalv2)
    call jelira(ca2//'.VALV', 'TYPELONG', cval=ctype2)
    call jeveuo(ca2//'.LIMA', 'E', jlima02)
    call jeveuo(jexatr(ca2//'.LIMA', 'LONCUM'), 'E', jlimac2)
    call jelira(ca2//'.LIMA', 'LONT', lontav2)
    gd2 = zi(jdesc2-1+1)
    call jeveuo(jexnum('&CATA.GD.NOMCMP', gd2), 'L', jnocmp2)
    call jelira(jexnum('&CATA.GD.NOMCMP', gd2), 'LONMAX', ncmpmx2)
    nec2 = nbec(gd2)

    numel = 0
    call jeveuo(ligrch//'.NBNO', 'E', vi=nbno)
    call dismoi('NB_MA_SUP', ligrch, 'LIGREL', repi=inema)
    call dismoi('NB_GREL', ligrch, 'LIGREL', repi=igrel)
!
!
    call jeveuo(lisrel//'.RLCO', 'L', jrlco)
    call jeveuo(lisrel//'.RLDD', 'L', jrldd)
    call jeveuo(lisrel//'.RLNO', 'L', jrlno)
    call jeveuo(lisrel//'.RLNT', 'L', vi=rlnt)
    call jeveuo(lisrel//'.RLPO', 'L', jrlpo)
    call jeveuo(lisrel//'.RLSU', 'L', vi=rlsu)
    call jeveuo(lisrel//'.RLBE', 'L', jrlbe)
!
    call jeexin(ligrch//'.LGNS', iexi)
    if (iexi .gt. 0) then
        call jeveuo(ligrch//'.LGNS', 'E', jlgns)
    else
        jlgns = 1
    end if

    call jeveuo(ligrch//'.LIEL', 'E', jliel0)
    call jeveuo(jexatr(ligrch//'.LIEL', 'LONCUM'), 'E', jlielc)
    call jeveuo(ligrch//'.NEMA', 'E', jnema0)
    call jeveuo(jexatr(ligrch//'.NEMA', 'LONCUM'), 'E', jnemac)

    do irela = 1, nbrela
        indsur = rlsu(irela)
        if (indsur .ne. 0) goto 60
!
        ipntrl = zi(jrlpo+irela-1)
        nbterm = rlnt(irela)
        if (typval .eq. 'REEL') then
            beta = zr(jrlbe+irela-1)
        else if (typval .eq. 'COMP') then
            betac = zc(jrlbe+irela-1)
        else if (typval .eq. 'FONC') then
            betaf = zk24(jrlbe+irela-1) (1:19)
        else
            ASSERT(.false.)
        end if
        idecal = ipntrl-nbterm
        jrlcof = jrlco+idecal
        idnoeu = jrlno+idecal
        iddl = jrldd+idecal
!
        numel = 0
        inema0 = inema
        if (l_lag1) then
            nbno(1) = nbno(1)+1
        else
            nbno(1) = nbno(1)+2
        end if
        do ino = 1, nbterm
            nomnoe = zk8(idnoeu+ino-1)
            in = char8_to_int(nomnoe, lcolle, noma, "NOEUD")
!
            cmp = zk8(iddl+ino-1)
!
            icmp = indik8(nomcmp, cmp, 1, nbcmp)
            if (.not. exisdg(zi(jprnm-1+(in-1)*nec+1), icmp)) then
                valk(1) = cmp
                valk(2) = nomnoe
                call utmess('F', 'CHARGES2_31', nk=2, valk=valk)
            else
                do kddl = 1, nbcmp
                    if (cmp .eq. nomcmp(kddl)) then
                        numel = ntypel(kddl)
                        goto 30
                    end if
                end do
            end if
30          continue
!
            if (numel .ne. 0) then
                igrel = igrel+1
                call noligr(ligrch, igrel, numel, in, &
                            3, inema, nbno(1), jlgns, &
                            rapide=rapide, jliel0=jliel0, jlielc=jlielc, &
                            jnema0=jnema0, jnemac=jnemac, l_lag1=l_lag1)
            else
                call utmess('F', 'CHARGES2_33', sk=nomnoe)
            end if
        end do
!
!
!       --  stockage dans les cartes ca1 et ca2
!       -----------------------------------------
        nbnema = inema-inema0
        ASSERT(nbnema .eq. nbterm)
        zk8(jncmp1-1+1) = 'A1'
        zk8(jncmp2-1+1) = 'C'
        do j = 1, nbnema
            if (typcoe .eq. 'COMP') then
                zc(jvalv1) = zc(jrlcof-1+j)
            else
                zr(jvalv1) = zr(jrlcof-1+j)
            end if
            nunewm = -(inema0+j)
            call nocart(ca1, -3, 1, ligrel=ligrch, nma=1, limanu=[nunewm], &
                        rapide=rapide, jdesc=jdesc1, jnoma=jnoma1, &
                        jncmp=jncmp1, jnoli=jnoli1, &
                        jvale=jvale1, jvalv=jvalv1, jnocmp=jnocmp1, &
                        ncmpmx=ncmpmx1, nec=nec1, ctype=ctype1, &
                        jlima0=jlima01, jlimac=jlimac1, lontav=lontav1)
            if (j .lt. nbnema) then
                if (typval .eq. 'REEL') then
                    zr(jvalv2) = 0.d0
                else if (typval .eq. 'COMP') then
                    zc(jvalv2) = (0.0d0, 0.0d0)
                else if (typval .eq. 'FONC') then
                    zk24(jvalv2) = '&FOZERO'
                else
                    ASSERT(.false.)
                end if
            else
                if (typval .eq. 'REEL') then
                    zr(jvalv2) = beta
                else if (typval .eq. 'COMP') then
                    zc(jvalv2) = betac
                else if (typval .eq. 'FONC') then
                    zk24(jvalv2) = betaf
                else
                    ASSERT(.false.)
                end if
            end if
            call nocart(ca2, -3, 1, ligrel=ligrch, nma=1, limanu=[nunewm], &
                        rapide=rapide, jdesc=jdesc2, jnoma=jnoma2, &
                        jncmp=jncmp2, jnoli=jnoli2, &
                        jvale=jvale2, jvalv=jvalv2, jnocmp=jnocmp2, &
                        ncmpmx=ncmpmx2, nec=nec2, ctype=ctype2, &
                        jlima0=jlima02, jlimac=jlimac2, lontav=lontav2)
        end do
60      continue
    end do

!   -- A cause de l'argument rapide='OUI' pour les routines nocart et noligr,
!      Il faut faire des appels a jeecra('NUTIOC') :
    if (rapide .eq. 'OUI') then
        nedit = zi(jdesc1-1+3)
        ASSERT(nedit .eq. zi(jdesc2-1+3))
        call jeecra(ca1//'.LIMA', 'NUTIOC', ival=nedit)
        call jeecra(ca2//'.LIMA', 'NUTIOC', ival=nedit)
        call jeecra(ligrch//'.LIEL', 'NUTIOC', ival=igrel)
        call jeecra(ligrch//'.NEMA', 'NUTIOC', ival=inema)
    end if

!   -- Mise a jour de la SD char_dual :
    nbrela2 = 0
    nbterm2 = 0
    do irela = 1, nbrela
        indsur = rlsu(irela)
        if (indsur .eq. 0) then
            nbrela2 = nbrela2+1
            nbterm2 = nbterm2+rlnt(irela)
        end if
    end do

    if (nomgd .eq. 'DEPL_R') then
        call agdual(charge, type_liai, nbrela2, nbterm2)
    end if

!
!   -- impression des relations redondantes et donc supprimees :
!   ------------------------------------------------------------
    if (niv .ge. 2) then
        call utmess('I', 'CHARGES2_34')
        do irela = 1, nbrela
            indsur = rlsu(irela)
            if (indsur .eq. 1) then
                ipntrl = zi(jrlpo+irela-1)
                nbterm = rlnt(irela)
                if (typval .eq. 'REEL') then
                    beta = zr(jrlbe+irela-1)
                else if (typval .eq. 'COMP') then
                    betac = zc(jrlbe+irela-1)
                else if (typval .eq. 'FONC') then
                    betaf = zk24(jrlbe+irela-1) (1:19)
                else
                    ASSERT(.false.)
                end if
                idecal = ipntrl-nbterm
                jrlcof = jrlco+idecal
                idnoeu = jrlno+idecal
                iddl = jrldd+idecal
                call impre2(lisrel//'.RLCO', lisrel//'.RLDD', lisrel//'.RLNO', lisrel//'.RLBE', &
                            rlsu(irela), zi(jrlpo+irela-1), rlnt(irela), typcoe, &
                            typval, irela)
            end if
        end do
    end if
!
    if (detr_lisrel) then
        call jedetr(lisrel//'.RLCO')
        call jedetr(lisrel//'.RLDD')
        call jedetr(lisrel//'.RLNO')
        call jedetr(lisrel//'.RLBE')
        call jedetr(lisrel//'.RLNT')
        call jedetr(lisrel//'.RLPO')
        call jedetr(lisrel//'.RLSU')
        call jedetr(lisrel//'.RLNR')
        call jedetr(lisrel//'.RLTC')
        call jedetr(lisrel//'.RLTV')
    end if
!
999 continue
    call jedema()
end subroutine
