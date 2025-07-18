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
! person_in_charge: mickael.abbas at edf.fr
!
module calcul_module
    implicit none

#include "asterf_types.h"

!    Ce module contient quelques variables "globales" utilisees
!    par la routine calcul.F90 ou pour les variables de commande

!======================================================================
    character(len=16) :: ca_option_, ca_nomte_, ca_nomtm_
    character(len=19) :: ca_ligrel_
!         ca_option_ : option a calculer
!         ca_ligrel_ : ligrel sur lequel il faut calculer l'option
!         ca_nomte_  : type_element courant
!         ca_nomtm_  : type_maille associe au type_element courant

!======================================================================
    integer(kind=8) :: ca_igd_, ca_nec_, ca_ncmpmx_, ca_iachin_, ca_iachlo_, ca_iichin_
    integer(kind=8) :: ca_ilchlo_, ca_itypgd_, ca_ianueq_, ca_lprno_
    character(len=8) :: ca_typegd_
!     ca_igd_ : numero de la grandeur associee au champ a extraire
!     ca_nec_ : nombre d'entiers codes de ca_igd_
!     ca_ncmpmx_: nombre max de cmps pour ca_igd_
!     ca_iachin_: adresse jeveux de chin//'.VALE'
!     ca_iachlo_: adresse jeveux de chloc//'.VALE'
!     ca_ilchlo_: adresse jeveux de chloc//'.EXIS'
!     ca_iichin_: numero du champ chin dans la liste lchin.
!     ca_ianueq_: adresse de l'objet .nueq du nume_equa associe eventuelle
!            -ment au champ chin. (si ca_lprno_=1).
!     ca_lprno_ : 1-> l'objet .nueq est a prendre en compte
!                 (cham_no a nume_equa)
!             0-> l'objet .nueq n'est pas a prendre en compte
!                 (cham_no a representation constante ou autre champ)
!     ca_typegd_: type scalaire de la grandeur ca_igd_ : 'R', 'I', 'K8', ...
!     ca_itypgd_: type scalaire de la grandeur ca_igd_ : 1,2,...
!             sert a eviter des comparaisons de chaines
!             la convention est ecrite dans extra1

!======================================================================
 integer(kind=8) :: ca_iaoptt_, ca_lgco_, ca_iaopmo_, ca_ilopmo_, ca_iaopno_, ca_ilopno_, ca_iaopds_
    integer(kind=8) :: ca_iaoppa_, ca_npario_, ca_nparin_, ca_iamloc_, ca_ilmloc_, ca_iadsgd_
!     ca_iaoptt_ : adresse de l'objet du catalogue : '&CATA.TE.OPTTE'
!     ca_lgco_   : longueur d'une colonne de '&CATA.TE.OPTTE'
!              ( nombre total d'options possibles du catalogue)
!     ca_iaopmo_ : adresse de '&CATA.TE.OPTMOD'
!     ca_ilopmo_ : adresse du pt_long de '&CATA.TE.OPTMOD'
!     ca_iaopno_ : adresse de '&CATA.TE.OPTNOM'
!     ca_ilopno_ : adresse du pt_long de '&CATA.TE.OPTNOM'
!     ca_iaopds_ : adresse de '&CATA.OP.DESCOPT(OPT)'
!     ca_iaoppa_ : adresse de '&CATA.OP.OPTPARA(OPT)'
!     ca_npario_ : longueur de '&CATA.OP.OPTPARA(OPT)'
!              = nb_param "in" + nb_param "out"
!     ca_nparin_ : nombre de parametres "in" pour l'ca_option_ opt.
!     ce nombre permet de savoir si un parametre est "in" ou "out"
!              (ipar <= ca_nparin_) <=> (ipar est "in")
!     ca_iamloc_ : adresse de '&CATA.TE.MODELOC'
!     ca_ilmloc_ : adresse du pt_long de '&CATA.TE.MODELOC'
!     ca_iadsgd_ : adresse de '&CATA.GD.DESCRIGD'

!======================================================================
    integer(kind=8) :: ca_iamaco_, ca_ilmaco_, ca_iamsco_, ca_ilmsco_, ca_ialiel_, ca_illiel_
!     ca_iamaco_  : adresse de la connectivite du maillage
!     ca_ilmaco_  : adresse du pointeur de longueur de ca_iamaco_
!     ca_iamsco_  : adresse de la connectivite des mailles suppl. d'1 ligrel
!     ca_ilmsco_  : adresse du pointeur de longueur de ca_iamsco_
!     ca_ialiel_  : adresse de l'objet '.liel' du ligrel.
!     ca_illiel_  : adresse du pointeur de longueur de '.liel'.

!======================================================================
    integer(kind=8), parameter :: ca_iachid_ = 12
    integer(kind=8) :: ca_iachii_, ca_iachik_, ca_iachix_
!     ca_iachii_ : adresse de '&&CALCUL.LCHIN_I'
!     ca_iachik_ : adresse de '&&CALCUL.LCHIN_K8'
!     ca_iachix_ : adresse de '&&CALCUL.LCHIN_EXI'

!     '&&CALCUL.LCHIN_EXI' ::= v(l)    (dim = nin)
!             v(1) :  .false.    : le champ parametre n'existe pas.

!     '&&CALCUL.LCHIN_K8'  ::= v(k8)    (dim = 2*nin)
!             v(1) :  type_champ : 'chno','cart','chml' ou 'resl'.
!             v(2) :  type_gd    : 'c', 'r', 'i', 'k8', ...

!     '&&CALCUL.LCHIN_I'  ::= v(i)     (dim = ca_iachid_*nin)
!             v(1) :  ca_igd_   grandeur associee a lchin(i)
!             v(2) :  ca_nec_   nombre d'entiers codes
!             v(3) :  ca_ncmpmx_ nombre max de cmp pour ca_igd_
!                           jelira(jexnum('&CATA.GD.NOMCMP', ca_igd_), 'LONMAX', ncmpmx)
!             v(4) :  iadesc adresse de .desc  (ou .celd)
!             v(5) :  iavale
!                     si cham_no ou carte : adresse de .vale
!                     si cham_elem        : adresse de .celv
!             v(6) :  iaptma adresse de .ptma (pour 1 carte)
!             v(7) :  iaptms adresse de .ptms (pour 1 carte)
!             v(8) :  iaprn1 adresse du prno($mailla) (pour 1 cham_no)
!             v(9) :  iaprn2 adresse du prno(ligrel)  (pour 1 cham_no)
!             v(10):  ca_ianueq_ adresse    .nueq         (pour 1 cham_no)
!             v(11):  ca_lprno_  (dit si ca_ianueq_ est utilise pour 1 cham_no)
!             v(12):  adresse des cmp pour ca_igd_
!                           jeveuo(jexnum('&CATA.GD.NOMCMP', ca_igd_), 'L', v(12) )

!======================================================================
    integer(kind=8) :: ca_ianoop_, ca_ianote_, ca_nbobtr_, ca_iaobtr_, ca_nbobmx_
!     ca_ianoop_ : adresse dans zk16 de '&&CALCUL.NOMOP' v(k16)
!          v(iop) --> nom de l'option iop
!     ca_ianote_ : adresse dans zk16 de '&&CALCUL.NOMTE' v(k16)
!          v(ite) --> nom du type_element ite
!          ca_nbobtr_ : nombre d'objets de travail '&&CALCUL....' qui
!                   devront etre detruits a la fin de calcul.
!          ca_iaobtr_ : adresse dans zk24 de l'objet '&&CALCUL.OBJETS_TRAV'
!          ca_nbobmx_ : longueur de l'objet '&&CALCUL.OBJETS_TRAV'

!======================================================================
    integer(kind=8) :: ca_nbgr_, ca_igr_, ca_nbelgr_, ca_nbelmx_, ca_jcteat_
    integer(kind=8) :: ca_lcteat_, ca_iawloc_, ca_iawlo2_, ca_iawtyp_

!     ca_nbgr_   : nombre de grel du ligrel
!     ca_igr_    : numero du grel courant
!     ca_nbelgr_ : nombre d'elements dans le grel ca_igr_
!     ca_nbelmx_ : nombre maximum d'elements dans un grel du ligrel

!     ca_jcteat_ : adresse dans zk16 de l'objet '&CATA.TE.CTE_ATTR(CA_NOMTE_)'
!     ca_lcteat_ : longueur de l'objet '&CATA.TE.CTE_ATTR(CA_NOMTE_)'
!       remarque : si ca_nomte_ n'a pas d'attribut : ca_jcteat_=ca_lcteat_=0

!  ca_iawloc_ : adresse dans zi de '&&CALCUL.IA_CHLOC' v(i)
!           cet objet contient des informations sur les champs locaux
!   v(3*(ipar-1)+1) : ca_iachlo_
!      adresse du champ_local '&&CALCUL.//NOMPAR(ipar)
!      =-1 <=> / le champ "in" n'existe pas :
!                 / nompar n'appartient pas a lpain
!                 / chin//'.desc' (ou .celd) n'existe pas
!              / nompar n'appartient pas a lpaout
!      =-2 <=> aucun type_elem du ligrel ne declare nompar

!   v(3*(ipar-1)+2) : ca_ilchlo_ :
!      adresse d'un vecteur de booleens ( // champ_local)
!      de nom : '&&CALCUL.//NOMPAR(ipar)//'.EXIS'
!      si ca_ilchlo_ = -1 :
!          =>  le champ local est "out"
!              et/ou le champ global n'existe pas
!              et/ou le parametre n'est pas utilise
!   v(3*(ipar-1)+3) : ich : numero du champ associe au parametre.
!      i.e : indice dans lchin (ou lchout selon le cas)
!      ich = 0 s'il n'y a pas de champ associe a ipar

!  ca_iawlo2_ : adresse dans zi de '&&CALCUL.IA_CHLO2' v(i)
!           cet objet contient des informations sur les champs locaux
!   v(5*(ca_nbgr_*(ipar-1)+ca_igr_-1)+1):
!      mode local pour (ipar,ca_igr_)
!   v(5*(ca_nbgr_*(ipar-1)+ca_igr_-1)+2):
!      longueur du champ_local pour 1 element (lue dans le catalogue).
!      cette longueur ne tient pas compte de nbspt et ncdyn.
!      =-1 <=> le parametre n'est pas utilise par le type_element
!   v(5*(ca_nbgr_*(ipar-1)+ca_igr_-1)+3):
!      nombre de points de discretis. du champ_local
!      0 si resuelem
!   v(5*(ca_nbgr_*(ipar-1)+ca_igr_-1)+4):
!      longueur du champ local pour le grel ca_igr_
!   v(5*(ca_nbgr_*(ipar-1)+ca_igr_-1)+5):
!      adresse du debut du grel dans le champ_local (= 1 si ca_calvoi_=0)

!   ca_iawtyp_ : adresse dans zk8 de '&&CALCUL.TYPE_SCA' v(k8)
!          v(ipar) --> type_scalaire du champ_local

!======================================================================
    integer(kind=8) :: ca_iachoi_, ca_iachok_
!     ca_iachoi_ : adresse de '&&CALCUL.LCHOU_I'
!     ca_iachok_ : adresse de '&&CALCUL.LCHOU_K8'

!     '&&CALCUL.LCHOU_K8'  ::= v(k8)    (dim = 2*nin)
!             v(1) :  type_champ : 'chml' ou 'resl'.
!             v(2) :  type_gd    : 'c', 'r'

!     '&&CALCUL.LCHOU_I'  ::= v(i)     (dim = 2*nout)
!         -- si chml :
!             v(1) :  adresse de l_chout(i).celd
!             v(2) :  adresse de l_chout(i).celv
!         -- si resl :
!             v(1) :  adresse de l_chout(i).desc

!======================================================================
    integer(kind=8) :: ca_iel_
!     ca_iel_    : numero de l'element courant (dans le grel ca_igr_)
!         (ca_iel_ est mis a jour par extrai,te0000,montee,...)

!======================================================================
    integer(kind=8) :: ca_nbobj_, ca_iainel_, ca_ininel_
!     ca_nbobj_  : nombre d'objets '&inel.xxxx' cree par l'initialisation
!              du type_elem
!     ca_ininel_ : adresse dans zk24 de l'objet '&&CALCUL.NOM_&INEL'
!              qui contient les noms des objets '&inel.xxxx'
!     ca_iainel_ : adresse dans zi de l'objet '&&CALCUL.IAD_&INEL'
!              qui contient les adresses des objets '&inel.xxxx'

!======================================================================
    integer(kind=8) :: ca_icaeli_, ca_icaelk_
!     ca_icaelk_ est l'adresse d'un vecteur de k24 contenant :
!       v(1) : nom du maillage  (k8)
!       v(2) : nom du ligrel    (k19)
!       v(3) : nom de la maille    (k8)
!       v(3+  1) : nom du 1er noeud de la maille  (k8)
!       v(3+  i) : ...
!       v(3+nbno) : nom du der noeud de la maille  (k8)
!       v(3+nbno+1) : nom du type_element (k16)
!       v(3+nbno+2) : nom de l'ca_option_     (k16)
!     ca_icaeli_ est l'adresse d'un vecteur de is contenant :
!       v(1) : numero de la maille
!       v(2) : nombre de noeuds de la maille (nbno)
!       v(2+   1) : numero du 1er noeud de la maille
!       v(2+nbno) : numero du der noeud de la maille
!       v(2+nbno +1) : numero du grel
!       v(2+nbno +2) : numero de l'element dans le grel

!======================================================================
    integer(kind=8) :: ca_nute_, ca_nuop_, ca_jnbelr_, ca_jnoelr_, ca_iactif_
    integer(kind=8) :: ca_jpnlfp_, ca_jnolfp_, ca_nblfpg_
!     ca_nute_ : numero du type_elem ca_nomte_.
!     ca_nuop_ : numero de l'option ca_option_.
!     ca_jnbelr_ : adresse dans zi  de '&CATA.TE.NBELREFE'
!     ca_jnoelr_ : adresse dans zk8 de '&CATA.TE.NOELREFE'
!     ca_iactif_ :
!                   1 : on est en "dessous" de la routine calcul
!                   2 : on est en "dessous" de la routine op0033 (CALC_POINT_MAT)
!                   3 : on est en "dessous" de la routine te0000
!                   0 : on n'est ni "sous" calcul ni "sous" op0033
!         Ce "booleen" (quaternaire !) permet par exemple d'arreter les utilisations des
!         routines qui ne doivent etre appellees que "sous" la routine te0000 :
!           jevech, tecach, tecael, elref1, ...
!         Il permet egalement une programmation differente de certains utilitaires :
!           utmess, lc0000, teattr, rcvarc, ...

!     ca_jpnlfp_ : adresse dans zk32 de '&CATA.TE.PNLOCFPG'
!     ca_jnolfp_ : adresse dans zi de '&CATA.TE.NOLOCFPG'
!     ca_nblfpg_ : dimension des objets pnlocfpg et nolocfpg
!======================================================================
    integer(kind=8) :: ca_caindz_(512), ca_capoiz_
!======================================================================
    integer(kind=8) :: ca_nbsav_
!======================================================================
    integer(kind=8) :: ca_nbcvrc_, ca_jvcnom_, ca_jvcfon_, ca_jvcval_
!     ca_nbcvrc_ : nombre de cvrc (variable de commande scalaire)
!     ca_jvcnom_ : adresse dans ZK8 des noms des cvrc
!     ca_jvcfon_ : adresse dans zk8 de '&&OP0033.TVCFON' (CALC_POINT_MAT)
!     ca_jvcval_ : adresse dans ZR  de '&&OP0033.TVCVAL' (CALC_POINT_MAT)
!======================================================================
    integer(kind=8) :: ca_ctempl_
    real(kind=8) :: ca_ctempr_, ca_ctempm_, ca_ctempp_
!     ca_ctempl_ : 1 if temperature is coupled variable (non external state variable)
!     ca_ctempr_ : for reference temperature when coupled variable (non external state variable)
!     ca_ctempm_ : for previous temperature when coupled variable (non external state variable)
!     ca_ctempp_ : for end temperature when coupled variable (non external state variable)
!======================================================================
    integer(kind=8) :: ca_cpcapl_
    real(kind=8) :: ca_cpcapm_, ca_cpcapp_
!     ca_cpcapl_ : 1 if capillar pressure is coupled variable (non external state variable)
!     ca_cpcapm_ : for previous pressure when coupled variable (non external state variable)
!     ca_cpcapp_ : for end pressure when coupled variable (non external state variable)
!======================================================================
    integer(kind=8) :: ca_nfpgmx_
    parameter(ca_nfpgmx_=10)
  integer(kind=8) :: ca_nfpg_, ca_jfpgl_, ca_decala_(ca_nfpgmx_), ca_km_, ca_kp_, ca_kr_, ca_iredec_
!     ca_nfpg_   : nombre de familles de la famille liste "mater"
!     ca_jfpgl_  : adresse dans zk8 de la liste des noms des familles
!     ca_decala_ : tableau des decalage des numeros des pg :
!         ca_decala_(1) : 0
!         ca_decala_(k) : nombre cumule des pg des familles (1:k-1)
!     ca_km_,ca_kp_,ca_kr_ : sauvegarde des numeros d'elements utilises dans
!                rcvarc. ils evitent des appels couteux a tecach.
!     ca_iredec_ : 0 : on n'est pas "sous" redece.f
!               1 : on est  "sous" redece.F90  (voir carr01 ci-dessous)
!======================================================================
    real(kind=8) :: ca_timed1_, ca_timef1_, ca_td1_, ca_tf1_
!     attention : il n'est pas valable si l'on n'est pas "sous" redece
!     ca_timed1_  : valeur de l'instant "-" du "gros" pas de temps
!     ca_timef1_  : valeur de l'instant "+" du "gros" pas de temps
!     td      : valeur de l'instant "-" du "petit" pas de temps
!     tf      : valeur de l'instant "+" du "petit" pas de temps
!======================================================================
    integer(kind=8) :: ca_calvoi_, ca_jrepe_, ca_jptvoi_, ca_jelvoi_

!     ca_calvoi_  = 1 : dans une routine te00ij, on veut pouvoir acceder
!                   aux champs des elements voisins (tecac2.f)
!     ca_calvoi_  = 0 : sinon

!     les 3 adresses suivantes sont remplies si ca_calvoi_=1
!      * ca_jrepe_   : adresse jeveux de ligrel.repe
!      * ca_jptvoi_  : adresse jeveux de maillage.vge.ptvois
!      * ca_jelvoi_  : adresse jeveux de maillage.vge.elvois
!======================================================================
!   Gestion du parallelisme :
!   =========================
    aster_logical :: ca_ldist_, ca_ldgrel_, ca_lparal_
    integer(kind=8) :: ca_nbproc_, ca_rang_
    integer(kind=8), pointer :: ca_numsd_(:) => null()
    aster_logical, pointer :: ca_paral_(:) => null()

!      * ca_nbproc_  : nombre de processeurs
!      * ca_rang_    : numero du processeur (0:nbproc-1)
!      * ca_ldist_   : .true. => le calcul est "distribue" :
!          Un element n'est calcule que par un processeur.
!      * ca_ldgrel_   : .true.  => le calcul est "distribue par grel"
!          Tous les elements d'un grel sont traites par le meme processeur.
!          Le processeur "rang" calcule le grel de numero igrel si :
!          mod(igrel,nbproc) .eq. rang
!      * ca_lparal_   : .true.  => le calcul est "distribue par element"
!      * ca_paral_    : pointeur sur l'objet '&&CALCUL.PARALLELE'
!          Le processeur "rang" calcule l'element iel si :
!          ca_paral_(iel)=.true.
!      * ca_numsd_   : ca_numsd_(ima) -> rang
!          rang est le numero du processeur qui doit traiter l'element porte
!          par la maille ima.

!      Toutes ces variables sont calculees dans la routine debca1.F90
!      Puis, calcul.F90 modifie ca_paral_(:) entre 2 grels.

!======================================================================
contains

!>  Initialise iactif (appel nécessaire en cas de sortie prématurée de calcul)
    subroutine calcul_init()
        ca_iactif_ = 0
    end subroutine

!>  Dit si on se trouve sous *calcul*  (ou SIMU_POINT_MAT)
    function calcul_status()
        integer(kind=8) :: calcul_status
        calcul_status = ca_iactif_
    end function calcul_status

end module
