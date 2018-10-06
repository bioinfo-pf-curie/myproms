/*========================================================================*/
/*                  myProMS database v3.5.16                              */
/* MySQL script for myproms v3.6 database                                 */
/* Requires MySQL 5+                                                      */
/* Patrick Poullet    25/06/2018                                          */
/* Command:                                                               */
/* >mysql -u <DB_USER> -p -D <DB_NAME> -h <DB_HOST> < myproms_db.sql      */
/*========================================================================*/


/*==============================================================*/
/* Table : ANALYSIS                                             */
/*==============================================================*/
create table ANALYSIS
(
   ID_ANALYSIS          int not null auto_increment,
   ID_SAMPLE            int not null,
   NAME                 varchar(100),
   DES                  varchar(100),
   START_DATE           datetime,
   MS_TYPE              varchar(10),
   INSTRUMENT           varchar(30),
   DATA_FILE            varchar(100),
   FILE_FORMAT          varchar(15),
   WIFF_FILE            varchar(100),
   LAB_CODE             varchar(50),
   LABELING             varchar(50),
   DECOY                varchar(30),
   FALSE_DISCOVERY      varchar(20),
   MIN_SCORE            float,
   MAX_RANK             smallint,
   TAXONOMY             varchar(50),
   COMMENTS             varchar(250),
   VALID_STATUS         smallint,
   VALID_USER           varchar(20),
   LOWER_SCORES         smallint,
   VERIFIED_MG          smallint,
   FIRST_VALID_DATE     datetime,
   VALID_DATE           datetime,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   DISPLAY_POS          smallint,
   primary key (ID_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANALYSIS_DATABANK                                    */
/*==============================================================*/
create table ANALYSIS_DATABANK
(
   ID_DATABANK          int not null,
   ID_ANALYSIS          int not null,
   DB_RANK              smallint,
   primary key (ID_DATABANK, ID_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANALYSIS_MODIFICATION                                */
/*==============================================================*/
create table ANALYSIS_MODIFICATION
(
   ID_ANALYSIS          int not null,
   ID_MODIFICATION      int not null,
   MODIF_TYPE           varchar(1) not null,
   SPECIFICITY          text,
   primary key (ID_ANALYSIS, ID_MODIFICATION, MODIF_TYPE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANALYSIS_PROTEIN                                     */
/*==============================================================*/
create table ANALYSIS_PROTEIN
(
   ID_ANALYSIS          int not null,
   ID_PROTEIN           int not null,
   DB_RANK              smallint,
   CONF_LEVEL           int,
   SCORE                float,
   NUM_PEP              int,
   NUM_MATCH            int,
   PEP_COVERAGE         float,
   MATCH_GROUP          int,
   PEP_SPECIFICITY      float,
   VISIBILITY           smallint,
   primary key (ID_ANALYSIS, ID_PROTEIN)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANALYSIS_REFRT                                       */
/*==============================================================*/
create table ANALYSIS_REFRT
(
   ID_REFERENCE_RT      int not null,
   ID_ANALYSIS          int not null,
   NUM_PEP              int,
   DATA                 text,
   primary key (ID_REFERENCE_RT, ID_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANALYSIS_SWATH_LIB                                   */
/*==============================================================*/
create table ANALYSIS_SWATH_LIB
(
   ID_SWATH_LIB         int not null,
   ID_ANALYSIS          int not null,
   VERSION_NAME         varchar(20),
   primary key (ID_SWATH_LIB, ID_ANALYSIS)
);

/*==============================================================*/
/* Table : ANA_COMPARISON                                       */
/*==============================================================*/
create table ANA_COMPARISON
(
   ID_COMPARISON        int not null,
   ID_ANALYSIS          int not null,
   COMP_GROUP           smallint not null,
   ANA_POS              smallint,
   primary key (ID_COMPARISON, ID_ANALYSIS, COMP_GROUP)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANA_QUANTIFICATION                                   */
/*==============================================================*/
create table ANA_QUANTIFICATION
(
   ID_QUANTIFICATION    int not null,
   ID_ANALYSIS          int not null,
   QUANTIF_FILE         varchar(50),
   IS_REFERENCE         smallint,
   primary key (ID_QUANTIFICATION, ID_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ANNOTATIONSET                                        */
/*==============================================================*/
create table ANNOTATIONSET
(
   ID_ANNOTATIONSET     int not null,
   ID_EXPLORANALYSIS    int not null,
   NAME                 varchar(50),
   RANK                 smallint,
   ANNOT_TYPE           varchar(50),
   ANNOT_LIST           text,
   primary key (ID_ANNOTATIONSET)
)
engine = InnoDB;

/*==============================================================*/
/* Table : BIOSAMPLE                                            */
/*==============================================================*/
create table BIOSAMPLE
(
   ID_BIOSAMPLE         int not null auto_increment,
   ID_SPECIES           int not null,
   ID_REFBIOSAMPLE      int,
   NAME                 varchar(50),
   DES                  varchar(100),
   IS_REFERENCE         smallint,
   RECORD_DATE          datetime,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_BIOSAMPLE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : BIOSAMPLE_PROPERTY                                   */
/*==============================================================*/
create table BIOSAMPLE_PROPERTY
(
   ID_BIOSAMPLE         int not null,
   ID_PROPERTY          int not null,
   RANK                 smallint,
   PROPERTY_VALUE       text,
   primary key (ID_BIOSAMPLE, ID_PROPERTY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : CATEGORY                                             */
/*==============================================================*/
create table CATEGORY
(
   ID_CATEGORY          int not null,
   ID_CLASSIFICATION    int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   COMMENTS             varchar(250),
   LIST_TYPE            varchar(4),
   DISPLAY_POS          smallint,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_CATEGORY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : CATEGORY_PROTEIN                                     */
/*==============================================================*/
create table CATEGORY_PROTEIN
(
   ID_CATEGORY_PROTEIN  int not null auto_increment,
   ID_CATEGORY          int not null,
   ID_PROTEIN           int not null,
   primary key (ID_CATEGORY_PROTEIN)
)
engine = InnoDB;

/*==============================================================*/
/* Table : CAT_COMPARISON                                       */
/*==============================================================*/
create table CAT_COMPARISON
(
   ID_COMPARISON        int not null,
   ID_CATEGORY          int not null,
   COMP_GROUP           smallint,
   CAT_POS              smallint,
   primary key (ID_COMPARISON, ID_CATEGORY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : CLASSIFICATION                                       */
/*==============================================================*/
create table CLASSIFICATION
(
   ID_CLASSIFICATION    int not null,
   ID_PROJECT           int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   COMMENTS             varchar(250),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_CLASSIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : COMPARISON                                           */
/*==============================================================*/
create table COMPARISON
(
   ID_COMPARISON        int not null,
   ID_CATEGORY          int,
   ID_PROJECT           int not null,
   NAME                 varchar(50),
   COMMENTS             varchar(250),
   NUM_GROUPS           smallint,
   CAT_EXCLUSION        smallint,
   PEPTIDE_PARAMS       text,
   SORT_ORDER           varchar(15),
   AUTOCHECK_PARAMS     text,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_COMPARISON)
)
engine = InnoDB;

/*==============================================================*/
/* Table : DATABANK                                             */
/*==============================================================*/
create table DATABANK
(
   ID_DATABANK          int not null,
   ID_DBTYPE            int not null,
   NAME                 varchar(50),
   VERSION_NAME         varchar(20),
   DES                  varchar(100),
   VERSION_DATE         date,
   FASTA_FILE           varchar(100),
   NUM_ENTRY            int,
   IDENTIFIER_TYPE      varchar(50),
   DECOY_TAG            varchar(100),
   IS_CRAP              smallint,
   COMMENTS             varchar(250),
   USE_STATUS           varchar(3),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   ORGANISM             varchar(100),
   primary key (ID_DATABANK)
)
engine = InnoDB;

/*==============================================================*/
/* Table : DATABANK_SWATHLIB                                    */
/*==============================================================*/
create table DATABANK_SWATHLIB
(
   ID_DATABANK          int not null,
   ID_SWATH_LIB         int not null,
   DB_RANK              smallint,
   primary key (ID_DATABANK, ID_SWATH_LIB)
)
engine = InnoDB;

/*==============================================================*/
/* Table : DATABANK_TYPE                                        */
/*==============================================================*/
create table DATABANK_TYPE
(
   ID_DBTYPE            int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   PARSE_RULES          varchar(100),
   DEF_IDENT_TYPE       varchar(50),
   IDENT_MOD_STRING     varchar(30),
   PARSE_STRING         varchar(500),
   primary key (ID_DBTYPE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : DESIGN                                               */
/*==============================================================*/
create table DESIGN
(
   ID_DESIGN            int not null,
   ID_EXPERIMENT        int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_DESIGN)
)
engine = InnoDB;

/*==============================================================*/
/* Table : EXPCONDITION                                         */
/*==============================================================*/
create table EXPCONDITION
(
   ID_EXPCONDITION      int not null,
   ID_DESIGN            int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   DISPLAY_POS          smallint,
   primary key (ID_EXPCONDITION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : EXPCONDITION_QUANTIF                                 */
/*==============================================================*/
create table EXPCONDITION_QUANTIF
(
   ID_EXPCONDITION      int not null,
   ID_QUANTIFICATION    int not null,
   COND_FUNCTION        varchar(15),
   QUANTIF_ELEMENT      varchar(25),
   primary key (ID_EXPCONDITION, ID_QUANTIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : EXPERIMENT                                           */
/*==============================================================*/
create table EXPERIMENT
(
   ID_EXPERIMENT        int not null,
   ID_PROJECT           int not null,
   ID_SPECIES           int,
   NAME                 varchar(50),
   DES                  varchar(100),
   START_DATE           datetime,
   COMMENTS             varchar(250),
   DISPLAY_POS          smallint,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_EXPERIMENT)
)
engine = InnoDB;

/*==============================================================*/
/* Table : EXPLORANALYSIS                                       */
/*==============================================================*/
create table EXPLORANALYSIS
(
   ID_EXPLORANALYSIS    int not null,
   ID_CATEGORY          int,
   ID_EXPERIMENT        int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   ANA_TYPE             varchar(10),
   PARAM_LIST           text,
   FILTER_LIST          text,
   CAT_EXCLUSION        smallint,
   STATUS               smallint,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_EXPLORANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : EXPLORANA_ANA                                        */
/*==============================================================*/
create table EXPLORANA_ANA
(
   ID_EXPLORANALYSIS    int not null,
   ID_ANALYSIS          int not null,
   GROUP_POS            smallint,
   primary key (ID_EXPLORANALYSIS, ID_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : EXPLORANA_QUANTIF                                    */
/*==============================================================*/
create table EXPLORANA_QUANTIF
(
   ID_QUANTIFICATION    int not null,
   ID_EXPLORANALYSIS    int not null,
   TARGET_POS           smallint not null,
   primary key (ID_QUANTIFICATION, ID_EXPLORANALYSIS, TARGET_POS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : GEL2D                                                */
/*==============================================================*/
create table GEL2D
(
   ID_GEL2D             int not null,
   ID_EXPERIMENT        int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   IMAGE_FILE           varchar(50),
   COMMENTS             varchar(250),
   DISPLAY_POS          smallint,
   START_DATE           datetime,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_GEL2D)
)
engine = InnoDB;

/*==============================================================*/
/* Table : GOANA_ANALYSIS                                       */
/*==============================================================*/
create table GOANA_ANALYSIS
(
   ID_GOANALYSIS        int not null,
   ID_ANALYSIS          int not null,
   primary key (ID_GOANALYSIS, ID_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : GOANA_QUANTIFICATION                                 */
/*==============================================================*/
create table GOANA_QUANTIFICATION
(
   ID_GOANALYSIS        int not null,
   ID_QUANTIFICATION    int not null,
   RATIO                varchar(5),
   primary key (ID_GOANALYSIS, ID_QUANTIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : GOANNOTATION                                         */
/*==============================================================*/
create table GOANNOTATION
(
   ID_GOANNOTATION      int not null auto_increment,
   ID_SPECIES           int not null,
   ID_IDENTIFIER        int,
   NAME                 varchar(50),
   DES                  varchar(100),
   ANNOT_FILE           varchar(100),
   STATUS               smallint,
   VERSION_DATE         date,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_GOANNOTATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : GO_ANALYSIS                                          */
/*==============================================================*/
create table GO_ANALYSIS
(
   ID_GOANALYSIS        int not null auto_increment,
   ID_ONTOLOGY          int not null,
   ID_EXPERIMENT        int not null,
   ID_GOANNOTATION      int not null,
   ID_CATEGORY          int,
   ID_PARENT_GOANA      int,
   NAME                 varchar(50),
   DES                  varchar(100),
   GOA_TYPE             varchar(10),
   ASPECT               varchar(3),
   PARAM_STRG           text,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_GOANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : IDENTIFIER                                           */
/*==============================================================*/
create table IDENTIFIER
(
   ID_IDENTIFIER        int not null,
   ID_SPECIES           int,
   NAME                 varchar(50),
   CODE                 varchar(15),
   RESRC_NAME           varchar(50),
   RESRC_URL            varchar(250),
   TARGET               varchar(10),
   primary key (ID_IDENTIFIER)
)
engine = InnoDB;

/*==============================================================*/
/* Table : INSTRUMENT                                           */
/*==============================================================*/
create table INSTRUMENT
(
   ID_INSTRUMENT        int not null,
   NAME                 varchar(50),
   ALIAS                text,
   DES                  varchar(100),
   COMMENTS             varchar(250),
   USE_STATUS           varchar(3),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   DELTA_PARENT         varchar(5),
   DELTA_FRAGMENT       varchar(5),
   TOL_FRAGMENT         float,
   NR_LEVEL             float,
   RULES                varchar(500),
   primary key (ID_INSTRUMENT)
)
engine = InnoDB;

/*==============================================================*/
/* Table : MASTERPROT_IDENTIFIER                                */
/*==============================================================*/
create table MASTERPROT_IDENTIFIER
(
   ID_IDENTIFIER        int not null,
   ID_MASTER_PROTEIN    int not null,
   RANK                 smallint not null,
   VALUE                text,
   primary key (ID_IDENTIFIER, ID_MASTER_PROTEIN, RANK)
)
engine = InnoDB;

/*==============================================================*/
/* Table : MASTER_PROTEIN                                       */
/*==============================================================*/
create table MASTER_PROTEIN
(
   ID_MASTER_PROTEIN    int not null auto_increment,
   ID_SPECIES           int,
   PROT_DES             varchar(250),
   PROT_LENGTH          int,
   MW                   float,
   PROT_SEQ             text,
   UPDATE_DATE          datetime,
   primary key (ID_MASTER_PROTEIN)
)
engine = InnoDB;

/*==============================================================*/
/* Table : MODIFICATION                                         */
/*==============================================================*/
create table MODIFICATION
(
   ID_MODIFICATION      int not null auto_increment,
   UNIMOD_ACC           int,
   PSI_MS_NAME          text,
   INTERIM_NAME         varchar(50),
   DES                  varchar(100),
   SYNONYMES            text,
   COMPOSITION          text,
   MONO_MASS            float,
   AVGE_MASS            float,
   SPECIFICITY          text,
   DISPLAY_CODE         varchar(5),
   DISPLAY_COLOR        varchar(7),
   IS_SUBST             smallint,
   IS_LABEL             smallint,
   PEAKVIEW_CODE        varchar(10),
   VALID_STATUS         smallint,
   primary key (ID_MODIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : MODIFIED_RESIDUE                                     */
/*==============================================================*/
create table MODIFIED_RESIDUE
(
   ID_MODIF_RES         int not null auto_increment,
   ID_QUANTIFICATION    int not null,
   RESIDUE              char(1),
   POSITION             int,
   primary key (ID_MODIF_RES)
)
engine = InnoDB;

/*==============================================================*/
/* Table : MODIFICATION_SITE                                    */
/*==============================================================*/
create table MODIFICATION_SITE
(
   ID_MODIFICATION_SITE int not null auto_increment,
   ID_MODIFICATION      int not null,
   ID_CATEGORY_PROTEIN  int not null,
   ID_CATEGORY          int not null,
   RESIDUE              char(1),
   POSITION             int,
   primary key (ID_MODIFICATION_SITE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : OBSERVATION                                          */
/*==============================================================*/
create table OBSERVATION
(
   ID_OBSERVATION       int not null auto_increment,
   ID_BIOSAMPLE         int,
   ID_ANALYSIS          int not null,
   TARGET_POS           smallint,
   primary key (ID_OBSERVATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : OBS_EXPCONDITION                                     */
/*==============================================================*/
create table OBS_EXPCONDITION
(
   ID_EXPCONDITION      int not null,
   ID_OBSERVATION       int not null,
   FRACTION_GROUP       smallint,
   TECH_REP_GROUP       smallint,
   primary key (ID_EXPCONDITION, ID_OBSERVATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : OBS_MODIFICATION                                     */
/*==============================================================*/
create table OBS_MODIFICATION
(
   ID_OBSERVATION       int not null,
   ID_MODIFICATION      int not null,
   primary key (ID_OBSERVATION, ID_MODIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : ONTOLOGY                                             */
/*==============================================================*/
create table ONTOLOGY
(
   ID_ONTOLOGY          int not null auto_increment,
   NAME                 varchar(50),
   OBO_FILE             varchar(100),
   DATA_VERSION         varchar(15),
   STATUS               smallint,
   VERSION_DATE         date,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_ONTOLOGY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PARENT_EXPLORANALYSIS                                */
/*==============================================================*/
create table PARENT_EXPLORANALYSIS
(
   ID_EXPLORANALYSIS    int not null,
   ID_PARENT_EXPLORANALYSIS int not null,
   DISPLAY_POS          smallint,
   primary key (ID_EXPLORANALYSIS, ID_PARENT_EXPLORANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PARENT_QUANTIFICATION                                */
/*==============================================================*/
create table PARENT_QUANTIFICATION
(
   ID_PARENT_QUANTIFICATION int not null,
   ID_QUANTIFICATION    int not null,
   PAR_FUNCTION         varchar(15),
   primary key (ID_PARENT_QUANTIFICATION, ID_QUANTIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PARENT_SWATH_LIB                                     */
/*==============================================================*/
create table PARENT_SWATH_LIB
(
   ID_SWATH_LIB         int not null,
   ID_PARENT_SWATH_LIB  int not null,
   VERSION_NAME         varchar(20),
   primary key (ID_SWATH_LIB, ID_PARENT_SWATH_LIB)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PATHWAYANA_ANALYSIS                                  */
/*==============================================================*/
create table PATHWAYANA_ANALYSIS
(
   ID_ANALYSIS          int not null,
   ID_PATHWAY_ANALYSIS  int not null,
   primary key (ID_ANALYSIS, ID_PATHWAY_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PATHWAYANA_QUANTIFICATION                            */
/*==============================================================*/
create table PATHWAYANA_QUANTIFICATION
(
   ID_QUANTIFICATION    int not null,
   ID_PATHWAY_ANALYSIS  int not null,
   RATIO                varchar(5),
   primary key (ID_QUANTIFICATION, ID_PATHWAY_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PATHWAY_ANALYSIS                                     */
/*==============================================================*/
create table PATHWAY_ANALYSIS
(
   ID_PATHWAY_ANALYSIS  int not null auto_increment,
   ID_PARENT_PATHWAYANA int,
   ID_EXPERIMENT        int not null,
   ID_CATEGORY          int,
   NAME                 varchar(50),
   DES                  varchar(100),
   PARAM_STRG           text,
   STATUS               smallint,
   RECORD_DATE          datetime,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_PATHWAY_ANALYSIS)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PEPTIDE                                              */
/*==============================================================*/
create table PEPTIDE
(
   ID_PEPTIDE           int not null auto_increment,
   ID_ANALYSIS          int not null,
   PEP_SEQ              varchar(150),
   PEP_LENGTH           smallint,
   QUERY_NUM            int,
   PEP_RANK             smallint,
   SEARCH_RANK          smallint,
   SCORE                float,
   MISS_CUT             smallint,
   MR_EXP               double,
   MR_CALC              double,
   MR_OBS               double,
   MR_DELTA             double,
   SUBST                text,
   CHARGE               smallint,
   ELUTION_TIME         varchar(50),
   VALID_STATUS         smallint,
   COMMENTS             varchar(250),
   DATA                 text,
   SPEC_COUNT           smallint,
   primary key (ID_PEPTIDE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PEPTIDE_MODIFICATION                                 */
/*==============================================================*/
create table PEPTIDE_MODIFICATION
(
   ID_MODIFICATION      int not null,
   ID_PEPTIDE           int not null,
   MODIF_TYPE           varchar(1) not null,
   POS_STRING           text,
   REF_POS_STRING       text,
   primary key (ID_MODIFICATION, ID_PEPTIDE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PEPTIDE_PROTEIN_ATTRIB                               */
/*==============================================================*/
create table PEPTIDE_PROTEIN_ATTRIB
(
   ID_PROTEIN           int not null,
   ID_PEPTIDE           int not null,
   ID_ANALYSIS          int not null,
   PEP_BEG              int not null,
   PEP_END              int,
   FLANKING_AA          varchar(2),
   IS_SPECIFIC          smallint,
   primary key (ID_PROTEIN, ID_PEPTIDE, PEP_BEG)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROJECT                                              */
/*==============================================================*/
create table PROJECT
(
   ID_PROJECT           int not null,
   ID_IDENTIFIER        int,
   NAME                 varchar(50),
   DES                  varchar(100),
   PROT_VISIBILITY      smallint,
   OWNER                varchar(50),
   WORK_GROUP           varchar(50),
   START_DATE           datetime,
   STATUS               smallint,
   COMMENTS             varchar(250),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   LAST_OPEN            datetime,
   LAST_USER            varchar(20),
   primary key (ID_PROJECT)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROJECT_ACCESS                                       */
/*==============================================================*/
create table PROJECT_ACCESS
(
   ID_USER              varchar(15) not null,
   ID_PROJECT           int not null,
   ID_PROFILE           int not null,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_USER, ID_PROJECT)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROJECT_BIOSAMPLE                                    */
/*==============================================================*/
create table PROJECT_BIOSAMPLE
(
   ID_PROJECT           int not null,
   ID_BIOSAMPLE         int not null,
   primary key (ID_PROJECT, ID_BIOSAMPLE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROJECT_MODIFICATION                                 */
/*==============================================================*/
create table PROJECT_MODIFICATION
(
   ID_PROJECT           int not null,
   ID_MODIFICATION      int not null,
   primary key (ID_PROJECT, ID_MODIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROJECT_PROPERTY                                     */
/*==============================================================*/
create table PROJECT_PROPERTY
(
   ID_PROJECT           int not null,
   ID_PROPERTY          int not null,
   primary key (ID_PROJECT, ID_PROPERTY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROPERTY                                             */
/*==============================================================*/
create table PROPERTY
(
   ID_PROPERTY          int not null auto_increment,
   NAME                 varchar(50),
   DES                  varchar(100),
   PROPERTY_TYPE        varchar(1),
   POSSIBLE_VALUES      varchar(150),
   USE_IN_ANALYSIS      smallint,
   IS_VERIFIED          smallint,
   primary key (ID_PROPERTY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROTEIN                                              */
/*==============================================================*/
create table PROTEIN
(
   ID_PROTEIN           int not null auto_increment,
   ID_PROJECT           int not null,
   ID_MASTER_PROTEIN    int,
   ALIAS                varchar(50),
   IDENTIFIER           varchar(50),
   PROT_DES             varchar(250),
   PROT_LENGTH          int,
   PROT_SEQ             text,
   MW                   float,
   ORGANISM             varchar(100),
   COMMENTS             varchar(250),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_PROTEIN)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROTEIN_QUANTIFICATION                               */
/*==============================================================*/
create table PROTEIN_QUANTIFICATION
(
   ID_PROT_QUANTIF      int not null auto_increment,
   ID_PROTEIN           int not null,
   ID_QUANTIFICATION    int not null,
   ID_QUANTIF_PARAMETER int not null,
   QUANTIF_VALUE        double,
   TARGET_POS           smallint,
   primary key (ID_PROT_QUANTIF)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROTEIN_VALIDATION                                   */
/*==============================================================*/
create table PROTEIN_VALIDATION
(
   ID_PROT_VALID        int not null auto_increment,
   ID_ANALYSIS          int not null,
   IDENTIFIER           varchar(50),
   DB_RANK              smallint,
   MW                   float,
   PROT_DES             varchar(250),
   SEL_STATUS           int,
   NUM_MATCH            int,
   MAX_MATCH            int,
   SCORE                float,
   MAX_SCORE            float,
   ORGANISM             varchar(100),
   CONF_LEVEL           int,
   PROT_LENGTH          int,
   MATCH_GROUP          int,
   primary key (ID_PROT_VALID)
)
engine = InnoDB;

/*==============================================================*/
/* Table : PROTQUANTIF_MODRES                                   */
/*==============================================================*/
create table PROTQUANTIF_MODRES
(
   ID_MODIF_RES         int not null,
   ID_PROT_QUANTIF      int not null,
   primary key (ID_MODIF_RES, ID_PROT_QUANTIF)
)
engine = InnoDB;

/*==============================================================*/
/* Table : QUANTIFICATION                                       */
/*==============================================================*/
create table QUANTIFICATION
(
   ID_QUANTIFICATION    int not null auto_increment,
   ID_MODIFICATION      int,
   ID_QUANTIFICATION_METHOD int not null,
   ID_DESIGN            int,
   NAME                 varchar(50),
   FOCUS                varchar(8),
   QUANTIF_ANNOT        text,
   STATUS               smallint,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_QUANTIFICATION)
)
engine = InnoDB;

/*==============================================================*/
/* Table : QUANTIFICATION_METHOD                                */
/*==============================================================*/
create table QUANTIFICATION_METHOD
(
   ID_QUANTIFICATION_METHOD int not null,
   NAME                 varchar(50),
   CODE                 varchar(15),
   DES                  varchar(100),
   primary key (ID_QUANTIFICATION_METHOD)
)
engine = InnoDB;

/*==============================================================*/
/* Table : QUANTIFICATION_PARAMETER                             */
/*==============================================================*/
create table QUANTIFICATION_PARAMETER
(
   ID_QUANTIF_PARAMETER int not null,
   ID_QUANTIFICATION_METHOD int not null,
   NAME                 varchar(50),
   CODE                 varchar(15),
   UNIT                 varchar(15),
   DES                  varchar(100),
   primary key (ID_QUANTIF_PARAMETER)
)
engine = InnoDB;

/*==============================================================*/
/* Table : QUANTIF_REFRT                                        */
/*==============================================================*/
create table QUANTIF_REFRT
(
   ID_REFERENCE_RT      int not null,
   ID_QUANTIFICATION    int not null,
   SOURCE_RANK          smallint not null,
   NUM_PEP              int,
   SLOPE                float,
   Y_INTERCEPT          float,
   CORRELATION          float,
   DATA                 text,
   primary key (ID_REFERENCE_RT, ID_QUANTIFICATION, SOURCE_RANK)
)
engine = InnoDB;

/*==============================================================*/
/* Table : QUERY_MODIFICATION                                   */
/*==============================================================*/
create table QUERY_MODIFICATION
(
   ID_MODIFICATION      int not null,
   ID_QUERY             int not null,
   PEP_RANK             smallint not null,
   POS_STRING           text,
   REF_POS_STRING       text,
   primary key (ID_MODIFICATION, ID_QUERY, PEP_RANK)
)
engine = InnoDB;

/*==============================================================*/
/* Table : QUERY_VALIDATION                                     */
/*==============================================================*/
create table QUERY_VALIDATION
(
   ID_QUERY             int not null auto_increment,
   ID_ANALYSIS          int not null,
   QUERY_NUM            int,
   EXT_SPECTRUMID       int,
   VALID_STATUS         smallint,
   MASS_DATA            varchar(50),
   MAX_SCORE            float,
   CHARGE               smallint,
   ELUTION_TIME         varchar(50),
   INFO_PEP0            text,
   INFO_PEP1            text,
   INFO_PEP2            text,
   INFO_PEP3            text,
   INFO_PEP4            text,
   INFO_PEP5            text,
   INFO_PEP6            text,
   INFO_PEP7            text,
   INFO_PEP8            text,
   INFO_PEP9            text,
   INFO_PEP10           text,
   primary key (ID_QUERY)
)
engine = InnoDB;

/*==============================================================*/
/* Index : QUERYVAL_QUERYNUM_IDX                                */
/*==============================================================*/
create index QUERYVAL_QUERYNUM_IDX on QUERY_VALIDATION
(
   QUERY_NUM
);

/*==============================================================*/
/* Table : RANK_PROTEIN_MATCH                                   */
/*==============================================================*/
create table RANK_PROTEIN_MATCH
(
   ID_ANALYSIS          int not null,
   QUERY_NUM            int not null,
   PEP_RANK             smallint not null,
   IDENTIFIER           varchar(50) not null,
   MATCH_INFO           varchar(255),
   MATCH_MULTI          int,
   primary key (ID_ANALYSIS, QUERY_NUM, PEP_RANK, IDENTIFIER)
)
engine = InnoDB;

/*==============================================================*/
/* Table : REFERENCE_RT                                         */
/*==============================================================*/
create table REFERENCE_RT
(
   ID_REFERENCE_RT      int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   COLUMN_TYPE          varchar(25),
   PROT_IDENTIFIER      varchar(50),
   PROT_SEQ             text,
   NUM_PEP              int,
   DATA                 text,
   primary key (ID_REFERENCE_RT)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SAMPLE                                               */
/*==============================================================*/
create table SAMPLE
(
   ID_SAMPLE            int not null,
   ID_SPOT              int,
   ID_EXPERIMENT        int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   START_DATE           datetime,
   DISPLAY_POS          smallint,
   COMMENTS             varchar(250),
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_SAMPLE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SPECIES                                              */
/*==============================================================*/
create table SPECIES
(
   ID_SPECIES           int not null auto_increment,
   ID_REF_SPECIES       int,
   COMMON_NAME          varchar(100),
   SCIENTIFIC_NAME      varchar(100),
   TAXONID              int,
   IS_REFERENCE         smallint,
   primary key (ID_SPECIES)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SPECTRUM                                             */
/*==============================================================*/
create table SPECTRUM
(
   ID_SPECTRUM          int not null,
   DATA_FILE            varchar(50),
   PEP_SEQ              varchar(150),
   QUERY_NUM            int,
   SPEC_FILE            varchar(11),
   SCORE                float,
   ELUTION_TIME         varchar(50),
   PEP_RANK             smallint,
   SEARCH_RANK          smallint,
   MISS_CUT             smallint,
   MR_EXP               double,
   MR_CALC              double,
   MR_OBS               double,
   MR_DELTA             double,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   VAR_MOD              text,
   SUBST                text,
   COMMENTS             varchar(250),
   INSTRUMENT           varchar(30),
   primary key (ID_SPECTRUM)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SPECTRUM_MODIFICATION                                */
/*==============================================================*/
create table SPECTRUM_MODIFICATION
(
   ID_MODIFICATION      int not null,
   ID_SPECTRUM          int not null,
   MODIF_TYPE           varchar(1) not null,
   SPECIFICITY          text,
   POS_STRING           text,
   primary key (ID_MODIFICATION, ID_SPECTRUM, MODIF_TYPE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SPOT                                                 */
/*==============================================================*/
create table SPOT
(
   ID_SPOT              int not null,
   ID_GEL2D             int not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   X_POS                smallint,
   Y_POS                smallint,
   ISOELECTRIC_POINT    float,
   MOLECULAR_WEIGHT     float,
   INTENSITY            float,
   EXTERNAL_ID          varchar(25),
   COMMENTS             varchar(250),
   START_DATE           datetime,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_SPOT)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SWATH_LIB                                            */
/*==============================================================*/
create table SWATH_LIB
(
   ID_SWATH_LIB         int not null auto_increment,
   ID_REFERENCE_RT      int,
   NAME                 varchar(50),
   DES                  varchar(100),
   PARAM_STRG           text,
   IDENTIFIER_TYPE      varchar(50),
   ORGANISM             varchar(100),
   INSTRUMENT           varchar(30),
   SPLIT                smallint,
   STATISTICS           text,
   USE_STATUS           varchar(3),
   VERSION_NAME         varchar(20),
   START_DATE           datetime,
   UPDATE_DATE          datetime,
   UPDATE_USER          varchar(20),
   primary key (ID_SWATH_LIB)
)
engine = InnoDB;

/*==============================================================*/
/* Table : SWATH_LIB_MODIFICATION                               */
/*==============================================================*/
create table SWATH_LIB_MODIFICATION
(
   ID_SWATH_LIB         int not null,
   ID_MODIFICATION      int not null,
   SPECIFICITY          text,
   primary key (ID_SWATH_LIB, ID_MODIFICATION)
);

/*==============================================================*/
/* Table : USER_LIST                                            */
/*==============================================================*/
create table USER_LIST
(
   ID_USER              varchar(15) not null,
   USER_STATUS          varchar(20) not null,
   USER_NAME            varchar(50),
   USER_EMAIL           varchar(100),
   USER_TEL             varchar(20),
   LABORATORY           varchar(150),
   USER_INFO            varchar(100),
   USER_PASSWORD        varchar(25),
   LAST_CONNECTION      datetime,
   IS_CLOSED            smallint,
   USER_PREF            varchar(50),
   MASCOT_IDS           varchar(25),
   WORK_GROUP           varchar(50),
   START_DATE           datetime,
   primary key (ID_USER)
)
engine = InnoDB;

/*==============================================================*/
/* Table : USER_PROFILE                                         */
/*==============================================================*/
create table USER_PROFILE
(
   ID_PROFILE           int not null,
   NAME                 varchar(50) not null,
   primary key (ID_PROFILE)
)
engine = InnoDB;

/*==============================================================*/
/* Table : VALIDATION_HISTORY                                   */
/*==============================================================*/
create table VALIDATION_HISTORY
(
   ID_VAL_HISTORY       int not null auto_increment,
   ID_ANALYSIS          int not null,
   STEP                 smallint,
   STATUS               smallint,
   VAL_TYPE             varchar(6),
   PARAM_STRG           text,
   QUERY_VAL_STRG       text,
   PEP_VAL_STRG         text,
   PROT_VAL_STRG        text,
   START_DATE           datetime,
   VALID_USER           varchar(20),
   primary key (ID_VAL_HISTORY)
)
engine = InnoDB;

/*==============================================================*/
/* Table : VALIDATION_TEMPLATE                                  */
/*==============================================================*/
create table VALIDATION_TEMPLATE
(
   ID_VAL_TEMPLATE      int not null auto_increment,
   ID_USER              varchar(15) not null,
   NAME                 varchar(50),
   DES                  varchar(100),
   VAL_TYPE             varchar(6),
   SEARCH_ENGINE        varchar(15),
   MS_TYPE              varchar(10),
   PARAM_STRG           text,
   IS_DEFAULT           smallint,
   primary key (ID_VAL_TEMPLATE)
)
engine = InnoDB;

alter table ANALYSIS add constraint FK_SAMPLE_ANALYSIS foreign key (ID_SAMPLE)
      references SAMPLE (ID_SAMPLE);

alter table ANALYSIS_DATABANK add constraint FK_ANALYSIS_DATABANK foreign key (ID_DATABANK)
      references DATABANK (ID_DATABANK);

alter table ANALYSIS_DATABANK add constraint FK_ANALYSIS_DATABANK2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANALYSIS_MODIFICATION add constraint FK_ANALYSIS_MODIFICATION foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANALYSIS_MODIFICATION add constraint FK_ANALYSIS_MODIFICATION2 foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table ANALYSIS_PROTEIN add constraint FK_ANALYSIS_PROTEIN foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANALYSIS_PROTEIN add constraint FK_ANALYSIS_PROTEIN2 foreign key (ID_PROTEIN)
      references PROTEIN (ID_PROTEIN);

alter table ANALYSIS_REFRT add constraint FK_ANALYSIS_REFRT foreign key (ID_REFERENCE_RT)
      references REFERENCE_RT (ID_REFERENCE_RT);

alter table ANALYSIS_REFRT add constraint FK_ANALYSIS_REFRT2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANALYSIS_SWATH_LIB add constraint FK_ANALYSIS_SWATH_LIB foreign key (ID_SWATH_LIB)
      references SWATH_LIB (ID_SWATH_LIB);

alter table ANALYSIS_SWATH_LIB add constraint FK_ANALYSIS_SWATH_LIB2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANA_COMPARISON add constraint FK_ANA_COMPARISON foreign key (ID_COMPARISON)
      references COMPARISON (ID_COMPARISON);

alter table ANA_COMPARISON add constraint FK_ANA_COMPARISON2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANA_QUANTIFICATION add constraint FK_ANA_QUANTIFICATION foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table ANA_QUANTIFICATION add constraint FK_ANA_QUANTIFICATION2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table ANNOTATIONSET add constraint FK_EXPLORANA_ANNOT foreign key (ID_EXPLORANALYSIS)
      references EXPLORANALYSIS (ID_EXPLORANALYSIS);

alter table BIOSAMPLE add constraint FK_BIOSAMPLESPECIES foreign key (ID_SPECIES)
      references SPECIES (ID_SPECIES);

alter table BIOSAMPLE add constraint FK_REF_BIOSAMPLE foreign key (ID_REFBIOSAMPLE)
      references BIOSAMPLE (ID_BIOSAMPLE);

alter table BIOSAMPLE_PROPERTY add constraint FK_BIOSAMPLE_PROPERTY foreign key (ID_BIOSAMPLE)
      references BIOSAMPLE (ID_BIOSAMPLE);

alter table BIOSAMPLE_PROPERTY add constraint FK_BIOSAMPLE_PROPERTY2 foreign key (ID_PROPERTY)
      references PROPERTY (ID_PROPERTY);

alter table CATEGORY add constraint FK_CLASSIFICATION_CATEGORY foreign key (ID_CLASSIFICATION)
      references CLASSIFICATION (ID_CLASSIFICATION);

alter table CATEGORY_PROTEIN add constraint FK_CAT_CATPROTEIN foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table CATEGORY_PROTEIN add constraint FK_PROT_CATPROTEIN foreign key (ID_PROTEIN)
      references PROTEIN (ID_PROTEIN);

alter table CAT_COMPARISON add constraint FK_CAT_COMPARISON foreign key (ID_COMPARISON)
      references COMPARISON (ID_COMPARISON);

alter table CAT_COMPARISON add constraint FK_CAT_COMPARISON2 foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table CLASSIFICATION add constraint FK_CLASSIFICATION_PROJECT foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table COMPARISON add constraint FK_COMPARISON_FILTER foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table COMPARISON add constraint FK_COMPARISON_PROJECT foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table DATABANK add constraint FK_DATABANK_TYPE_DATABANK foreign key (ID_DBTYPE)
      references DATABANK_TYPE (ID_DBTYPE);

alter table DATABANK_SWATHLIB add constraint FK_DATABANK_SWATHLIB foreign key (ID_DATABANK)
      references DATABANK (ID_DATABANK);

alter table DATABANK_SWATHLIB add constraint FK_DATABANK_SWATHLIB2 foreign key (ID_SWATH_LIB)
      references SWATH_LIB (ID_SWATH_LIB);

alter table DESIGN add constraint FK_EXPDESIGN foreign key (ID_EXPERIMENT)
      references EXPERIMENT (ID_EXPERIMENT);

alter table EXPCONDITION add constraint FK_DESIGNEXPCONDITION foreign key (ID_DESIGN)
      references DESIGN (ID_DESIGN);

alter table EXPCONDITION_QUANTIF add constraint FK_EXPCONDITION_QUANTIF foreign key (ID_EXPCONDITION)
      references EXPCONDITION (ID_EXPCONDITION);

alter table EXPCONDITION_QUANTIF add constraint FK_EXPCONDITION_QUANTIF2 foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table EXPERIMENT add constraint FK_PROJECT_EXPERIMENT foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table EXPERIMENT add constraint FK_SPECIES_EXPERIMENT foreign key (ID_SPECIES)
      references SPECIES (ID_SPECIES);

alter table EXPLORANALYSIS add constraint FK_EXPLORANAFILTER foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table EXPLORANALYSIS add constraint FK_EXP_EXPLORANA foreign key (ID_EXPERIMENT)
      references EXPERIMENT (ID_EXPERIMENT);

alter table EXPLORANA_ANA add constraint FK_EXPLORANA_ANA foreign key (ID_EXPLORANALYSIS)
      references EXPLORANALYSIS (ID_EXPLORANALYSIS);

alter table EXPLORANA_ANA add constraint FK_EXPLORANA_ANA2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table EXPLORANA_QUANTIF add constraint FK_EXPLORANA_QUANTIF foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table EXPLORANA_QUANTIF add constraint FK_EXPLORANA_QUANTIF2 foreign key (ID_EXPLORANALYSIS)
      references EXPLORANALYSIS (ID_EXPLORANALYSIS);

alter table GEL2D add constraint FK_EXPERIMENT_GEL2D foreign key (ID_EXPERIMENT)
      references EXPERIMENT (ID_EXPERIMENT);

alter table GOANA_ANALYSIS add constraint FK_GOANA_ANALYSIS foreign key (ID_GOANALYSIS)
      references GO_ANALYSIS (ID_GOANALYSIS);

alter table GOANA_ANALYSIS add constraint FK_GOANA_ANALYSIS2 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table GOANA_QUANTIFICATION add constraint FK_GOANA_QUANTIFICATION foreign key (ID_GOANALYSIS)
      references GO_ANALYSIS (ID_GOANALYSIS);

alter table GOANA_QUANTIFICATION add constraint FK_GOANA_QUANTIFICATION2 foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table GOANNOTATION add constraint FK_GOANNOT_IDENTIFIER foreign key (ID_IDENTIFIER)
      references IDENTIFIER (ID_IDENTIFIER);

alter table GOANNOTATION add constraint FK_GOANNOT_SPECIES foreign key (ID_SPECIES)
      references SPECIES (ID_SPECIES);

alter table GO_ANALYSIS add constraint FK_EXPGOANALYSIS foreign key (ID_EXPERIMENT)
      references EXPERIMENT (ID_EXPERIMENT);

alter table GO_ANALYSIS add constraint FK_GOANA_CATEGORY foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table GO_ANALYSIS add constraint FK_GOANA_GOANA foreign key (ID_PARENT_GOANA)
      references GO_ANALYSIS (ID_GOANALYSIS);

alter table GO_ANALYSIS add constraint FK_GOANA_GOANNOT foreign key (ID_GOANNOTATION)
      references GOANNOTATION (ID_GOANNOTATION);

alter table GO_ANALYSIS add constraint FK_GOANA_ONTOLOGY foreign key (ID_ONTOLOGY)
      references ONTOLOGY (ID_ONTOLOGY);

alter table IDENTIFIER add constraint FK_IDENT_SPECIES foreign key (ID_SPECIES)
      references SPECIES (ID_SPECIES);

alter table MASTERPROT_IDENTIFIER add constraint FK_MASTERPROT_IDENTIFIER foreign key (ID_IDENTIFIER)
      references IDENTIFIER (ID_IDENTIFIER);

alter table MASTERPROT_IDENTIFIER add constraint FK_MASTERPROT_IDENTIFIER2 foreign key (ID_MASTER_PROTEIN)
      references MASTER_PROTEIN (ID_MASTER_PROTEIN);

alter table MASTER_PROTEIN add constraint FK_MASTERPROT_SPECIES foreign key (ID_SPECIES)
      references SPECIES (ID_SPECIES);

alter table MODIFIED_RESIDUE add constraint FK_QUANTIF_MODRES foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table MODIFICATION_SITE add constraint FK_CATPROTEIN_SITE foreign key (ID_CATEGORY_PROTEIN)
      references CATEGORY_PROTEIN (ID_CATEGORY_PROTEIN);

alter table MODIFICATION_SITE add constraint FK_CAT_MODIFSITE foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table MODIFICATION_SITE add constraint FK_MODIF_MODIFSITE foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table OBSERVATION add constraint FK_ANA_OBS foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table OBSERVATION add constraint FK_BIOSAMPLEOBS foreign key (ID_BIOSAMPLE)
      references BIOSAMPLE (ID_BIOSAMPLE);

alter table OBS_EXPCONDITION add constraint FK_OBS_EXPCONDITION foreign key (ID_EXPCONDITION)
      references EXPCONDITION (ID_EXPCONDITION);

alter table OBS_EXPCONDITION add constraint FK_OBS_EXPCONDITION2 foreign key (ID_OBSERVATION)
      references OBSERVATION (ID_OBSERVATION);

alter table OBS_MODIFICATION add constraint FK_OBS_MODIFICATION foreign key (ID_OBSERVATION)
      references OBSERVATION (ID_OBSERVATION);

alter table OBS_MODIFICATION add constraint FK_OBS_MODIFICATION2 foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table PARENT_EXPLORANALYSIS add constraint FK_PARENT_EXPLORANALYSIS foreign key (ID_EXPLORANALYSIS)
      references EXPLORANALYSIS (ID_EXPLORANALYSIS);

alter table PARENT_EXPLORANALYSIS add constraint FK_PARENT_EXPLORANALYSIS2 foreign key (ID_PARENT_EXPLORANALYSIS)
      references EXPLORANALYSIS (ID_EXPLORANALYSIS);

alter table PARENT_QUANTIFICATION add constraint FK_PARENT_QUANTIFICATION foreign key (ID_PARENT_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table PARENT_QUANTIFICATION add constraint FK_PARENT_QUANTIFICATION2 foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table PARENT_SWATH_LIB add constraint FK_PARENT_SWATH_LIB foreign key (ID_SWATH_LIB)
      references SWATH_LIB (ID_SWATH_LIB);

alter table PARENT_SWATH_LIB add constraint FK_PARENT_SWATH_LIB2 foreign key (ID_PARENT_SWATH_LIB)
      references SWATH_LIB (ID_SWATH_LIB);

alter table PATHWAYANA_ANALYSIS add constraint FK_PATHWAYANA_ANALYSIS foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table PATHWAYANA_ANALYSIS add constraint FK_PATHWAYANA_ANALYSIS2 foreign key (ID_PATHWAY_ANALYSIS)
      references PATHWAY_ANALYSIS (ID_PATHWAY_ANALYSIS);

alter table PATHWAYANA_QUANTIFICATION add constraint FK_PATHWAYANA_QUANTIFICATION foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table PATHWAYANA_QUANTIFICATION add constraint FK_PATHWAYANA_QUANTIFICATION2 foreign key (ID_PATHWAY_ANALYSIS)
      references PATHWAY_ANALYSIS (ID_PATHWAY_ANALYSIS);

alter table PATHWAY_ANALYSIS add constraint FK_EXPPATHWAYANA foreign key (ID_EXPERIMENT)
      references EXPERIMENT (ID_EXPERIMENT);

alter table PATHWAY_ANALYSIS add constraint FK_PATHWAYANA_CATEGORY foreign key (ID_CATEGORY)
      references CATEGORY (ID_CATEGORY);

alter table PATHWAY_ANALYSIS add constraint FK_PATHWAYANA_PATHWAYANA foreign key (ID_PARENT_PATHWAYANA)
      references PATHWAY_ANALYSIS (ID_PATHWAY_ANALYSIS);

alter table PEPTIDE add constraint FK_ANALYSIS_PEPTIDE foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table PEPTIDE_MODIFICATION add constraint FK_PEPTIDE_MODIFICATION foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table PEPTIDE_MODIFICATION add constraint FK_PEPTIDE_MODIFICATION2 foreign key (ID_PEPTIDE)
      references PEPTIDE (ID_PEPTIDE);

alter table PEPTIDE_PROTEIN_ATTRIB add constraint FK_PEPTIDE_PROTEIN_ATTRIB foreign key (ID_PROTEIN)
      references PROTEIN (ID_PROTEIN);

alter table PEPTIDE_PROTEIN_ATTRIB add constraint FK_PEPTIDE_PROTEIN_ATTRIB2 foreign key (ID_PEPTIDE)
      references PEPTIDE (ID_PEPTIDE);

alter table PEPTIDE_PROTEIN_ATTRIB add constraint FK_PEPTIDE_PROTEIN_ATTRIB3 foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table PROJECT add constraint FK_PROJECTIDENTIFIER foreign key (ID_IDENTIFIER)
      references IDENTIFIER (ID_IDENTIFIER);

alter table PROJECT_ACCESS add constraint FK_PROJECT_ACCESS foreign key (ID_USER)
      references USER_LIST (ID_USER);

alter table PROJECT_ACCESS add constraint FK_PROJECT_ACCESS2 foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table PROJECT_ACCESS add constraint FK_PROJECT_ACCESS3 foreign key (ID_PROFILE)
      references USER_PROFILE (ID_PROFILE);

alter table PROJECT_BIOSAMPLE add constraint FK_PROJECT_BIOSAMPLE foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table PROJECT_BIOSAMPLE add constraint FK_PROJECT_BIOSAMPLE2 foreign key (ID_BIOSAMPLE)
      references BIOSAMPLE (ID_BIOSAMPLE);

alter table PROJECT_MODIFICATION add constraint FK_PROJECT_MODIFICATION foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table PROJECT_MODIFICATION add constraint FK_PROJECT_MODIFICATION2 foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table PROJECT_PROPERTY add constraint FK_PROJECT_PROPERTY foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table PROJECT_PROPERTY add constraint FK_PROJECT_PROPERTY2 foreign key (ID_PROPERTY)
      references PROPERTY (ID_PROPERTY);

alter table PROTEIN add constraint FK_PROJECT_PROTEIN foreign key (ID_PROJECT)
      references PROJECT (ID_PROJECT);

alter table PROTEIN add constraint FK_PROT_MASTERPROT foreign key (ID_MASTER_PROTEIN)
      references MASTER_PROTEIN (ID_MASTER_PROTEIN);

alter table PROTEIN_QUANTIFICATION add constraint FK_PROTPROTQUANTIF foreign key (ID_PROTEIN)
      references PROTEIN (ID_PROTEIN);

alter table PROTEIN_QUANTIFICATION add constraint FK_PROTQUANTIFPARAM foreign key (ID_QUANTIF_PARAMETER)
      references QUANTIFICATION_PARAMETER (ID_QUANTIF_PARAMETER);

alter table PROTEIN_QUANTIFICATION add constraint FK_PROTQUANTIFQUANTIF foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table PROTEIN_VALIDATION add constraint FK_ANALYSIS_PROTVAL foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table PROTQUANTIF_MODRES add constraint FK_PROTQUANTIF_MODRES foreign key (ID_MODIF_RES)
      references MODIFIED_RESIDUE (ID_MODIF_RES);

alter table PROTQUANTIF_MODRES add constraint FK_PROTQUANTIF_MODRES2 foreign key (ID_PROT_QUANTIF)
      references PROTEIN_QUANTIFICATION (ID_PROT_QUANTIF);

alter table QUANTIFICATION add constraint FK_ANAQUANTIFMETH foreign key (ID_QUANTIFICATION_METHOD)
      references QUANTIFICATION_METHOD (ID_QUANTIFICATION_METHOD);

alter table QUANTIFICATION add constraint FK_DESIGNQUANTIF foreign key (ID_DESIGN)
      references DESIGN (ID_DESIGN);

alter table QUANTIFICATION add constraint FK_MODIF_QUANTIFICATION foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table QUANTIFICATION_PARAMETER add constraint FK_QUANTIFMETHPARAM foreign key (ID_QUANTIFICATION_METHOD)
      references QUANTIFICATION_METHOD (ID_QUANTIFICATION_METHOD);

alter table QUANTIF_REFRT add constraint FK_QUANTIF_REFRT foreign key (ID_REFERENCE_RT)
      references REFERENCE_RT (ID_REFERENCE_RT);

alter table QUANTIF_REFRT add constraint FK_QUANTIF_REFRT2 foreign key (ID_QUANTIFICATION)
      references QUANTIFICATION (ID_QUANTIFICATION);

alter table QUERY_MODIFICATION add constraint FK_QUERY_MODIFICATION foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table QUERY_MODIFICATION add constraint FK_QUERY_MODIFICATION2 foreign key (ID_QUERY)
      references QUERY_VALIDATION (ID_QUERY);

alter table QUERY_VALIDATION add constraint FK_ANALYSIS_VALIDATION foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table RANK_PROTEIN_MATCH add constraint FK_ANALYSIS_RANK_MATCH foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table SAMPLE add constraint FK_EXPERIMENT_SAMPLE foreign key (ID_EXPERIMENT)
      references EXPERIMENT (ID_EXPERIMENT);

alter table SAMPLE add constraint FK_SPOT_SAMPLE foreign key (ID_SPOT)
      references SPOT (ID_SPOT);

alter table SPECIES add constraint FK_REFSPECIES foreign key (ID_REF_SPECIES)
      references SPECIES (ID_SPECIES);

alter table SPECTRUM_MODIFICATION add constraint FK_SPECTRUM_MODIFICATION foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table SPECTRUM_MODIFICATION add constraint FK_SPECTRUM_MODIFICATION2 foreign key (ID_SPECTRUM)
      references SPECTRUM (ID_SPECTRUM);

alter table SPOT add constraint FK_GEL2D_SPOT foreign key (ID_GEL2D)
      references GEL2D (ID_GEL2D);

alter table SWATH_LIB add constraint FK_SWATHLIB_REFRT foreign key (ID_REFERENCE_RT)
      references REFERENCE_RT (ID_REFERENCE_RT);

alter table SWATH_LIB_MODIFICATION add constraint FK_SWATH_LIB_MODIFICATION foreign key (ID_SWATH_LIB)
      references SWATH_LIB (ID_SWATH_LIB);

alter table SWATH_LIB_MODIFICATION add constraint FK_SWATH_LIB_MODIFICATION2 foreign key (ID_MODIFICATION)
      references MODIFICATION (ID_MODIFICATION);

alter table VALIDATION_HISTORY add constraint FK_ANALYSIS_VALHIS foreign key (ID_ANALYSIS)
      references ANALYSIS (ID_ANALYSIS);

alter table VALIDATION_TEMPLATE add constraint FK_USER_VALTEMPLATE foreign key (ID_USER)
      references USER_LIST (ID_USER);
