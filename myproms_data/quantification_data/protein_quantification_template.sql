/*=================================================================*/
/* SQLite code to generate protein_quantification_template.db file */
/* Version 1.0.0                                                   */
/*=================================================================*/

/*==============================================================*/
/* Table : PROTEIN_QUANTIFICATION                               */
/*==============================================================*/
create table PROTEIN_QUANTIFICATION
(
	ID_PROT_QUANTIF      INTEGER PRIMARY KEY,
	ID_PROTEIN           INTEGER NOT NULL,
	SITE_CODE            TEXT,
	ID_QUANTIF_PARAMETER INTEGER NOT NULL,
	TARGET_POS           INTEGER,
	QUANTIF_VALUE        REAL NOT NULL
);
create index PROTPROTQUANTIF_IDX on PROTEIN_QUANTIFICATION (ID_PROTEIN);
create index PROTQUANTIFPARAM_IDX on PROTEIN_QUANTIFICATION (ID_QUANTIF_PARAMETER);
create index PROTQUANTIFTGTPOS_IDX on PROTEIN_QUANTIFICATION (TARGET_POS);


