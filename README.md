# WNV_QTL_Resource
## Description

 This repostitory contains the workflow associated with the manuscript titled'Adaptive immune response to West Nile virus infection in the Collaborative Cross mouse model: A database of cellular phenotypes and Quantitative Trait Loci'. There are four Jupyter notebooks with step by step instructions, along with R programs that contain the code in these notebooks. These processes were used to generate all of the analyese presented in this manuscript.

## Getting started

### Prerequisites/Dependencies

 

### Inputs and Outputs

The Redcap file extract, data dictionary and OMOP concept mapping files for Redcap questions and answers need to be uploaded to OMOP_ETL. These will be used as inputs for the ETL. This process will directly update Postgres tables in the OMOP database on OMOP_DB.

The following six Athena OMOP vocabulary file extracts also need to be uploaded to OMOP_ETL: DOMAIN, RELATIONSHIP, VOCABULARY, CONCEPT, CONCEPT_SYNONYM, CONCEPT_CLASS. 

The DQD process on OMOP_DQD will output a JSON file to be used for reporting.

### Installing: Steps to run all proccess:

1.	Create a working OMOP database (OMOP) and schema (aireadi_omop) in Postgres. Create the tables in this database/schema with the following SQL table creation statements. This is done on the OMOP_DB VM:
- OMOPCDM_postgresql_5.4_aireadi_ddl.sql

2.	Download vocabulary tables from ATHENA, if necessary, and transfer the csv files to OMOP_DB(concept_ancestor, concept_relationship) and OMOP_ETL (concept_class, concept_synonym, concept, domain, relationship, vocabulary). These files will be tab delimited. The following is the ATHENA website:
- https://athena.ohdsi.org/search-terms/start 

4.	Load 6 of the 8 vocabulary tables into Postgres with the following Python scripts from the OMOP_ETL VM:
-	OMOP_Datamart_Step_2_Azure_Concept_Class.py
-	OMOP_Datamart_Step_2_Azure_Concept_Synonym.py
-	OMOP_Datamart_Step_2_Azure_Concept.py
-	OMOP_Datamart_Step_2_Azure_Domain.py
-	OMOP_Datamart_Step_2_Azure_Relationship.py
-	OMOP_Datamart_Step_2_Azure_Vocabulary.py
