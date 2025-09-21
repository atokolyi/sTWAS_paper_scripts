library(bigrquery)

sql = "select concept_id, concept_code from concept where vocabulary_id='ICD10CM'"
tb = bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), sql, billing = Sys.getenv("GOOGLE_PROJECT"))
df = bq_table_download(tb)
saveRDS(df,"icd10cm_concept_map.rds")

sql = "select concept_id, concept_code from concept where vocabulary_id='ICD9CM'"
tb = bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), sql, billing = Sys.getenv("GOOGLE_PROJECT"))
df = bq_table_download(tb)
saveRDS(df,"icd9cm_concept_map.rds")