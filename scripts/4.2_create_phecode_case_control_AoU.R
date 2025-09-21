options(stringsAsFactors=F)
library(tidyverse)
library(bigrquery)
library('stringr')

# Load ICD-CM/PheCode mappings
icdcm_ur = read.table("ICD-CM_to_phecode_unrolled.txt",sep="\t",header=T,colClasses=c(rep("character",3)))
icdcm9 = icdcm_ur[which(icdcm_ur$flag=="9"),]
icdcm10 = icdcm_ur[which(icdcm_ur$flag=="10"),]
cmap_9 = readRDS("icd9cm_concept_map.rds")
cmap_10 = readRDS("icd10cm_concept_map.rds")

# Expand PheCode ranges for exclusion group
# From function.expandPhecodes.r (https://github.com/umich-cphds/createUKBphenome)
expandPhecodes <- function(x,addIntegers=T){
    if(is.na(x) | x == "") return(NA)
    if(grepl("\\-",x)){
        # split range
        # character prefix
        i1 <- strsplit(x,"-")[[1]]

        # numeric length of digits before "."
        nprefix <- max(nchar(gsub("\\..+","",i1)))
        # numbers of digits
        ndigits <- max(c(nchar(gsub("^[0-9]+[\\.]{0,1}","",i1)),0))
        # add "." to length of formatted number if present
        addDot <- max(as.numeric(grepl("\\.",i1)))
        # create sequence of ICD codes
        seq1 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^ndigits))
        # format sequence to match intput
        seq1 <- formatC(seq1, format='f', digits=ndigits,width=nprefix+ndigits+addDot,flag=0)
        # add integers if within range
        if(addIntegers) seq1 <- unique(sort(c(seq1,gsub("\\..+","",seq1[which(round(as.numeric(seq1)) == as.numeric(seq1))]))))

        if(ndigits == 2){
            seq2 <- seq(as.numeric(i1[1]),as.numeric(i1[2]),(1/10^(ndigits-1)))
            seq2 <- formatC(seq2, format='f', digits=ndigits-1,width=nprefix+ndigits+addDot-1,flag=0)
            seq1 <- unique(sort(c(seq1,seq2)))
        }
        return(seq1)
    } else {
        return(x)
    }
}

# Load UKBB PheCodes with ≥ 200 cases
phens = read.table("UKB_PHENOME_DESCRIPTION_20240807_MIN200.txt",
                   header=TRUE,colClasses=c("character","character",rep(NA,8)),sep="\t")
rownames(phens) = paste0("X",phens$phecode)

# Parallelize per PheCode
args = commandArgs(trailingOnly=TRUE)
N = as.integer(args[1])

# Get the ICD9-CM and ICD10-CM codes for the PheCode
include_cm9 = unique(icdcm9[which(icdcm9$phecode==phens$phecode[N]),]$ICD)
include_cm10 = unique(icdcm10[which(icdcm10$phecode==phens$phecode[N]),]$ICD)

# Get the ICD9-CM and ICD10-CM codes for the exclusion criteria of the PheCode
ranges = unlist(str_split(phens$phecode_exclude_range[N],", "))
exclude_phecode = c()
for (r in ranges) {
    exclude_phecode = c(exclude_phecode, expandPhecodes(r))
}
exclude_phecode = unique(exclude_phecode)
exclude_phecode9 = exclude_phecode[which(exclude_phecode %in% unique(icdcm9$phecode))]
exclude_phecode10 = exclude_phecode[which(exclude_phecode %in% unique(icdcm10$phecode))]
exclude_cm9 = unique(icdcm9[which(icdcm9$phecode %in% exclude_phecode9),]$ICD)
exclude_cm10 = unique(icdcm10[which(icdcm10$phecode %in% exclude_phecode10),]$ICD)

# Map the include and exclude codes to concepts (used by AoU)
include_concepts = unique(c(cmap_9[which(cmap_9$concept_code %in% include_cm9),]$concept_id,
                           cmap_10[which(cmap_10$concept_code %in% include_cm10),]$concept_id))
exclude_concepts = unique(c(cmap_9[which(cmap_9$concept_code %in% exclude_cm9),]$concept_id,
                           cmap_10[which(cmap_10$concept_code %in% exclude_cm10),]$concept_id))

# Query AoU EHR database for all matches to the include or exclude concepts
both_concepts = unique(c(include_concepts,exclude_concepts))
both_concepts_str = paste0(both_concepts,collapse=", ")
get_sqtl = function(concept) {
    return(paste0("
    SELECT
        c_occurrence.person_id,
        c_occurrence.condition_concept_id,
        c_standard_concept.concept_name as standard_concept_name,
        c_standard_concept.concept_code as standard_concept_code,
        c_standard_concept.vocabulary_id as standard_vocabulary,
        c_occurrence.condition_start_datetime,
        c_occurrence.condition_end_datetime,
        c_occurrence.condition_type_concept_id,
        c_type.concept_name as condition_type_concept_name,
        c_occurrence.stop_reason,
        c_occurrence.visit_occurrence_id,
        visit.concept_name as visit_occurrence_concept_name,
        c_occurrence.condition_source_value,
        c_occurrence.condition_source_concept_id,
        c_source_concept.concept_name as source_concept_name,
        c_source_concept.concept_code as source_concept_code,
        c_source_concept.vocabulary_id as source_vocabulary,
        c_occurrence.condition_status_source_value,
        c_occurrence.condition_status_concept_id,
        c_status.concept_name as condition_status_concept_name
    FROM
        ( SELECT
            *
        FROM
            `condition_occurrence` c_occurrence
        WHERE
            (
                condition_source_concept_id IN (SELECT
                    DISTINCT c.concept_id
                FROM
                    `cb_criteria` c
                JOIN
                    (SELECT
                        CAST(cr.id as string) AS id
                    FROM
                        `cb_criteria` cr
                    WHERE
                        concept_id IN (",concept,")
                        AND full_text LIKE '%_rank1]%'      ) a
                        ON (c.path LIKE CONCAT('%.', a.id, '.%')
                        OR c.path LIKE CONCAT('%.', a.id)
                        OR c.path LIKE CONCAT(a.id, '.%')
                        OR c.path = a.id)
                WHERE
                    is_standard = 0
                    AND is_selectable = 1)
            )
            AND (
                c_occurrence.PERSON_ID IN (SELECT
                    distinct person_id
                FROM
                    `cb_search_person` cb_search_person
                WHERE
cb_search_person.person_id IN (SELECT
                        person_id
                    FROM
                        `cb_search_person` p
                    WHERE
                        has_ehr_data = 1 )
                    AND cb_search_person.person_id IN (SELECT
                        person_id
                    FROM
                        `cb_search_person` p
                    WHERE
                        has_array_data = 1 ) )
            )) c_occurrence
    LEFT JOIN
        `concept` c_standard_concept
            ON c_occurrence.condition_concept_id = c_standard_concept.concept_id
    LEFT JOIN
        `concept` c_type
            ON c_occurrence.condition_type_concept_id = c_type.concept_id
    LEFT JOIN
        `visit_occurrence` v
            ON c_occurrence.visit_occurrence_id = v.visit_occurrence_id
    LEFT JOIN
        `concept` visit
            ON v.visit_concept_id = visit.concept_id
    LEFT JOIN
        `concept` c_source_concept
            ON c_occurrence.condition_source_concept_id = c_source_concept.concept_id
    LEFT JOIN
        `concept` c_status
            ON c_occurrence.condition_status_concept_id = c_status.concept_id"))
}
sql = get_sqtl(both_concepts_str)
bq_auth()
tb = bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), sql, billing = Sys.getenv("GOOGLE_PROJECT"))
df = bq_table_download(tb,page_size = 10000)
df$person_id = as.character(df$person_id)

# Get the person IDs of the people with an include or an exclude concept match
include_people = unique(df[which(df$condition_source_concept_id %in% include_concepts),]$person_id)
exclude_people = unique(df[which(df$condition_source_concept_id %in% exclude_concepts),]$person_id)
exclude_people = exclude_people[which(!exclude_people %in% include_people)]

df_include = df[which(df$person_id %in% include_people & df$condition_source_concept_id %in% include_concepts),]
df_include = df_include[order(df_include$condition_start_datetime),]
df_include = df_include[!duplicated(df_include$person_id),]
df_include = as.data.frame(df_include)
rownames(df_include) = as.character(df_include$person_id)

# Load cohort covariates
demo = read.table("demographics_201958.tsv",header=T,sep="\t")
demo$person_id = as.character(demo$person_id)

# Limit sex-specific PheCodes
if (phens$sex[N]=="Male") {
    exclude_people = unique(c(exclude_people,demo[which(demo$sex_at_birth=="Female"),]$person_id))
}
if (phens$sex[N]=="Female") {
    exclude_people = unique(c(exclude_people,demo[which(demo$sex_at_birth=="Male"),]$person_id))
}

# Compile PheCode case/control matrix
demo$has_vp = demo$person_id %in% include_people
if (length(exclude_people)>0) {
        demo[which(demo$person_id %in% exclude_people),]$has_vp = NA
}
demo$first_date_phecode = NA
demo[which(demo$has_vp==1),]$first_date_phecode = as.Date(df_include[demo[which(demo$has_vp==1),]$person_id,]$condition_start_datetime)
demo$age_first_phecode = NA
demo[which(demo$has_vp==1),]$age_first_phecode = time_length(difftime(
   as.Date(demo[which(demo$has_vp==1),]$first_date_phecode),
    as.Date(demo[which(demo$has_vp==1),]$date_of_birth)
    ),unit="years")
demo$age_at_freeze_inc_death_inc_phecode = demo$age_at_freeze_inc_death
demo[which(demo$has_vp==1),]$age_at_freeze_inc_death_inc_phecode = demo[which(demo$has_vp==1),]$age_first_phecode
demo$sex_at_birth = as.integer(demo$sex_at_birth=="Male")
demo$resid = -1

write.table(demo,paste0("demo_phecodes_out/",N,".tsv"),col.names=T,row.names=F,sep="\t",quote=F)