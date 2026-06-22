library(rsample)

dim(metadata.ER_POS_SURV)


set.seed(123)

metabric_surv_split <- initial_split(
  metadata.ER_POS_SURV,
  prop = 0.8,
  strata = EVENT_STAT)

metadata_surv_train <- training(metabric_surv_split)

metadata_surv_test  <- testing(metabric_surv_split)

train_surv.id <- metadata_surv_train$PATIENT_ID
test_surv.id <- metadata_surv_test$PATIENT_ID

length(train_surv.id)
length(test_surv.id)

table(metadata_surv_train$EVENT_STAT)
table(metadata_surv_test$EVENT_STAT) 

intersect(train_surv.id, test_surv.id)

saveRDS(counts_data, "counts_data.rds")
saveRDS(metadata_surv_train, "metadata_surv_train.rds")
saveRDS(metadata_surv_test, "metadata_surv_test.rds")
saveRDS(train_surv.id, "train_surv_id.rds")
saveRDS(test_surv.id, "test_surv_id.rds")
