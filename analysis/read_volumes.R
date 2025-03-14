#!/usr/bin/env Rscript

# Extract SAMSEG stats from processed LR and HR images
# Diana Giraldo, Sept 2024
# Last update: March 2025

library(dplyr)
library(lubridate)

#File with session info and MRI processing pipeline
SESfile <- "/home/vlab/MS_proj/info_files/session_info_pelt_2024.csv"
# Directory with processed data
PRO_DIR = "/home/vlab/MS_proj/processed_MRI"

###############################################################################
# Read session info
DS <- read.csv(SESfile, header = TRUE, colClasses = "character") %>%
  mutate(Date = as.Date(Date)) %>%
  select(-t0, -Month) %>%
  rename(Subject = Subject.folder,
         Session = Session.folder)

###############################################################################
# Data frame with volume estimations from SAMSEG
DVOL <- data.frame()

for (i in 1:nrow(DS)){
  subcode <- DS$Subject[i]
  sesscode <- DS$Session[i]
  ssfname = ifelse(DS$proc_pipe[i] %in% c("B", "A"), "prettier_samseg", "samseg")
  # Files to read
  inssfile <- sprintf("%s/%s/%s/anat/%s/samseg.stats", PRO_DIR, subcode, sesscode, ssfname)
  intivfile <- sprintf("%s/%s/%s/anat/%s/sbtiv.stats", PRO_DIR, subcode, sesscode, ssfname)
  inunkfile <- sprintf("%s/%s/%s/anat/%s/count_unknownswithinbrain.txt", PRO_DIR, subcode, sesscode, ssfname)
  # Read files
  if (file.exists(inssfile) & file.exists(intivfile) & file.exists(inunkfile)) {
    # Estimated volumes
    segvols <- read.table(inssfile, 
                          header = FALSE, sep = ",",  dec =".", comment.char = "") 
    # Intracraneal volume
    tiv <- read.table(intivfile, 
                      header = FALSE, sep = ",",  dec =".", comment.char = "") 
    # Unknowns in segmentation
    count_unknowns <- read.table(inunkfile)$V1
    tmpvols <- rbind(segvols, tiv) 
    tmpdf <- as.data.frame(t(tmpvols$V2))
    names(tmpdf) <- make.names(paste0("samseg.", gsub("# Measure ", "", tmpvols$V1)))
    tmpdf <- tmpdf %>%
      mutate(Subject = subcode, Session = sesscode) %>%
      select(Subject, Session, everything()) %>%
      mutate(samseg.Unknowns = count_unknowns)
    DVOL <- bind_rows(DVOL, tmpdf)
    rm(segvols, tiv, tmpvols,  tmpdf, inssfile, intivfile, inunkfile, count_unknowns)
  }
}

# Organize SAMSEG estimation of volumes (in mm^3)
DVOL <- DVOL %>%
  select(Subject, Session,
         samseg.Intra.Cranial, samseg.Lesions,
         ends_with("Cerebral.Cortex"), ends_with("Cerebral.White.Matter"),
         ends_with("Cerebellum.Cortex"), ends_with("Cerebellum.White.Matter"),
         ends_with("Amygdala"), ends_with("Hippocampus"),
         ends_with("Accumbens.area"), ends_with("Putamen"), ends_with("Pallidum"),
         ends_with("Caudate"), ends_with("Thalamus"),
         ends_with("choroid.plexus"), ends_with("VentralDC"),
         ends_with("Inf.Lat.Vent"), ends_with("Ventricle"), 
         samseg.CSF, samseg.Brain.Stem, samseg.Unknowns) %>%
  # SAMSEG estimations are in mm^3 -> convert to cm^3 
  mutate(across(samseg.Intra.Cranial:samseg.Brain.Stem, ~ .x/1000))

A <- left_join(DS,DVOL)
rm(DVOL)

###############################################################################
# Data frame with lesion estimations from LST-lpa
thLST <- 0.1
DLS  <- data.frame()

for (i in 1:nrow(DS)){
  subcode <- DS$Subject[i]
  sesscode <- DS$Session[i]
  ssdname = ifelse(DS$proc_pipe[i] %in% c("B", "A"), "prettier_LST", "LST")
  # File to read
  inlstfile <- sprintf("%s/%s/%s/anat/%s/LST_lpa_%0.1f.csv", PRO_DIR, subcode, sesscode, ssdname, thLST)
  # Read file
  if (file.exists(inlstfile)) {
    tmpdf <- read.table(inlstfile,
                        header = TRUE, sep = ",",  dec =".", comment.char = "") %>%
      select(TLV, N) %>%
      rename(lstlpa.Lesion.Volume = TLV, lstlpa.nLesions = N) %>%
      mutate(Subject = subcode, Session = sesscode) %>%
      select(Subject, Session, everything())
    DLS <- bind_rows(DLS, tmpdf)
    rm(tmpdf, inlstfile)
  }
}

A <- left_join(A, DLS)
rm(DLS)

###############################################################################
# Get t0 and Month
A <- A %>%
  rename(MRIdate = Date, MRIpipeline = proc_pipe) %>%
  group_by(Subject) %>%
  mutate(t0 = min(MRIdate),
         Month = round(time_length(interval(t0, MRIdate), "month"))) %>%
  arrange(t0, MRIdate) %>%
  ungroup()

# Calculate extra features
A <- A %>%
  mutate(samseg.Cerebral.GMCortex = samseg.Left.Cerebral.Cortex + samseg.Right.Cerebral.Cortex,
         samseg.Cerebral.WM = samseg.Right.Cerebral.White.Matter + samseg.Left.Cerebral.White.Matter,
         samseg.Ventricles = samseg.Left.Lateral.Ventricle + samseg.Right.Lateral.Ventricle +
           samseg.Right.Inf.Lat.Vent + samseg.Left.Inf.Lat.Vent + samseg.3rd.Ventricle + 
           samseg.4th.Ventricle + samseg.5th.Ventricle) %>%
  select(Subject:Set, t0:Month, lstlpa.nLesions, samseg.Intra.Cranial, lstlpa.Lesion.Volume, everything()) %>%
  # Normalise volumes with estimated intra-cranial volume (TIV)
  mutate(across(lstlpa.Lesion.Volume:samseg.Ventricles,
                ~ (.x/samseg.Intra.Cranial),
                .names = "{col}.normTIV"))
