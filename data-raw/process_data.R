library(tidyverse)

load("data-raw/runDb.RData") # nolint: nonportable_path_linter.
db_ptk <- DB |> ungroup()

load("data-raw/runDb_stk.RData") # nolint: nonportable_path_linter.
db_stk <- DB |> ungroup()

uka_db_full <- rbind(db_ptk, db_stk)
saveRDS(uka_db_full, "data/uka_map_db.RDS") # nolint: nonportable_path_linter.


# kinome mapping file
kinome_mp_file_v1 <- read_delim("data-raw/KinomeMapping/kinase_mapping.txt", # nolint: nonportable_path_linter.
  delim = "\t",
  na = c(
    "N/A", "Not found by Ethan", "Not found by Jake", "Not found",
    "Not Found By Ethan", "Not Found By Jake"
  )
)

kinome_mp_file_v2 <- read_delim(
  "data-raw/KinomeMapping/2021_05_20-creedenzymatic_map.txt", # nolint: nonportable_path_linter.
  delim = "\t",
  na = c(
    "N/A", "Not found by Ethan", "Not found by Jake", "Not found",
    "Not Found By Ethan", "Not Found By Jake"
  )
) |>
  filter(!is.na(hgnc_symbol)) |>
  mutate_at(
    c(
      "class", "group", "family", "subfamily",
      "krsa_id", "uka_id", "kea3_id", "ptmsea_id"
    ), toupper
  ) |>
  select(1L:12L)

kinome_mp_file_v3 <- read_delim(
  "data-raw/KinomeMapping/2021_06_11-creedenzymatic_map.txt", # nolint: nonportable_path_linter.
  delim = "\t"
) |>
  mutate_at(c(
    "group", "family", "subfamily",
    "krsa_id", "uka_id", "kea3_id", "ptmsea_id"
  ), toupper) |>
  filter(!is.na(hgnc_symbol)) |>
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) |>
  select(1L:26L)


kinome_mp_file_v4 <- read_delim(
  "data-raw/KinomeMapping/2021_07_09-creedenzymatic_map.txt", # nolint: nonportable_path_linter.
  delim = "\t",
  na = "N/A"
) |>
  mutate_at(c(
    "group", "family", "subfamily", "krsa_id",
    "uka_id", "kea3_id", "ptmsea_id"
  ), toupper) |>
  filter(!is.na(hgnc_symbol)) |>
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) |>
  select(1L:26L)

kinome_mp_file_v5 <- read_delim(
  "data-raw/KinomeMapping/2021_07_09-creedenzymatic_map_KA_edits.txt", # nolint: nonportable_path_linter.
  delim = "\t",
  na = "N/A"
) |>
  mutate_at(c(
    "group", "family", "subfamily",
    "krsa_id", "uka_id", "kea3_id", "ptmsea_id"
  ), toupper) |>
  filter(!is.na(hgnc_symbol)) |>
  mutate(subfamily = ifelse(subfamily == "N/A", family, subfamily)) |>
  mutate(
    kea3_id = ifelse(kea3_id == "NOT FOUND" |
      is.na(kea3_id), hgnc_symbol, kea3_id) # nolint: indentation_linter.
  )

kinome_mp_file_v6 <- read_delim(
  "data-raw/KinomeMapping/2024_03_07-creedenzymatic_map_AH-ASI_edits.txt", # nolint: nonportable_path_linter.
  delim = "\t",
  na = "N/A"
) |>
  mutate_at(c(
    "group", "family", "subfamily", "krsa_id",
    "uka_id", "kea3_id", "ptmsea_id"
  ), toupper) |>
  filter(!is.na(hgnc_symbol)) |>
  mutate(subfamily = ifelse(is.na(subfamily), family, subfamily)) |>
  mutate(kea3_id = ifelse(kea3_id == "NOT FOUND" | is.na(kea3_id), hgnc_symbol, kea3_id)) |>
  filter(!is.na(hgnc_id))


kinome_mp_file <- kinome_mp_file_v6


# peptide to HGNC

stk_pamchip_87102_mapping <- read_delim(
  "data-raw/Pamchips_Layout/2021_10_13-JFC_complete-stk_peptides_mapping-KA_edit.txt", # nolint: nonportable_path_linter, line_length_linter.
  delim = "\t"
) |>
  select(ID, HGNC)

ptk_pamchip_86402_mapping <- read_delim(
  "data-raw/Pamchips_Layout/2021_10_13-JFC_complete-ptk_peptides_mapping.txt", # nolint: nonportable_path_linter, line_length_linter.
  delim = "\t"
) |>
  select(ID, HGNC)

add_stk <- tibble::tibble(
  ID = c("H2B1B_ 27_40", "E1A_ADE05_212_224"),
  HGNC = c("H2BC3", NA)
)

add_ptk <- tibble::tibble(
  ID = "CD3Z_147_159",
  HGNC = "CD247"
)

stk_pamchip_87102_mapping <- rbind(stk_pamchip_87102_mapping, add_stk)
ptk_pamchip_86402_mapping <- rbind(ptk_pamchip_86402_mapping, add_ptk)


stk_pamchip_87102_array_layout <- readxl::read_xlsx(
  "data-raw/Pamchips_Layout/87102_ArrayLayout.xlsx" # nolint: nonportable_path_linter.
) |>
  filter(ID != "#REF", Ser != "NA", Thr != "NA") |>
  select(ID, Ser, Thr, UniprotAccession) |>
  mutate(
    Ser = gsub("[", "", Ser, fixed = TRUE),
    Ser = gsub("]", "", Ser, fixed = TRUE),
    Ser = gsub("\\s+", "", Ser),
    Thr = gsub("[", "", Thr, fixed = TRUE),
    Thr = gsub("]", "", Thr, fixed = TRUE),
    Thr = gsub("\\s+", "", Thr),
    UniprotAccession = gsub("†", "", UniprotAccession, fixed = TRUE),
    ID = gsub("\\s+", "", ID)
  )

ptk_pamchip_86402_array_layout <- readxl::read_xlsx(
  "data-raw/Pamchips_Layout/86402_ArrayLayout.xlsx" # nolint: nonportable_path_linter.
) |>
  filter(ID != "#REF", Tyr != "NA") |>
  select(ID, Tyr, UniprotAccession) |>
  mutate(
    Tyr = gsub("[", "", Tyr, fixed = TRUE),
    Tyr = gsub("]", "", Tyr, fixed = TRUE),
    Tyr = gsub("\\s+", "", Tyr),
    UniprotAccession = gsub("†", "", UniprotAccession, fixed = TRUE),
    ID = gsub("\\s+", "", ID)
  )


iptm_format_ptk <- ptk_pamchip_86402_array_layout |>
  separate_rows(Tyr, sep = ",") |>
  mutate(site_residue = "Y") |>
  select(
    substrate_ac = UniprotAccession,
    site_residue,
    site_position = Tyr
  ) |>
  filter(site_position != "")

utils::write.table(iptm_format_ptk |> filter(site_position != ""),
  "data-raw/iptm_net_input_ptk_chip.csv", # nolint: nonportable_path_linter.
  col.names = FALSE, sep = "\t", quote = FALSE, row.names = FALSE
)

# iptm_format_list <- split(iptm_format, seq(nrow(iptm_format))) |> map(as.list)
# iptmnetr::get_ptm_enzymes_from_list(sites) -> iptm_res

iptm_res_ptk <- iptmnetr::get_ptm_enzymes_from_file(
  "data-raw/iptm_net_input_ptk_chip.csv" # nolint: nonportable_path_linter.
)

iptm_res_ptk <- iptm_res_ptk |>
  filter(!enz_name %in% c(
    "conventional protein kinase C", "IkappaB kinase complex (human)",
    "AKT kinase", "aurora kinase", "IKBKG", "Abl fusion"
  )) |>
  mutate(enz_name = case_when(
    enz_name == "FAK" ~ "PTK2",
    TRUE ~ enz_name
  ))

iptmnet_hgnc_missing <- setdiff(
  {
    iptm_res_ptk[["enz_name"]] |> unique()
  },
  kinome_mp_file[["hgnc_symbol"]]
)

iptm_res_ptk <- iptm_res_ptk |>
  filter(!enz_name %in% iptmnet_hgnc_missing)


iptm_res_ptk <- iptm_res_ptk |>
  select(head = enz_name, sub_id, site) |>
  mutate(entry = paste0(sub_id, ";", site, "-p;u")) |>
  select(head, entry) |>
  group_by(head) |>
  mutate(entry = list(unique(entry)), len = n()) |>
  ungroup() |>
  distinct()


iptm_res_ptk <- setNames(as.list(as.data.frame(t(iptm_res_ptk))), iptm_res_ptk[["head"]])

cmapR::write_gmt(
  iptm_res_ptk,
  "data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_ptk.gmt" # nolint: nonportable_path_linter.
)

ptm_sea_iptmnet_mapping_ptk <- readLines(
  "data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_ptk.gmt" # nolint: nonportable_path_linter.
)


## STK iptment
stk_pamchip_87102_array_layout_ptmsea <- readxl::read_xlsx(
  "data-raw/Pamchips_Layout/87102_ArrayLayout.xlsx" # nolint: nonportable_path_linter.
) |>
  filter(ID != "#REF") |>
  select(Ser, Thr, substrate_ac = UniprotAccession) |>
  mutate(
    Ser = ifelse(Ser == "[]", NA, Ser),
    Thr = ifelse(Thr == "[]", NA, Thr),
    Ser = str_replace(Ser, fixed("["), ""),
    Ser = str_replace(Ser, fixed("]"), ""),
    Thr = str_replace(Thr, fixed("["), ""),
    Thr = str_replace(Thr, fixed("]"), "")
  )

iptm_format_stk <- stk_pamchip_87102_array_layout_ptmsea |>
  separate_rows(Ser, sep = ", ") |>
  separate_rows(Thr, sep = ", ") |>
  rename(S = Ser, T = Thr) |>
  pivot_longer(
    cols = 1L:2L,
    names_to = "site_residue",
    values_to = "site_position"
  ) |>
  filter(!is.na(site_position), substrate_ac != "NA")


utils::write.table(iptm_format_stk,
  "data-raw/iptm_net_input_stk_chip.csv", # nolint: nonportable_path_linter.
  col.names = FALSE, sep = "\t",
  quote = FALSE, row.names = FALSE
)


iptm_res_stk <- iptmnetr::get_ptm_enzymes_from_file(
  "data-raw/iptm_net_input_stk_chip.csv" # nolint: nonportable_path_linter.
)

iptm_res_stk <- iptm_res_stk |>
  filter(!enz_name %in% c(
    "conventional protein kinase C", "IkappaB kinase complex (human)",
    "AKT kinase", "aurora kinase", "IKBKG"
  )) |>
  mutate(enz_name = case_when(
    enz_name == "hMAP3K7/Phos:1" ~ "MAP3K7",
    enz_name == "hAURKA" ~ "AURKA",
    enz_name == "hNUAK1" ~ "NUAK1",
    enz_name == "hCDK9" ~ "CDK9",
    enz_name == "hATM" ~ "ATM",
    enz_name == "hIKBKE" ~ "IKBKE",
    enz_name == "hAURKB" ~ "AURKB",
    enz_name == "Pdk1" ~ "PDK1",
    enz_name == "CHK2" ~ "CHEK2",
    TRUE ~ enz_name
  ))

iptmnet_hgnc_missing <- setdiff(
  {
    iptm_res_stk[["enz_name"]] |> unique()
  },
  kinome_mp_file[["hgnc_symbol"]]
)
iptm_res_stk <- iptm_res_stk |>
  filter(!enz_name %in% iptmnet_hgnc_missing)


iptm_res_stk <- iptm_res_stk |>
  select(head = enz_name, sub_id, site) |>
  mutate(entry = paste0(sub_id, ";", site, "-p;u")) |>
  select(head, entry) |>
  group_by(head) |>
  mutate(entry = list(unique(entry)), len = n()) |>
  ungroup() |>
  distinct()

iptm_res_stk <- setNames(as.list(as.data.frame(t(iptm_res_stk))), iptm_res_stk[["head"]])

cmapR::write_gmt(
  iptm_res_stk,
  "data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_stk.gmt" # nolint: nonportable_path_linter.
)

ptm_sea_iptmnet_mapping_stk <- readLines(
  "data-raw/PTM_SEA_DB/ptm_sea_iptmnet_mapping_stk.gmt" # nolint: nonportable_path_linter.
)



ptk_pamchip_86402_array_layout_ptmsea <- ptk_pamchip_86402_array_layout |>
  separate_rows(Tyr, sep = ",") |>
  mutate(
    Position = paste0("Y", Tyr),
    PTM_SEA_ID = paste0(UniprotAccession, ";", Position, "-p")
  ) |>
  dplyr::select(ID, PTM_SEA_ID)


stk_pamchip_87102_array_layout_ptmsea <- stk_pamchip_87102_array_layout |>
  separate_rows(Ser, sep = ",") |>
  separate_rows(Thr, sep = ",") |>
  mutate(
    Ser = paste0("S", Ser),
    Thr = paste0("T", Thr)
  ) |>
  pivot_longer(2L:3L, names_to = "AA", values_to = "Position") |>
  filter(nchar(Position) != 1L) |>
  select(-AA) |>
  mutate(PTM_SEA_ID = paste0(UniprotAccession, ";", Position, "-p")) |>
  dplyr::select(ID, PTM_SEA_ID)

# PTM DB
ptmsea_all <- cmapR::parse_gmt(
  "data-raw/PTM_SEA_DB/ptm.sig.db.all.uniprot.human.v1.9.0.gmt" # nolint: nonportable_path_linter.
) # nolint: nonportable_path_linter.

covertList <- function(x) {
  tibble(
    entry = x[["entry"]],
    head = x[["head"]]
  )
}

processedDB <- map_df(ptmsea_all, covertList)

filterDBS_ptk <- processedDB |> # nolint: object_name_linter.
  filter(grepl("Kinase", head, ignore.case = TRUE)) |>
  filter(grepl(";Y", entry, fixed = TRUE)) |>
  group_by(head) |>
  mutate(len = n()) |>
  ungroup()


onlyChipPeps_dbs_ptk <- filterDBS_ptk |> # nolint: object_name_linter.
  mutate(
    ids = entry,
    ids = str_remove(ids, fixed(";u")),
    ids = str_remove(ids, fixed("-2"))
  ) |>
  filter(ids %in% ptk_pamchip_86402_array_layout_ptmsea[["PTM_SEA_ID"]]) |>
  select(head, entry) |>
  group_by(head) |>
  mutate(len = n()) |>
  ungroup() |>
  group_by(head) |>
  mutate(entry = list(unique(entry))) |>
  ungroup() |>
  distinct() |>
  select(head, entry, len)

onlyChipPeps_dbs_ptk <- setNames( # nolint: object_name_linter.
  as.list(
    as.data.frame(
      t(onlyChipPeps_dbs_ptk)
    )
  ), onlyChipPeps_dbs_ptk[["head"]]
)

cmapR::write_gmt(
  onlyChipPeps_dbs_ptk,
  "data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt" # nolint: nonportable_path_linter.
)

ptk_pamchip_86402_onlyChipPeps_dbs <- readLines( # nolint: object_name_linter.
  "data-raw/PTM_SEA_DB/ptk_pamchip_86402_onlyChipPeps_dbs.gmt" # nolint: nonportable_path_linter.
)

filterDBS_stk <- processedDB |> # nolint: object_name_linter.
  filter(grepl("Kinase", head, ignore.case = TRUE)) |>
  filter(grepl(";T|S", entry)) |>
  group_by(head) |>
  mutate(len = n()) |>
  ungroup()


onlyChipPeps_dbs_stk <- filterDBS_stk |> # nolint: object_name_linter.
  mutate(ids = entry, ids = str_remove(ids, fixed(";u"))) |>
  filter(ids %in% stk_pamchip_87102_array_layout_ptmsea[["PTM_SEA_ID"]]) |>
  select(head, entry) |>
  group_by(head) |>
  mutate(len = n()) |>
  ungroup() |>
  group_by(head) |>
  mutate(entry = list(unique(entry))) |>
  ungroup() |>
  distinct() |>
  select(head, entry, len)

onlyChipPeps_dbs_stk <- setNames( # nolint: object_name_linter.
  as.list(
    as.data.frame(
      t(onlyChipPeps_dbs_stk)
    )
  ), onlyChipPeps_dbs_stk[["head"]]
)

cmapR::write_gmt(
  onlyChipPeps_dbs_stk,
  "data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt" # nolint: object_name_linter, nonportable_path_linter.
)

stk_pamchip_87102_onlyChipPeps_dbs <- readLines( # nolint: object_name_linter.
  "data-raw/PTM_SEA_DB/stk_pamchip_87102_onlyChipPeps_dbs.gmt" # nolint: nonportable_path_linter.
)


ptmsea_all_dbs <- readLines(
  "data-raw/PTM_SEA_DB/ptm.sig.db.all.uniprot.human.v1.9.0.gmt" # nolint: nonportable_path_linter.
)

stk_peps <- KRSA::KRSA_coverage_STK_PamChip_87102_v2[["Substrates"]] |>
  unique() |>
  as.character()
ptk_peps <- KRSA::KRSA_coverage_PTK_PamChip_86402_v1[["Substrates"]] |>
  unique() |>
  as.character()

setdiff(stk_peps, stk_pamchip_87102_mapping[["ID"]])



usethis::use_data(uka_db_full,
  kinome_mp_file,
  kinome_mp_file_v1,
  kinome_mp_file_v2,
  kinome_mp_file_v3,
  kinome_mp_file_v4,
  kinome_mp_file_v5,
  kinome_mp_file_v6,
  stk_pamchip_87102_mapping,
  stk_pamchip_87102_array_layout_ptmsea,
  ptk_pamchip_86402_mapping,
  ptk_pamchip_86402_array_layout_ptmsea,
  ptk_pamchip_86402_onlyChipPeps_dbs,
  stk_pamchip_87102_onlyChipPeps_dbs,
  ptm_sea_iptmnet_mapping_ptk,
  ptm_sea_iptmnet_mapping_stk,
  ptmsea_all_dbs,
  overwrite = TRUE
)
