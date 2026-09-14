# Mapping from the data's tissue groups to the anatomical organs gganatogram
# knows how to draw. One tissue group fans out to several organs.

# Tissue groups whose organ set is the same for both sexes.
BASE_ORGAN_MAP <- list(
    "blood" = "leukocyte",
    "large intestine" = c("caecum", "rectum", "colon"),
    "respiratory system" = c("lung", "bronchus", "pulmonary_valve", "diaphragm"),
    "spleen" = "spleen",
    # The 2026-09 data splits what used to be one "gastrointestinal accessory
    # organs" group into three tissue groups of its own.
    "liver" = "liver",
    "pancreas" = "pancreas",
    "gallbladder" = "gall_bladder",
    "lymphatic system" = c("lymph_node", "tonsil"),
    "bone marrow" = "bone_marrow",
    "small intestine" = c("small_intestine", "ileum", "duodenum"),
    "renal system" = c("kidney", "renal_cortex", "urinary_bladder"),
    "integumentary system (skin)" = "skin",
    "cerebral lobes and cortical areas" = c(
        "brain", "frontal_cortex", "temporal_lobe", "prefrontal_cortex"
    ),
    "endocrine system" = c("thyroid_gland", "pituitary_gland", "adrenal_gland"),
    "trachea" = "trachea",
    "stomach" = "stomach",
    "oesophagus" = c("esophagus", "gastroesophageal_junction"),
    "cardiovascular system" = c(
        "heart", "aorta", "left_ventricle", "mitral_valve",
        "tricuspid_valve", "atrial_appendage", "coronary_artery"
    ),
    "adipose tissue" = "adipose_tissue",
    "brainstem and cerebellar structures" = c(
        "spinal_cord", "cerebellum", "cerebellar_hemisphere"
    )
)

# "tongue" and "amygdala" are male-only because gganatogram has no female
# version of either: setdiff(hgMale_key$organ, hgFemale_key$organ) is tongue,
# penis, vas_deferens, seminal_vesicle, epididymis, amygdala, testis, prostate.
# The asymmetry is the package's, not an oversight here.
MALE_ORGAN_MAP <- c(BASE_ORGAN_MAP, list(
    "nasal, oral, and pharyngeal regions" = c(
        "nose", "nasal_pharynx", "throat", "tongue",
        "salivary_gland", "parotid_gland", "submandibular_gland"
    ),
    "sensory-related structures" = c("nerve", "amygdala"),
    "prostate" = "prostate"
))

FEMALE_ORGAN_MAP <- c(BASE_ORGAN_MAP, list(
    "nasal, oral, and pharyngeal regions" = c(
        "nose", "nasal_pharynx", "nasal_septum", "throat",
        "salivary_gland", "parotid_gland", "submandibular_gland"
    ),
    "sensory-related structures" = "nerve",
    "breast" = "breast",
    "female reproductive system" = c(
        "ectocervix", "uterine_cervix", "vagina", "ovary",
        "uterus", "fallopian_tube", "endometrium"
    )
))

# Tissue groups present in the data that no body map can show. "thymus" simply
# does not exist in gganatogram; it is in neither hgMale_key nor hgFemale_key.
# Every organ that could stand in for "digestive system (general)" (stomach,
# colon, small_intestine) is already claimed by its own tissue group, and one
# organ in two groups would duplicate rows in the plot frame and shade the same
# shape with two different values.
#
# They are not dropped: make_plot_data() carries them through with organ = NA so
# they still reach the tables, and the About tab says why they are never drawn.
ORGAN_LESS_GROUPS <- c("thymus", "digestive system (general)")

#' Organ map for a sex.
organ_map_for <- function(sex) {
    if (identical(sex, "male")) MALE_ORGAN_MAP else FEMALE_ORGAN_MAP
}

#' gganatogram's organ ordering for a sex.
#'
#' Drawing in this order puts smaller structures on top of the larger ones that
#' would otherwise occlude them.
organ_levels_for <- function(sex) {
    if (identical(sex, "male")) hgMale_key$organ else hgFemale_key$organ
}
