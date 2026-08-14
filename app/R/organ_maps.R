# Mapping from the data's tissue groups to the anatomical organs gganatogram
# knows how to draw. One tissue group fans out to several organs.
#
# The male and female maps used to be maintained as two independent lists that
# were 90% identical; they are built from a shared base here so the common
# entries cannot drift apart.

# Tissue groups whose organ set is the same for both sexes.
BASE_ORGAN_MAP <- list(
    "blood" = "leukocyte",
    "epithelium and mucosal tissues" = "pleura",
    "large intestine" = c("caecum", "rectum", "colon"),
    "respiratory system" = c("lung", "bronchus", "pulmonary_valve", "diaphragm"),
    "spleen" = "spleen",
    "gastrointestinal accessory organs" = c("liver", "pancreas", "gall_bladder"),
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

# NOTE: "tongue" (male only) and "amygdala" (male only) are carried over
# verbatim from the original lists. gganatogram does provide female versions of
# both, so this asymmetry may be an oversight rather than a deliberate choice --
# left as-is to preserve existing output.
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
