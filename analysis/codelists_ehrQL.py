from ehrql import codelist_from_csv

# Demographic codelists
ethnicity_codes = codelist_from_csv(
    "codelists/opensafely-ethnicity-snomed-0removed.csv",
    column="code",
    category_column="Grouping_6",
)

# Primary diagnosis (incident disease codelist)
gout_snomed = codelist_from_csv(
    "codelists/user-markdrussell-gout.csv", column="code",
)

gout_icd = codelist_from_csv(
    "codelists/user-markdrussell-gout-admissions.csv", column="code",
)

# Primary diagnosis (prevalent disease codelist)
gout_prev_snomed = codelist_from_csv(
    "codelists/user-markdrussell-gout-prevalent.csv", column="code",
)

gout_prev_icd = codelist_from_csv(
    "codelists/user-markdrussell-gout-admissions.csv", column="code",
)

# Comorbidities
chd_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-chd_cod.csv", column="code",
)

diabetes_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-dmtype2audit_cod.csv", column="code",
)

stroke_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-strk_cod.csv", column="code",
)

tia_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-tia_cod.csv", column="code",
)

cva_codes = (
    stroke_codes + tia_codes
)

ckd_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-ckdatrisk2_cod.csv", column="code",
)

depression_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-depr_cod.csv", column="code",
)

heart_failure_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-hflvsd_cod.csv", column="code",
)

liver_disease_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-cldatrisk1_cod.csv", column="code",
)

transplant_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-orgtransp_cod.csv", column="code",
)

alcohol_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-excessalc_cod.csv", column="code",
)

hypertension_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-hyp_cod.csv", column="code",
)

# BMI and smoking status
bmi_codes = ["60621009", "846931000000101"]

clear_smoking_codes = codelist_from_csv(
    "codelists/opensafely-smoking-clear.csv",
    column="CTV3Code",
    category_column="Category",
)

# Disease-specific characteristics and relevant clinical events
tophi_codes = codelist_from_csv(
    "codelists/user-markdrussell-gouty-tophi.csv", column="code",
)

chronic_gout_codes = codelist_from_csv(
    "codelists/user-markdrussell-gout-chronic-arthritis.csv", column="code",
)

flare_codes = codelist_from_csv(
    "codelists/user-markdrussell-gout-flaresattacks.csv", column="code",
)

gout_admission_codes = codelist_from_csv(
    "codelists/user-markdrussell-gout-admissions.csv", column="code",
)

# Blood tests
urate_codes = codelist_from_csv(
    "codelists/user-markdrussell-serum-urate.csv", column="code",
)

creatinine_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-cre_cod.csv", column="code",
)

cholesterol_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-chol2_cod.csv", column="code",
)

hba1c_codes = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-ifcchbam_cod.csv", column="code",
)

# Medications
allopurinol_codes = codelist_from_csv(       
    "codelists/user-markdrussell-allopurinol-dmd.csv", column="code"
)

allopurinol_high_codes = codelist_from_csv(       
    "codelists/user-markdrussell-allopurinol-300mg-doses-dmd.csv", column="code"
)

febuxostat_codes = codelist_from_csv(       
    "codelists/user-markdrussell-febuxostat-dmd.csv", column="code"
)

febuxostat_high_codes = codelist_from_csv(       
    "codelists/user-markdrussell-febuxostat-120mg-doses-dmd.csv", column="code"
)

benzbromarone_codes = codelist_from_csv(       
    "codelists/user-markdrussell-benzbromarone-dmd.csv", column="code"
)

probenecid_codes = codelist_from_csv(       
    "codelists/user-markdrussell-probenecid-dmd.csv", column="code"
)

ult_codes = (
    allopurinol_codes + febuxostat_codes + benzbromarone_codes + probenecid_codes
)

colchicine_codes = codelist_from_csv(       
    "codelists/user-markdrussell-colchicine-dmd.csv", column="code"
)

steroid_codes = codelist_from_csv(       
    "codelists/nhs-drug-refsets-c19corstedrug_cod.csv", column="code"
)

nsaid_codes = codelist_from_csv(       
    "codelists/nhs-drug-refsets-oralnsaiddrug_cod.csv", column="code"
)

diuretic_codes = codelist_from_csv(       
    "codelists/user-markdrussell-diuretics-dmd.csv", column="code"
)

sglt2_codes = codelist_from_csv(
    "codelists/nhs-drug-refsets-sglt2idrug_cod.csv", column="code"
)

ace_codes = codelist_from_csv(
    "codelists/nhs-drug-refsets-ace_cod.csv", column="code"
)

arb_codes = codelist_from_csv(
    "codelists/nhs-drug-refsets-aii_cod.csv", column="code"
)

ace_arb_codes = (
    ace_codes + arb_codes
)


# Referral codes
rheumatology_referral = codelist_from_csv(
    "codelists/user-markdrussell-referral-to-rheumatology-only.csv", column = "code"
)

# Outpatient appointment codes (treatment function code for specialty)
rheumatology_outpatient = ["410"]