from ehrql import create_dataset, days, months, years, case, when, get_parameter
from ehrql.tables.tpp import patients, medications, practice_registrations, clinical_events, apcs, addresses, ethnicity_from_sus, opa, emergency_care_attendances
from ehrql.codes import ICD10Code
from datetime import date, datetime
import codelists_ehrQL as codelists
import json
from analysis.dataset_definition_incidence import create_dataset_with_variables, get_population, preceding_registration

dataset = create_dataset_with_variables()

# Read parameters from project.yaml
studystart_date = get_parameter("studystart_date")
studyend_date = get_parameter("studyend_date")
studyfup_date = get_parameter("studyfup_date")
primary_disease = get_parameter("primary_disease")

# Read other parameters and ensure they are lists
def ensure_list(x):
    if isinstance(x, list):
        return x
    if x is None:
        return []
    if isinstance(x, str):
        s = x.strip()
        if not s:
            return []
        if s[0] == "[" and s[-1] == "]":
            try:
                v = json.loads(s)
                if isinstance(v, list):
                    return [str(t).strip() for t in v if str(t).strip()]
            except Exception:
                s = s[1:-1]
        if "," in s:
            return [p.strip().strip('\'"') for p in s.split(",") if p.strip().strip('\'"')]
        if " " in s:
            return [p for p in s.split() if p]
        return [s]
    s = str(x).strip()
    return [s] if s else []

diseases_list = ensure_list(get_parameter("diseases_list"))
comorbidities_list = ensure_list(get_parameter("comorbidities_list"))
disease_features_list = ensure_list(get_parameter("disease_features_list"))
events_list = ensure_list(get_parameter("events_list"))
admissions_list = ensure_list(get_parameter("admissions_list"))
bloods_list = ensure_list(get_parameter("bloods_list"))
medications_list = ensure_list(get_parameter("medications_list"))
outpatients_list = ensure_list(get_parameter("outpatients_list"))

# Store date of primary diagnosis
dx_date = getattr(dataset, primary_disease + "_inc_date")

# Function to identify code on or before end of follow-up period (can define as first or last code in period)
def code_before_studyend(dx_codelist):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        clinical_events.date.is_on_or_before(studyfup_date)
    ).sort_by(
        clinical_events.date
    )

# Function to identify last code in primary care up to X [specify] months before primary diagnosis date
def code_before_diagnosis(dx_codelist, pre_time_window):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        clinical_events.date.is_on_or_between((dx_date - months(pre_time_window)), dx_date)
    ).sort_by(
        clinical_events.date
    ).last_for_patient()

# Function to identify recurrent clinical events in primary care after diagnosis (separated by at least X days of lag)
def recurrent_events(dx_codelist, anchor_date, lag):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        (clinical_events.date > (anchor_date + days(lag))) & (clinical_events.date <= studyfup_date)
    ).sort_by(
        clinical_events.date
    ).first_for_patient()

# Function to identify admissions to secondary care after diagnosis date (using primary admission diagnosis code; separated by at least X days of lag)
def admission_events(dx_codelist, anchor_date, lag):
    return apcs.where(
        apcs.primary_diagnosis.is_in(dx_codelist)
    ).where(
        (apcs.admission_date > (anchor_date + days(lag))) & (apcs.admission_date <= studyfup_date)
    ).sort_by(
        apcs.admission_date
    ).first_for_patient()

# For all primary admission events need to expand 3-character ICD10 codes
def expand_three_char_icd10_codes(dx_codelist):
    return dx_codelist + [f"{code}X" for code in dx_codelist if len(code) == 3]

# Function to identify emergency department attendances after diagnosis date (using SNOMED diagnosis code in first position only; separated by at least X days of lag)
def ed_attendance_events(dx_codelist, anchor_date, lag):
    return emergency_care_attendances.where(
        emergency_care_attendances.diagnosis_01.is_in(dx_codelist)
    ).where(
        (emergency_care_attendances.arrival_date > (anchor_date + days(lag))) & (emergency_care_attendances.arrival_date <= studyfup_date)
    ).sort_by(
        emergency_care_attendances.arrival_date
    ).first_for_patient()

# Function to identify recurrent blood tests after a specified time point (only those with associated numeric values)
def recurrent_bloods(dx_codelist, anchor_date):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        (clinical_events.date > anchor_date) & (clinical_events.date <= studyfup_date)
    ).where(
        clinical_events.numeric_value.is_not_null() & (clinical_events.numeric_value > 0) & (clinical_events.numeric_value < 3000)
    ).sort_by(
        clinical_events.date
    ).first_for_patient()

# Function to identify prescription issued on or before study end date
def medication_dates(dmd_codelist):
    return medications.where(
        medications.dmd_code.is_in(dmd_codelist)
    ).where(
        medications.date.is_on_or_before(studyfup_date)
    ).sort_by(
        medications.date
    )

# Function to identify last prescription issued up to X (specify) months before diagnosis
def medication_diag_dates(dmd_codelist, pre_time_window):
    return medications.where(
        medications.dmd_code.is_in(dmd_codelist)
    ).where(
        medications.date.is_on_or_between((dx_date - months(pre_time_window)), dx_date)
    ).sort_by(
        medications.date
    ).last_for_patient()

# Function to identify prescription count for a drug within X (specify) months of diagnosis
def medication_count(dmd_codelist, time_window):
    return medications.where(
        medications.dmd_code.is_in(dmd_codelist)
    ).where(
        medications.date.is_on_or_between(dx_date, (dx_date + months(time_window))) & medications.date.is_on_or_before(studyfup_date)
    ).sort_by(
        medications.date
    ).count_for_patient()

# Function to identify recurrent prescriptions in primary care after diagnosis (separated by at least X days of lag)
def recurrent_meds(dmd_codelist, anchor_date, lag):
    return medications.where(
        medications.dmd_code.is_in(dmd_codelist)
    ).where(
        (medications.date > (anchor_date + days(lag))) & (medications.date <= studyfup_date)
    ).sort_by(
        medications.date
    ).first_for_patient()

# Function to identify first specialty outpatient appointment in the X (specify) months before before diagnosis and up to study end date
def appointment_dates(tfc_codelist, time_window):
    return opa.where(
        opa.treatment_function_code.is_in(tfc_codelist)
    ).where(
        opa.appointment_date.is_on_or_between((dx_date - months(time_window)), studyfup_date)
    ).sort_by(
        opa.appointment_date
    ).first_for_patient()

# Function to identify first specialty referral in the X (specify) months before diagnosis and up to study end date
def referral_dates(ref_codelist, time_window):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(ref_codelist)
    ).where(
        clinical_events.date.is_on_or_between((dx_date - months(time_window)), studyfup_date)
    ).sort_by(
        clinical_events.date
    ).first_for_patient()

# Practice pseudoid, region and registration end date (most recent practice at time of primary diagnosis)
dataset.practice_id = preceding_registration(dx_date).practice_pseudo_id
dataset.region = preceding_registration(dx_date).practice_nuts1_region_name
dataset.reg_end_date = preceding_registration(dx_date).end_date

# IMD at time of primary diagnosis
address_per_patient = addresses.for_patient_on(dx_date)
imd_rounded = address_per_patient.imd_rounded
dataset.imd_quintile = case(
    when((imd_rounded >= 0) & (imd_rounded < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="Unknown",
)

# Date of diagnosis for comorbidities (first recorded code before study end date)
for comorbidity in comorbidities_list:
    comorbidity_codelist = getattr(codelists, f"{comorbidity}_codes")
    dataset.add_column(f"{comorbidity}_date", code_before_studyend(comorbidity_codelist).first_for_patient().date)

# Body mass index (last recorded code up to 120 [specify] months before diagnosis)
dataset.bmi_value = code_before_diagnosis(codelists.bmi_codes, 120).numeric_value
dataset.bmi_date = code_before_diagnosis(codelists.bmi_codes, 120).date

# Smoking status (last recorded code before primary diagnosis)
dataset.most_recent_smoking_code=clinical_events.where(
        clinical_events.ctv3_code.is_in(codelists.clear_smoking_codes)
    ).where(
        clinical_events.date.is_on_or_before(dx_date)
    ).sort_by(
        clinical_events.date
    ).last_for_patient().ctv3_code.to_category(codelists.clear_smoking_codes)

def filter_codes_by_category(codelist, include):
    return {k:v for k,v in codelist.items() if v in include}

dataset.ever_smoked=clinical_events.where(
        clinical_events.ctv3_code.is_in(filter_codes_by_category(codelists.clear_smoking_codes, include=["S", "E"]))
    ).where(
        clinical_events.date.is_on_or_before(dx_date)
    ).exists_for_patient()

dataset.smoking_status=case(
    when(dataset.most_recent_smoking_code == "S").then("S"),
    when((dataset.most_recent_smoking_code == "E") | ((dataset.most_recent_smoking_code == "N") & (dataset.ever_smoked == True))).then("E"),
    when((dataset.most_recent_smoking_code == "N") & (dataset.ever_smoked == False)).then("N"),
    otherwise="M"
)

# Disease-specific characteristics: first recorded code in primary care record
for disease_feature in disease_features_list:
    disease_feature_codelist = getattr(codelists, f"{disease_feature}_codes")
    dataset.add_column(f"{disease_feature}_date", code_before_studyend(disease_feature_codelist).first_for_patient().date)

# Disease-related clinical events: first X [specify] events after diagnosis, separated by at least 14 [specify] days
for event in events_list:
    event_codelist = getattr(codelists, f"{event}_codes")
    anchor = dx_date
    for i in range(1, 50 + 1):
        event_date = recurrent_events(event_codelist, anchor, 14).date
        dataset.add_column(f"{event}_date_{i}", event_date)
        anchor = event_date

# Recurrent primary diagnosis codes/consults: first X [specify] recorded diagnostic codes after diagnosis, separated by at least 1 [specify] days
primarydx_codelist = getattr(codelists, f"{primary_disease}_snomed")
anchor = dx_date
for i in range(1, 50 + 1):
    consult_date = recurrent_events(primarydx_codelist, anchor, 1).date
    dataset.add_column(f"{primary_disease}_cons_date_{i}", consult_date)
    anchor = consult_date        

# Disease-related admission events: first X [specify] events after diagnosis, separated by at least 14 [specify] days
for diagnosis in admissions_list:
    admission_codelist = getattr(codelists, f"{diagnosis}_admission_codes")
    admission_codelist = expand_three_char_icd10_codes(admission_codelist)   
    anchor = dx_date
    for i in range(1, 50 + 1):
        admission_date = admission_events(admission_codelist, anchor, 14).admission_date
        dataset.add_column(f"{diagnosis}_adm_date_{i}", admission_date)
        anchor = admission_date

# Disease-related emergency department attendances: first X [specify] events after diagnosis, separated by at least 14 [specify] days
for diagnosis in admissions_list:
    emergency_codelist = getattr(codelists, f"{diagnosis}_snomed")
    anchor = dx_date
    for i in range(1, 50 + 1):
        emergency_date = ed_attendance_events(emergency_codelist, anchor, 14).arrival_date
        dataset.add_column(f"{diagnosis}_ed_date_{i}", emergency_date)
        anchor = emergency_date    

# Blood tests: first X [specify] recorded blood tests from 24 [specify] months before diagnosis date
for blood in bloods_list:
    blood_codelist = getattr(codelists, f"{blood}_codes")
    anchor = (dx_date - months(24))
    for i in range(1, 50 + 1):
        blood_event = recurrent_bloods(blood_codelist, anchor)
        dataset.add_column(f"{blood}_value_{i}", blood_event.numeric_value)
        dataset.add_column(f"{blood}_date_{i}", blood_event.date)
        anchor = blood_event.date

# Medications at baseline and after diagnosis
for medication in medications_list:
    medication_codelist = getattr(codelists, f"{medication}_codes")
    ## First prescription date before study end date
    dataset.add_column(f"{medication}_first_date", medication_dates(medication_codelist).first_for_patient().date)
    ## Last prescription date before study end date
    dataset.add_column(f"{medication}_last_date", medication_dates(medication_codelist).last_for_patient().date)
    ## Last prescription in the X (specify) months before diagnosis date
    dataset.add_column(f"{medication}_bl_date", medication_diag_dates(medication_codelist, 6).date)
    ## Count of prescriptions within X (specify) months after diagnosis
    dataset.add_column(f"{medication}_count_6m", medication_count(medication_codelist, 6))
    dataset.add_column(f"{medication}_count_12m", medication_count(medication_codelist, 12))
    ## First X [specify] prescriptions after diagnosis, separated by at least 1 [specify] day; do further cleaning later in pipeline
    anchor = dx_date
    for i in range(1, 50 + 1):
        med_date = recurrent_meds(medication_codelist, anchor, 1).date
        dataset.add_column(f"{medication}_date_{i}", med_date)
        anchor = med_date

# Outpatient appointments and referrals in relevant specialties (first recorded appointment/referral from [specify] months before diagnosis to study end)
for specialty in outpatients_list:
    ## Outpatient appointments
    outpatient_codelist = getattr(codelists, f"{specialty}_outpatient")
    dataset.add_column(f"{specialty}_opa_date", appointment_dates(outpatient_codelist, 12).appointment_date)

    ## Secondary care referrals
    referral_codelist = getattr(codelists, f"{specialty}_referral")
    dataset.add_column(f"{specialty}_ref_date", referral_dates(referral_codelist, 12).date)

# Study definition criteria
incidence_dataset_population = get_population(dataset)

# Define study population (has primary diagnosis of interest; aged 18-110 at diagnosis; registered with practice for at least 12 months before diagnosis; alive at diagnosis)
dataset.define_population(
    incidence_dataset_population &
    (getattr(dataset, primary_disease + "_inc_case")) &
    ((getattr(dataset, primary_disease + "_age") >= 18) & (getattr(dataset, primary_disease + "_age") <= 110)) &
    (getattr(dataset, primary_disease + "_pre_reg")) &
    (getattr(dataset, primary_disease + "_alive_inc"))
)