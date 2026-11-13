from ehrql import create_dataset, days, months, case, when, get_parameter
from ehrql.tables.tpp import patients, medications, practice_registrations, clinical_events
import codelists_ehrQL as codelists

# Read parameters from project.yaml
studystart_date = get_parameter("studystart_date")
studyend_date = get_parameter("studyend_date")
studyfup_date = get_parameter("studyfup_date")
intervalstart_date = get_parameter("intervalstart_date")
intervalend_date = intervalstart_date + months(12) - days(1)

# Definitions
disease = "gout"
disease_codes = getattr(codelists, f"{disease}_prev_snomed")

medication = "ult"
medication_codes = getattr(codelists, f"{medication}_codes")

blood = "urate"
blood_codes = getattr(codelists, f"{blood}_codes")

dataset = create_dataset()
dataset.configure_dummy_data(population_size=1000)

# First primary care code during the 12-month consultation period
disease_events = clinical_events.where(
    clinical_events.snomedct_code.is_in(disease_codes)
)
consult_date = disease_events.where(
    clinical_events.date.is_on_or_between(intervalstart_date, intervalend_date)
).sort_by(clinical_events.date).first_for_patient().date

dataset.add_column(f"{disease}_consult_date", consult_date)

followup_end_date = consult_date + months(12)
dataset.followup_end_date = followup_end_date

# Practice registration
def practice_registration(dx_date):
    return practice_registrations.where(
        practice_registrations.start_date.is_on_or_before(dx_date)
    ).except_where(
        practice_registrations.end_date.is_on_or_before(dx_date)
    ).sort_by(
        practice_registrations.start_date,
        practice_registrations.end_date,
        practice_registrations.practice_pseudo_id,
    ).last_for_patient()

dataset.index_practice_id = practice_registration(intervalstart_date).practice_pseudo_id
dataset.consult_practice_id = practice_registration(consult_date).practice_pseudo_id
dataset.consult_reg_end_date = practice_registration(consult_date).end_date

# Demographics
dataset.age = patients.age_on(consult_date)
dataset.sex = patients.sex
dataset.date_of_death = patients.date_of_death

# Last prescription within six calendar months before a specified date
def medication_before_date(med_codelist, anchor_date):
    return medications.where(
        medications.dmd_code.is_in(med_codelist)
    ).where(
        medications.date.is_on_or_between(anchor_date - months(6), anchor_date - days(1))
    ).sort_by(medications.date).last_for_patient().date

dataset.add_column(
    f"{medication}_bl_date",
    medication_before_date(medication_codes, consult_date),
)

# First prescription on or after a specified date
def medication_after_date(med_codelist, anchor_date):
    return medications.where(
        medications.dmd_code.is_in(med_codelist)
    ).where(
        medications.date.is_on_or_between(anchor_date, anchor_date + months(12))
        & (medications.date <= studyfup_date)
    ).sort_by(medications.date).first_for_patient().date

dataset.add_column(
    f"{medication}_after_date",
    medication_after_date(medication_codes, consult_date),
)

# Last blood test within 24 calendar months before a specified date
def blood_before_date(blood_codelist, anchor_date):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(blood_codelist)
    ).where(
        clinical_events.date.is_on_or_between(anchor_date - months(24), anchor_date)
    ).sort_by(clinical_events.date).last_for_patient()

baseline_blood = blood_before_date(blood_codes, consult_date)
dataset.add_column(f"{blood}_bl_date", baseline_blood.date)
dataset.add_column(f"{blood}_bl_value", baseline_blood.numeric_value)

# Recurrent blood tests within period after event
def recurrent_bloods(dx_codelist, anchor_date):
    return clinical_events.where(
        clinical_events.snomedct_code.is_in(dx_codelist)
    ).where(
        (clinical_events.date > anchor_date)
        & (clinical_events.date <= followup_end_date)
        & (clinical_events.date <= studyfup_date)
    ).where(
        clinical_events.numeric_value.is_not_null()
        & (clinical_events.numeric_value > 0)
        & (clinical_events.numeric_value < 3000)
    ).sort_by(
        clinical_events.date
    ).first_for_patient()

anchor = consult_date
for i in range(1, 20 + 1):
    blood_event = recurrent_bloods(blood_codes, anchor)
    dataset.add_column(f"{blood}_value_{i}", blood_event.numeric_value)
    dataset.add_column(f"{blood}_date_{i}", blood_event.date)
    dataset.add_column(
        f"{medication}_before_{blood}_{i}_date",
        medication_before_date(medication_codes, blood_event.date),
    )
    anchor = blood_event.date

# Additional date to identify histories exceeding 20 blood-test days
dataset.add_column(
    f"{blood}_overflow_date",
    recurrent_bloods(blood_codes, anchor).date,
)

# Study population: incident or prevalent gout consulting within the accrual period
dataset.define_population(
    consult_date.is_not_null()
    & (dataset.age >= 18) & (dataset.age <= 110)
    & dataset.sex.is_in(["male", "female"])
    & practice_registration(consult_date).exists_for_patient()
    & (patients.date_of_death.is_null() | (patients.date_of_death > consult_date))
)