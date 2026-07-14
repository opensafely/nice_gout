from ehrql import create_dataset, months, years, case, when, minimum_of, maximum_of, get_parameter
from ehrql.tables.tpp import patients, practice_registrations, clinical_events, apcs, addresses, ethnicity_from_sus
from ehrql.codes import ICD10Code
from datetime import date, datetime
from functools import reduce
import codelists_ehrQL as codelists

# Read parameters from project.yaml
studystart_date = get_parameter("studystart_date")
studyend_date = get_parameter("studyend_date")
studyfup_date = get_parameter("studyfup_date")
diseases_list = get_parameter("diseases_list")

diseases = diseases_list if isinstance(diseases_list, list) else [diseases_list]
print("Diseases:", diseases)

# Define codelist types
codelist_types = ["snomed", "icd"]

# Any practice registration before study end date
any_registration = practice_registrations.where(
        practice_registrations.start_date <= studyend_date
    ).except_where(
        practice_registrations.end_date < studystart_date    
    ).exists_for_patient()

# Practice registration for at least 12 months prior to diagnosis date
def preceding_registration(dx_date):
    return practice_registrations.where(
        practice_registrations.start_date.is_on_or_before(dx_date - months(12))
    ).except_where(
        practice_registrations.end_date.is_on_or_before(dx_date)
    ).sort_by(
        practice_registrations.start_date,
        practice_registrations.end_date,
        practice_registrations.practice_pseudo_id,
    ).last_for_patient()

def create_dataset_with_variables():
    dataset = create_dataset()
    dataset.configure_dummy_data(population_size=10000)

    # Incident diagnostic code in primary care record (SNOMED), assuming before study end date
    def first_code_in_period_snomed(dx_codelist):
        return clinical_events.where(
            clinical_events.snomedct_code.is_in(dx_codelist)
        ).where(
            clinical_events.date.is_on_or_before(studyend_date)
        ).sort_by(
            clinical_events.date
        ).first_for_patient()

    # Incident diagnostic code in secondary care record (ICD10 primary diagnoses), assuming before study end date
    def first_code_in_period_icd(dx_codelist):
        return apcs.where(
            apcs.primary_diagnosis.is_in(dx_codelist)
        ).where(
            apcs.admission_date.is_on_or_before(studyend_date)
        ).sort_by(
            apcs.admission_date
        ).first_for_patient()

    # Count of diagnostic codes in primary care record
    def count_code_in_period_snomed(dx_codelist):
        return clinical_events.where(
            clinical_events.snomedct_code.is_in(dx_codelist)
        ).where(
            clinical_events.date.is_on_or_before(studyfup_date)
        ).except_where(
            clinical_events.date.is_before(studystart_date)
        ).count_for_patient()

    # Count of diagnostic codes in secondary care record
    def count_code_in_period_icd(dx_codelist):
        return apcs.where(
            apcs.primary_diagnosis.is_in(dx_codelist)
        ).where(
            apcs.admission_date.is_on_or_before(studyfup_date)
        ).except_where(
            apcs.admission_date.is_before(studystart_date)
        ).count_for_patient()

    # Expand 3-character ICD10 codes
    def expand_three_char_icd10_codes(dx_codelist):
        return dx_codelist + [f"{code}X" for code in dx_codelist if len(code) == 3]
    
    # Define sex
    dataset.sex = patients.sex

    # Date of death
    dataset.date_of_death = patients.date_of_death

    # Define patient ethnicity
    latest_ethnicity_code = (
        clinical_events.where(clinical_events.snomedct_code.is_in(codelists.ethnicity_codes))
        .where(clinical_events.date.is_on_or_before(studyend_date))
        .sort_by(clinical_events.date)
        .last_for_patient().snomedct_code.to_category(codelists.ethnicity_codes)
    )

    # Extract ethnicity from SUS records if it isn't present in primary care data
    ethnicity_sus = ethnicity_from_sus.code

    dataset.ethnicity = case(
        when((latest_ethnicity_code == "1") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["A", "B", "C"])))).then("White"),
        when((latest_ethnicity_code == "2") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["D", "E", "F", "G"])))).then("Mixed"),
        when((latest_ethnicity_code == "3") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["H", "J", "K", "L"])))).then("Asian or Asian British"),
        when((latest_ethnicity_code == "4") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["M", "N", "P"])))).then("Black or Black British"),
        when((latest_ethnicity_code == "5") | ((latest_ethnicity_code.is_null()) & (ethnicity_sus.is_in(["R", "S"])))).then("Chinese or Other Ethnic Groups"),
        otherwise="Unknown",
    )
    
    # Identify incident diagnosis date
    for disease in diseases:

        # Incident codes
        for codelist_type in codelist_types:

            if (f"{codelist_type}" == "snomed"):
                if hasattr(codelists, f"{disease}_snomed"):
                    disease_codelist = getattr(codelists, f"{disease}_snomed")
                    dataset.add_column(f"{disease}_prim_date", first_code_in_period_snomed(disease_codelist).date)
                    dataset.add_column(f"{disease}_prim_count", count_code_in_period_snomed(disease_codelist))
                else:
                    dataset.add_column(f"{disease}_prim_date", first_code_in_period_snomed([]).date)
                    dataset.add_column(f"{disease}_prim_count", count_code_in_period_snomed([]))
            elif (f"{codelist_type}" == "icd"):
                if hasattr(codelists, f"{disease}_icd"):
                    disease_codelist = getattr(codelists, f"{disease}_icd")
                    disease_codelist = expand_three_char_icd10_codes(disease_codelist)
                    dataset.add_column(f"{disease}_sec_date", first_code_in_period_icd(disease_codelist).admission_date)
                    dataset.add_column(f"{disease}_sec_count", count_code_in_period_icd(disease_codelist))
                else:
                    dataset.add_column(f"{disease}_sec_date", first_code_in_period_icd([]).admission_date)
                    dataset.add_column(f"{disease}_sec_count", count_code_in_period_icd([]))
            else:
                dataset.add_column(f"{disease}_{codelist_type}_inc_date", None)

        # Incident date for each disease (combined primary and secondary care)
        dataset.add_column(f"{disease}_inc_date",
            minimum_of(*[date for date in [
                (getattr(dataset, f"{disease}_prim_date", None)),
                (getattr(dataset, f"{disease}_sec_date", None))
                ] if date is not None]),
        )

        # Incident date within window (combined primary and secondary care)
        dataset.add_column(f"{disease}_inc_case",
            (getattr(dataset, disease + "_inc_date").is_on_or_between(studystart_date, studyend_date)
            ).when_null_then(False)
        )

        # 12 months registration preceding incident diagnosis date (combined primary and secondary care)
        dataset.add_column(f"{disease}_pre_reg",
            preceding_registration(getattr(dataset, f"{disease}_inc_date")
            ).exists_for_patient()
        )

        # Age at diagnosis (combined primary and secondary care)
        dataset.add_column(f"{disease}_age",
            (patients.age_on(getattr(dataset, f"{disease}_inc_date"))
            )
        )

        # Alive at incident diagnosis date (combined primary and secondary care)
        dataset.add_column(f"{disease}_alive_inc",
            ((dataset.date_of_death.is_after(getattr(dataset, f"{disease}_inc_date"))) | dataset.date_of_death.is_null()
            ).when_null_then(False)
        )

        # Incident date within window (primary care only)
        dataset.add_column(f"{disease}_inc_case_p",
            (getattr(dataset, disease + "_prim_date").is_on_or_between(studystart_date, studyend_date)
            ).when_null_then(False)
        )

        # 12 months registration preceding incident diagnosis date (primary care only)
        dataset.add_column(f"{disease}_pre_reg_p",
            preceding_registration(getattr(dataset, f"{disease}_prim_date")
            ).exists_for_patient()
        )

        # Age at diagnosis (primary care only)
        dataset.add_column(f"{disease}_age_p",
            (patients.age_on(getattr(dataset, f"{disease}_prim_date"))
            )
        )

        # Alive at incident diagnosis date (primary care only)
        dataset.add_column(f"{disease}_alive_inc_p",
            ((dataset.date_of_death.is_after(getattr(dataset, f"{disease}_prim_date"))) | dataset.date_of_death.is_null()
            ).when_null_then(False)
        )

        # Prevalent diagnoses
        for codelist_type in codelist_types:

            if (f"{codelist_type}" == "snomed"):
                if hasattr(codelists, f"{disease}_prev_snomed"):
                    disease_codelist = getattr(codelists, f"{disease}_prev_snomed")
                    dataset.add_column(f"{disease}_prev_prim_date", first_code_in_period_snomed(disease_codelist).date)
                else:
                    dataset.add_column(f"{disease}_prev_prim_date", first_code_in_period_snomed([]).date)
            elif (f"{codelist_type}" == "icd"):
                if hasattr(codelists, f"{disease}_prev_icd"):
                    disease_codelist = getattr(codelists, f"{disease}_prev_icd")
                    disease_codelist = expand_three_char_icd10_codes(disease_codelist)
                    dataset.add_column(f"{disease}_prev_sec_date", first_code_in_period_icd(disease_codelist).admission_date)
                else:
                    dataset.add_column(f"{disease}_prev_sec_date", first_code_in_period_icd([]).admission_date)
            else:
                dataset.add_column(f"{disease}_{codelist_type}_prev_date", None)

        dataset.add_column(f"{disease}_prev_date",
            minimum_of(*[date for date in [
                (getattr(dataset, f"{disease}_prev_prim_date", None)),
                (getattr(dataset, f"{disease}_prev_sec_date", None))
                ] if date is not None]),
        )

    return dataset

def get_population(dataset):
    # Create variable for anyone with at least one diagnostic code
    any_inc_case = reduce(lambda x, y: x | y, [
        getattr(dataset, f"{d}_inc_case") for d in diseases
    ])

    # Define population as any patient with at least one diagnostic code, registered after index date - then apply further restrictions later (age, death and preceding registration)
    return (any_inc_case
        & any_registration 
        & dataset.sex.is_in(["male", "female"]))

dataset = create_dataset_with_variables()
dataset.define_population(get_population(dataset))  