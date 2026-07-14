from ehrql import create_dataset, days, months, years, case, when, create_measures, INTERVAL, minimum_of, maximum_of, get_parameter
from ehrql.tables.tpp import patients, medications, practice_registrations, clinical_events, apcs, addresses, ons_deaths, appointments
from datetime import date, datetime
import codelists_ehrQL as codelists
from analysis.dataset_definition_incidence import dataset, preceding_registration

# Read parameters from project.yaml
studystart_date = get_parameter("studystart_date")
studyend_date = get_parameter("studyend_date")
studyfup_date = get_parameter("studyfup_date")
diseases_list = get_parameter("diseases_list")
disease = get_parameter("measure_disease")
measure_start_date = get_parameter("measure_start_date")
intervals = int(get_parameter("measure_intervals"))

intervals_years = int(intervals/12)

interval_start = INTERVAL.start_date
interval_end = INTERVAL.end_date

# Currently registered with a practice
curr_registered = practice_registrations.for_patient_on(interval_start).exists_for_patient()

# Registration for at least 12 months before index date
preceding_reg_int = preceding_registration(interval_start).exists_for_patient()

# Practice region and pseudoid at interval start
region_int = preceding_registration(interval_start).practice_nuts1_region_name
practice_id_int = preceding_registration(interval_start).practice_pseudo_id

# IMD quintile at interval start
address_per_patient_int = addresses.for_patient_on(interval_start)
imd_rounded_int = address_per_patient_int.imd_rounded
imd_quintile_int = case(
    when((imd_rounded_int >= 0) & (imd_rounded_int < int(32844 * 1 / 5))).then("1 (most deprived)"),
    when(imd_rounded_int < int(32844 * 2 / 5)).then("2"),
    when(imd_rounded_int < int(32844 * 3 / 5)).then("3"),
    when(imd_rounded_int < int(32844 * 4 / 5)).then("4"),
    when(imd_rounded_int < int(32844 * 5 / 5)).then("5 (least deprived)"),
    otherwise="Unknown",
)

# Age at interval start
age = patients.age_on(interval_start)

age_band = case(  
    when((age >= 18) & (age <= 29)).then("age_18_29"),
    when((age >= 30) & (age <= 39)).then("age_30_39"),
    when((age >= 40) & (age <= 49)).then("age_40_49"),
    when((age >= 50) & (age <= 59)).then("age_50_59"),
    when((age >= 60) & (age <= 69)).then("age_60_69"),
    when((age >= 70) & (age <= 79)).then("age_70_79"),
    when((age >= 80)).then("age_greater_equal_80"),
)

measures = create_measures()
measures.configure_dummy_data(population_size=10000, legacy=True)
measures.configure_disclosure_control(enabled=False)
measures.define_defaults(intervals=months(intervals).starting_on(measure_start_date))

# Prevalence denominator
prev_denominator = (
    ((age >= 18) & (age <= 110))
    & dataset.sex.is_in(["male", "female"])
    & (dataset.date_of_death.is_after(interval_start) | dataset.date_of_death.is_null())
    & curr_registered
)

# Dictionaries to store values
prev = {}
prev_numerators = {} 
inc_case = {}
inc_case_12m_alive = {} 
incidence_numerators = {}
incidence_denominators = {}

# Prevalent diagnosis (at interval start)
prev[disease + "_prev"] = (
    (getattr(dataset, disease + "_prev_date") < interval_start)
).when_null_then(False)

# Prevalence numerator - people registered for more than one year on index date who have an diagnostic code on or before index date
prev_numerators[disease + "_prev_num"] = (
    prev[disease + "_prev"] & prev_denominator
)

# Incident case (incident date within interval window)
inc_case[disease + "_inc_case"] = (
    (getattr(dataset, disease + "_inc_date")).is_on_or_between(interval_start, interval_end)
).when_null_then(False)

# Preceding registration and alive at incident diagnosis date
inc_case_12m_alive[disease + "_inc_case_12m_alive"] = ( 
    inc_case[disease + "_inc_case"]
    & getattr(dataset, disease + "_pre_reg")
    & getattr(dataset, disease + "_alive_inc")
).when_null_then(False)

# Incidence numerator - people with new diagnostic codes in the 1 month after index date who have 12m+ preceding registration and alive 
incidence_numerators[disease + "_inc_num"] = (
    inc_case_12m_alive[disease + "_inc_case_12m_alive"]
    & ((age >= 18) & (age <= 110))
    & dataset.sex.is_in(["male", "female"])
)

# Incidence denominator - people with 12m+ registration prior to index date who do not have a diagnostic code on or before index date
incidence_denominators[disease + "_inc_denom"] = (
    (~prev[disease + "_prev"])
    & ((age >= 18) & (age <= 110))
    & dataset.sex.is_in(["male", "female"])
    & (dataset.date_of_death.is_after(interval_start) | dataset.date_of_death.is_null())
    & preceding_reg_int
)

# Incidence by age and sex
measures.define_measure(
    name=disease + "_incidence",
    numerator=incidence_numerators[disease + "_inc_num"],
    denominator=incidence_denominators[disease + "_inc_denom"],
    group_by={
        "sex": dataset.sex,
        "age": age_band,  
    },
)

# Incidence by ethnicity
measures.define_measure(
    name=disease + "_inc_ethn",
    numerator=incidence_numerators[disease + "_inc_num"],
    denominator=incidence_denominators[disease + "_inc_denom"],
    group_by={
        "ethnicity": dataset.ethnicity,
    },
)

# Incidence by IMD quintile (at interval start)
measures.define_measure(
    name=disease + "_inc_imd",
    numerator=incidence_numerators[disease + "_inc_num"],
    denominator=incidence_denominators[disease + "_inc_denom"],
    group_by={
        "imd": imd_quintile_int,
    },
)

# Incidence by region (at interval start)
measures.define_measure(
    name=disease + "_inc_region",
    numerator=incidence_numerators[disease + "_inc_num"],
    denominator=incidence_denominators[disease + "_inc_denom"],
    group_by={
        "region": region_int,
    },
)

# Prevalence by age and sex
measures.define_measure(
    name=disease + "_prevalence",
    numerator=prev_numerators[disease + "_prev_num"],
    denominator=prev_denominator,
    intervals=years(intervals_years).starting_on(measure_start_date),
    group_by={
        "sex": dataset.sex,
        "age": age_band,  
    },
)

# Prevalence by ethnicity
measures.define_measure(
    name=disease + "_prev_ethn",
    numerator=prev_numerators[disease + "_prev_num"],
    denominator=prev_denominator,
    group_by={
        "ethnicity": dataset.ethnicity,
    },
)

# Prevalence by IMD quintile (at interval start)
measures.define_measure(
    name=disease + "_prev_imd",
    numerator=prev_numerators[disease + "_prev_num"],
    denominator=prev_denominator,
    group_by={
        "imd": imd_quintile_int,
    },
)

# Prevalence by region (at interval start)
measures.define_measure(
    name=disease + "_prev_region",
    numerator=prev_numerators[disease + "_prev_num"],
    denominator=prev_denominator,
    group_by={
        "region": region_int,
    },
)