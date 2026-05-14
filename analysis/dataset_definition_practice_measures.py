from ehrql import create_dataset, days, months, years, case, when, create_measures, INTERVAL, minimum_of, maximum_of, get_parameter
from ehrql.tables.tpp import patients, medications, practice_registrations, clinical_events, apcs, addresses, ons_deaths, appointments
from datetime import date, datetime
import codelists_ehrQL as codelists
from analysis.dataset_definition_incidence import dataset

# Read parameters from project.yaml
studystart_date = get_parameter("studystart_date")
studyend_date = get_parameter("studyend_date")
studyfup_date = get_parameter("studyfup_date")
diseases_list = get_parameter("diseases_list")
disease = get_parameter("measure_disease")
measure_start_date = get_parameter("measure_start_date")
intervals = int(get_parameter("measure_intervals"))

interval_start = INTERVAL.start_date
interval_end = INTERVAL.end_date

# Currently registered with a practice
curr_registered = practice_registrations.for_patient_on(interval_start).exists_for_patient()

# Practice pseudoid at interval start
practice_id_int = practice_registrations.for_patient_on(interval_start).practice_pseudo_id

# Age at interval start
age = patients.age_on(interval_start)

measures = create_measures()
measures.configure_dummy_data(population_size=10000)
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

# Prevalent diagnosis (at interval start)
prev[disease + "_prev"] = (
    (getattr(dataset, disease + "_prev_date") < interval_start)
).when_null_then(False)

# Prevalence numerator - people registered for more than one year on index date who have an diagnostic code on or before index date
prev_numerators[disease + "_prev_num"] = (
    prev[disease + "_prev"] & prev_denominator
)

# Prevalence and list size by practice
measures.define_measure(
    name=disease + "_prev_practice",
    numerator=prev_numerators[disease + "_prev_num"],
    denominator=prev_denominator,
    intervals=years(intervals).starting_on(measure_start_date),
    group_by={
        "practice_id": practice_id_int,
    },
)