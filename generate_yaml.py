from datetime import datetime, date

# Define primary disease, comparison diseases (if any), comorbidities, demographic comparisons of interest, disease features, clinical events, blood tests and medications of interest (need to ensure each has a codelist)
primary_disease = "gout"
comparison_diseases = ""
demographic = "agegroup sex ethnicity imd region"
comorbidities = "chd diabetes cva ckd hypertension depression heart_failure liver_disease transplant alcohol"
disease_features = "tophi chronic_gout"
events = "flare"
admissions = "gout"
bloods = "urate creatinine cholesterol hba1c"
medications = "ult allopurinol allopurinol_high febuxostat febuxostat_high benzbromarone probenecid colchicine steroid nsaid diuretic sglt2 ace_arb"
outpatients = "rheumatology"

# Define study period dates
studystart_date = "2016-07-01"
studyend_date = "2025-06-30"
studyfup_date = "2026-06-30"

# Measure disease incidence (yes or no)
incidence = "yes"

# Define intervention date(s) of interest for intervention analyses
intervention_date_1 = "2020-03-01"
intervention_date_2 = "2022-06-01"

# Store start years, number of monthly intervals, and merge disease lists
start_dt = datetime.strptime(studystart_date, "%Y-%m-%d").date()
end_dt = datetime.strptime(studyend_date, "%Y-%m-%d").date()
start_year = start_dt.year
end_year = end_dt.year

# Disease lists for Python and R scripts
primary_list = [primary_disease] if primary_disease else []
comparison_list = [d.strip() for d in comparison_diseases.split(",")] if comparison_diseases else []
diseases_combined = primary_list + comparison_list
diseases_list = " ".join(diseases_combined)
demographic_str = demographic.split()
demographic_list = " ".join(demographic_str)
comorbidities_str = comorbidities.split()
comorbidities_list = " ".join(comorbidities_str)
disease_features_str = disease_features.split()
disease_features_list = " ".join(disease_features_str)
events_str = events.split()
events_list = " ".join(events_str)
admissions_str = admissions.split()
admissions_list = " ".join(admissions_str)
bloods_str = bloods.split()
bloods_list = " ".join(bloods_str)
medications_str = medications.split()
medications_list = " ".join(medications_str)
outpatients_str = outpatients.split()
outpatients_list = " ".join(outpatients_str)

# Disease lists for Stata scripts (pipe-separated)
diseases_list_stata = "|".join(diseases_combined)
demographic_list_stata = "|".join(demographic_str)
comorbidities_list_stata = "|".join(comorbidities_str)
disease_features_list_stata = "|".join(disease_features_str)
events_list_stata = "|".join(events_str)
admissions_list_stata = "|".join(admissions_str)
bloods_list_stata = "|".join(bloods_str)
medications_list_stata = "|".join(medications_str)
outpatients_list_stata = "|".join(outpatients_str)

yaml_header = f"""
version: '4.0'

actions:    
  generate_dataset_incidence:
    run: ehrql:v1 generate-dataset analysis/dataset_definition_incidence.py
      --output output/dataset_incidence.csv
      --
      --studystart_date "{studystart_date}"
      --studyend_date "{studyend_date}"
      --studyfup_date "{studyfup_date}"
      --diseases_list {diseases_list}
    outputs:
      highly_sensitive:
        cohort: output/dataset_incidence.csv          

  generate_dataset_primary:
    run: ehrql:v1 generate-dataset analysis/dataset_definition_primary.py 
      --output output/dataset_primary.csv
      --
      --studystart_date "{studystart_date}"
      --studyend_date "{studyend_date}"
      --studyfup_date "{studyfup_date}"
      --primary_disease {primary_disease}
      --diseases_list {diseases_list}
      --comorbidities_list {comorbidities_list}
      --disease_features_list {disease_features_list}
      --events_list {events_list}
      --admissions_list {admissions_list}
      --bloods_list {bloods_list}
      --medications_list {medications_list}
      --outpatients_list {outpatients_list}
    needs: [generate_dataset_incidence]
    outputs:
      highly_sensitive:
        cohort: output/dataset_primary.csv         
"""

yaml_template = ""
if incidence == "yes":
    yaml_template = """
  measures_dataset_{disease}_{year}:
    run: ehrql:v1 generate-measures analysis/dataset_definition_incidence_measures.py
      --output output/measures/measures_dataset_{disease}_{year}.csv
      --
      --studystart_date "{studystart_date}"
      --studyend_date "{studyend_date}"
      --studyfup_date "{studyfup_date}"
      --diseases_list {diseases_list}
      --measure_disease {disease}
      --measure_start_date "{year}-07-01"
      --measure_intervals {intervals}
    needs: [generate_dataset_incidence]
    outputs:
      highly_sensitive:
        measure_csv: output/measures/measures_dataset_{disease}_{year}.csv
"""

yaml_body = ""
all_needs = []

# Monthly intervals calculated for each study year (change these dates if study months change)
for year in range(start_year, end_year):
    win_start = date(year, 7, 1)
    win_end = date(year + 1, 6, 30)

    if win_end <= end_dt:
        intervals = 12
    else:
        intervals = (end_dt.year - year) * 12 + (end_dt.month - 7) + 1

    for dis in diseases_combined:
        yaml_body += yaml_template.format(
            studystart_date=studystart_date,
            studyend_date=studyend_date,
            studyfup_date=studyfup_date,
            diseases_list=diseases_list,
            disease=dis,
            year=year,
            intervals=intervals,
        )
        all_needs.append(f"measures_dataset_{dis}_{year}")

needs_list = ", ".join(all_needs)

yaml_incidence_template = ""
if incidence == "yes":
  yaml_incidence_template = f"""
  incidence_cleaning:
    run: stata-mp:latest analysis/001_incidence_cleaning.do "{diseases_list_stata}" "{studystart_date}" "{studyend_date}"
    needs: [generate_dataset_incidence, {needs_list}]
    outputs:
      moderately_sensitive:
        log1: logs/incidence_cleaning.log   
        table1: output/tables/redacted_counts_*.csv
        
  incidence_graphs:
    run: stata-mp:latest analysis/002_incidence_graphs.do "{diseases_list_stata}" "{intervention_date_1}"
    needs: [incidence_cleaning]
    outputs:
      moderately_sensitive:
        log1: logs/incidence_graphs.log
        table1: output/tables/arima_standardised.csv
        figure1: output/figures/inc_*.svg
        figure2: output/figures/prev_*.svg

  sarima:
    run: r:latest analysis/100_sarima.R "{intervention_date_1}"
    needs: [incidence_graphs]
    outputs:
      moderately_sensitive:
        log1: logs/sarima_log.txt   
        figure1: output/figures/auto_residuals_*.png
        figure2: output/figures/obs_pred_*.png
        table1: output/tables/change_incidence_byyear.csv
  """

yaml_incidence = yaml_incidence_template.format(needs_list=needs_list)

# Practice-level measures at study start; could make this yearly depending on computational requirements
yaml_practice = f"""
  measures_practice_{primary_disease}:
    run: ehrql:v1 generate-measures analysis/dataset_definition_practice_measures.py
      --output output/measures/measures_practice_{primary_disease}.csv
      --
      --studystart_date "{studystart_date}"
      --studyend_date "{studyend_date}"
      --studyfup_date "{studyfup_date}"
      --diseases_list {diseases_list}
      --measure_disease {primary_disease}
      --measure_start_date "{studystart_date}"
      --measure_intervals 1
    needs: [generate_dataset_incidence]
    outputs:
      highly_sensitive:
        measure_csv: output/measures/measures_practice_{primary_disease}.csv
"""

yaml_footer = f"""
  cohort_cleaning:
    run: stata-mp:latest analysis/200_cohort_cleaning.do "{primary_disease}" "{studystart_date}" "{studyend_date}" "{studyfup_date}" "{intervention_date_2}" "{demographic_list_stata}" "{comorbidities_list_stata}" "{disease_features_list_stata}" "{events_list_stata}" "{admissions_list_stata}" "{bloods_list_stata}" "{medications_list_stata}" "{outpatients_list_stata}"
    needs: [generate_dataset_primary, measures_practice_{primary_disease}]
    outputs:
      highly_sensitive:
        log1: logs/cohort_cleaning.log   
        data1: output/data/cohort_processed.dta
        data2: output/data/flares_long.dta

  data_tables:
    run: stata-mp:latest analysis/300_data_tables.do "{primary_disease}" "{demographic_list_stata}" "{outpatients_list_stata}"
    needs: [cohort_cleaning]
    outputs:
      moderately_sensitive:
        log1: logs/data_tables.log   
        table1: output/tables/data_table_*.csv

  summary_tables:
    run: stata-mp:latest analysis/400_summary_tables.do "{primary_disease}" "{comorbidities_list_stata}" "{disease_features_list_stata}" "{events_list_stata}" "{admissions_list_stata}" "{bloods_list_stata}" "{medications_list_stata}" "{outpatients_list_stata}"
    needs: [cohort_cleaning]
    outputs:
      moderately_sensitive:
        log1: logs/summary_tables.log   
        table1: output/tables/summary_table_*.csv

  temporal_plots:
    run: stata-mp:latest analysis/500_temporal_plots.do "{primary_disease}" "{demographic_list_stata}" "{studystart_date}" "{studyend_date}" "{studyfup_date}" "{intervention_date_2}"
    needs: [data_tables]
    outputs:
      moderately_sensitive:
        log1: logs/temporal_plots.log   
        figure1: output/figures/temporal_plot_*.svg

  logistic_models:
    run: stata-mp:latest analysis/600_logistic_models.do "{primary_disease}"
    needs: [cohort_cleaning]
    outputs:
      moderately_sensitive:
        log1: logs/logistic_models.log   
        table1: output/tables/melogit_summary.csv
        table2: output/tables/logistic_summary.csv

  survival_models:
    run: stata-mp:latest analysis/700_survival_models.do "{primary_disease}" "{studyfup_date}"
    needs: [cohort_cleaning]
    outputs:
      moderately_sensitive:
        log1: logs/survival_models.log   
        table1: output/tables/landmark_cox_summary.csv
        figure1: output/figures/km_*.svg
        figure2: output/figures/loglog_*.svg

  generate_notebook:
    run: jupyter:latest jupyter nbconvert /workspace/analysis/report.ipynb --execute --to html --template basic --output-dir=/workspace/output --ExecutePreprocessor.timeout=86400 --no-input
    needs: [temporal_plots, incidence_graphs, sarima]
    outputs:
      moderately_sensitive:
        notebook: output/report.html                                
  """

# Combine header, body, and footer
generated_yaml = yaml_header + yaml_body + yaml_incidence + yaml_practice + yaml_footer

# Save to a file
with open("project.yaml", "w") as file:
    file.write(generated_yaml)