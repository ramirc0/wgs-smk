# TODO: Add some kind of check to return a boolean if duplication or chimerism values exceed these thresholds:
# Float max_duplication_in_reasonable_sample = 0.30
# Float max_chimerism_in_reasonable_sample = 0.15

rule check_prevalidation:
    input:
        duplication_metrics = f"{config.results_dir}/mark_duplicates/{{sample}}.metrics",
        agg_alignment_summary_metrics = f"{config.results_dir}/collect_aggregation_metrics/{{sample}}.alignment_summary_metrics",
    output:
        duplication_value = f"{config.results_dir}/check_prevalidation/{{sample}}.duplication_value.txt",
        chimerism_value = f"{config.results_dir}/check_prevalidation/{{sample}}.chimerism_value.txt"
    params:
        tmp_dir = config.tmp_dir
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/check_prevalidation/{{sample}}.log"
    run:
        import sys
        import csv
        import subprocess
        import uuid
        import tempfile
        from pathlib import Path
        
        sys.stderr = open(log[0], "w")

        duplication_csv = Path(f"{params.tmp_dir}/{str(uuid.uuid4())}.csv") 
        chimerism_csv = Path(f"{params.tmp_dir}/{str(uuid.uuid4())}.csv")

        grep_dup_met = f"grep -A 1 PERCENT_DUPLICATION {input.duplication_metrics} > {duplication_csv.resolve()}"
        grep_align_met = f"grep -A 3 PCT_CHIMERAS {input.agg_alignment_summary_metrics} | grep -v OF_PAIR > {chimerism_csv.resolve()}"
        subprocess.check_call(grep_dup_met, shell=True)
        subprocess.check_call(grep_align_met, shell=True)
 
        with duplication_csv.open("r") as dupfile:
            reader = csv.DictReader(dupfile, delimiter='\t')
            for row in reader:
                with open(output.duplication_value, "w") as file:
                    file.write(row['PERCENT_DUPLICATION'])
                    file.close()

        with chimerism_csv.open("r") as chimfile:
            reader = csv.DictReader(chimfile, delimiter='\t')
            for row in reader:
                with open(output.chimerism_value, "w") as file:
                    file.write(row['PCT_CHIMERAS'])
                    file.close()
