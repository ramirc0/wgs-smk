rule scatter_interval_list:
    input:
        interval_list = config.input_files.calling_interval_list
    output:
        outdir = directory(f"{config.results_dir}/scatter_interval_list"),
        scattered_intervals = expand(
            f"{config.results_dir}/scatter_interval_list/temp_{{scatter:04d}}_of_{config.scatter_interval_list.scatter_count}/{{scatter}}scattered.interval_list",
            scatter = list(range(1, config.scatter_interval_list.scatter_count + 1))
        )
    params:
        scatter_count = config.scatter_interval_list.scatter_count,
        break_bands_at_multiples_of = config.scatter_interval_list.break_bands_at_multiples_of
    software: SOFTWARE_ENV
    log: f"{config.log_dir}/scatter_interval_list.log"
    shell:
        """
        exec >> {log} 2>&1
        python workflow/rules/scripts/scatter_interval_list.py \
            {resources.mem_mb} \
            {params.scatter_count} \
            {params.break_bands_at_multiples_of} \
            {input.interval_list} \
            {output.outdir}
        """
