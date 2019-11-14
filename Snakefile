PYTHON = "PYTHONPATH=./ python"
PANDOC = "pandoc --filter pantable --filter pandoc-fignos --filter pandoc-tablenos --filter pandoc-citeproc"
PYTHON_SCRIPT = "PYTHONPATH=./ python {input} {output}"
PYTHON_SCRIPT_WITH_CONFIG = PYTHON_SCRIPT + " {CONFIG_FILE}"

CONFIG_FILE = "config/default.yaml"

configfile: CONFIG_FILE
include: "rules/data-preprocessing.smk"
include: "rules/sonnendach.smk"
include: "rules/capacityfactors.smk"
include: "rules/potential.smk"
include: "rules/analysis.smk"
include: "rules/publish-data.smk"
include: "rules/sync.smk"

localrules: all, paper, supplementary_material, clean

wildcard_constraints:
    layer = "({layer_list})".format(layer_list="|".join((f"({layer})" for layer in config["layers"]))),
    scenario = "({scenario_list})".format(scenario_list="|".join((f"({scenario})" for scenario in config["scenarios"])))

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'out-of-sight succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'out-of-sight crashed' {config[email]}")


rule all:
    message: "Run entire analysis and compile report."
    input:
        rules.wind_capacity_per_distance_plot.output,
        rules.wind_capacity_per_distance_map.output


GENERAL_DOCUMENT_DEPENDENCIES = [
    "report/literature.bib",
    "report/pandoc-metadata.yml",
    "report/energy-strategy-reviews.csl",
    "report/template.html",
    "report/report.css",
    "report/fonts/KlinicSlabBook.otf",
    "report/fonts/KlinicSlabBookIt.otf",
    "report/fonts/KlinicSlabMedium.otf",
    "report/fonts/KlinicSlabMediumIt.otf",
]


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--self-contained --css=report.css --number-sections --template template.html --to html5"
    elif suffix == "pdf":
        return "--css=report.css --number-sections --template template.html --pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule paper:
    message: "Compile {output}."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/paper.md",
        "build/technical-potential/normed-potentials-boxplots.png",
        "build/technical-potential/sufficient-potentials-map.png",
        "build/technical-social-potential/sufficient-potentials-map.png",
        "build/technical-social-potential/normed-potentials-boxplots.png",
        "build/necessary-land/overview-necessary-land-when-pv-100%.csv",
        "build/necessary-land/overview-necessary-land-when-pv-40%.csv",
        "build/necessary-land/necessary-land-map-when-pv-40%.png",
        "build/necessary-land/necessary-land-all-layers.png",
        "build/exclusion-layers-ROU.png",
        rules.layer_overview.output
    params: options = pandoc_options
    output: "build/paper-{dpi}dpi.{suffix}"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        rsync -rpLtgoD ../build .
        if [ {wildcards.dpi} != 600 ];then
            echo "Resizing images."
            for i in {input}; do
                if [[ $i == *.png ]];then
                    echo "Resizing $i."
                    convert -resample {wildcards.dpi}x{wildcards.dpi} -units PixelsPerInch $i $i
                fi
            done
        fi
        {PANDOC} paper.md pandoc-metadata.yml {params.options} -o ../build/paper-{wildcards.dpi}dpi.{wildcards.suffix}
        """


rule supplementary_material:
    message: "Compile the supplementary material."
    input:
        GENERAL_DOCUMENT_DEPENDENCIES,
        "report/supplementary.md",
        "report/data-sources.csv",
        rules.sonnendach_statistics.output.publish,
        expand(
            "build/exclusion-layers-{country_code}.png",
            country_code=[pycountry.countries.lookup(country).alpha_3
                          for country in config["scope"]["countries"]]
        ),
        "build/national/technical-potential/potentials-polished.csv",
        "build/national/technical-social-potential/potentials-polished.csv"
    params: options = pandoc_options
    output: "build/supplementary-material.{suffix}"
    shadow: "minimal"
    conda: "envs/report.yaml"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} supplementary.md {params.options} -o {output} --table-of-contents
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r ./build/*
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    message: "Run tests."
    input:
        expand("build/{layer}/technical-potential/potentials.csv", layer=config["layers"]),
        "build/technically-eligible-land.tif",
        "build/technically-eligible-area-km2.tif",
        "build/technically-eligible-electricity-yield-pv-prio-twh.tif",
        "build/administrative-borders-nuts.gpkg",
        "build/swiss/total-rooftop-area-according-to-sonnendach-data-km2.txt",
        "build/swiss/total-yield-according-to-sonnendach-data-twh.txt"
    output: "build/logs/test-report.html"
    conda: "envs/default.yaml"
    shell:
        "py.test --html={output} --self-contained-html"
