from pele_platform.analysis import Analysis

analysis = Analysis(resname="LIG", chain="L", simulation_output="output", report="mod_report", cpus=48)
analysis.generate(path="analysis", clustering_type="meanshift")
