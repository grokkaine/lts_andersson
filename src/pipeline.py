"""
Small pipeline for the purpose of running the project on UPPMAX on a
workstation and/or UPPMAX.
"""

_config_file = ""
_run_dict = {}


def downloadSRA(job):
    import preprocessing
    template_fpath = job["template_fpath"]
    run_dir = job["run_dir"]
    download_loc = job["download_loc"]
    sample_names = job["sample_names"]
    location = job["location"]
    preprocessing.SRA_download(template_fpath, run_dir,
                               download_loc, sample_names,
                               location)
    return


def merge_tsamples(job):
    with open("merge_log.txt", 'a') as f:
        f.write("Merging the time samples...\n")
    import preprocessing
    fname = job["merged_filename"]
    dir = job["sample_dir"]
    files = [dir + file for file in job["sample_filenames"]]
    preprocessing.merge_fq(fname, files, job["sample_names"])
    with open("merge_log.txt", 'a') as f:
        f.write("Done merging the time samples!\n")
    return


def read_config():
    # TODO: add sudo apt-get install python-yaml
    import yaml
    # config = yaml.load(_config_file)
    print("Reading the YAML config file...")
    with open(_config_file, 'r') as f:
        jobs = yaml.load_all(f)
        for job in jobs:
            # print("New job:")
            # for k, v in job.items(): print k, "->", v
            if job["run"]:
                # print "Running this one!"
                # for task in job["tasks"]:
                #     print "task: ", task["description"]
                #     _run_dict[task["name"]](job)
                print(job["name"])
                # call the function using a yaml supplied function name
                globals()[job["name"]](job)
            else:
                print("This doc isn't run.")
    return


def main():
    import sys
    if len(sys.argv) > 1:
        global _config_file
        _config_file = sys.argv[1]
        read_config()
    else:
        global _config_file
        _config_file = "configure.yaml"
        read_config()


if __name__ == "__main__":
    # execute only if run as a script
    main()
