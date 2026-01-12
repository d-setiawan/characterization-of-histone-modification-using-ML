import subprocess

# List of scripts in order
scripts = [
    {"script": "1) Pre-Processing.py", "input": None, "output": "pre_processing.txt"},
    {"script": "2) EDA.py", "input": None, "output": "EDA.txt"},
    {"script": "3) Prelim_Analysis.py", "input": None, "output": "PreAnalysis.txt"},
    {"script": "5) KNNR_testing.py", "input": None, "output": "KNNR.txt"},
    {"script": "5) LightGBM_testing.py", "input": None, "output": "LightGBM.txt"},
    {"script": "5) MLP_testing.py", "input": None, "output": "MLP.txt"},
    {"script": "5) SVR_testing.py", "input": None, "output": "SVR.txt"},
    {"script": "6) Optuna_Optimization.py", "input": None, "output": "Optuna.txt"},
    {"script": "7) Feature_Validation.py", "input": None, "output": "final_report.txt"},
]

def run_pipeline(scripts):
    for i, step in enumerate(scripts):
        script = step["script"]
        input_file = step["input"]
        output_file = step["output"]

        print(f"Running Step {i+1}: {script}...")

        try:
            # Build the command
            command = ["python", script]
            if input_file:
                command += ["--input", input_file]
            if output_file:
                command += ["--output", output_file]

            # Redirect stdout and stderr to the output file
            with open(output_file, "w") as log_file:
                subprocess.run(command, stdout=log_file, stderr=log_file, check=True)

        except subprocess.CalledProcessError as e:
            print(f"Error in {script}: {e}")
            return False

    return True

if __name__ == "__main__":
    success = run_pipeline(scripts)
    if success:
        print("Pipeline completed successfully!")
    else:
        print("Pipeline failed.")
