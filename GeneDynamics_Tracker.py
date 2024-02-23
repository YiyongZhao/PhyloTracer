import subprocess
software_path = "./software/dolloparsimony"  


mcl_file = "input.mcl"
other_parameters = ["--param1", "value1", "--param2", "value2"]  


command = [software_path, mcl_file] + other_parameters

try:
    subprocess.run(command, check=True)
    print("Program executed successfully！")
except subprocess.CalledProcessError as e:
    print("Program execution failed!", e)
