import sys
import subprocess

def run_dolloparsimony(file1, file2, file3, output_file):
    command_list = ["./software/dolloparsimony" , file1, file2, file3, output_file]
    try:
        subprocess.run(command_list, check=True)
        print("Program executed successfully！")
    except subprocess.CalledProcessError as e:
        print("Program execution failed:", e)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python GeneDynamics_Tracker.py file1 file2 file3 output_file")
        sys.exit(1)
    
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    file3 = sys.argv[3]
    output_file = sys.argv[4]
    
    run_dolloparsimony(file1, file2, file3, output_file)
