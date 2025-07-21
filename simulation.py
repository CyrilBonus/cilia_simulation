import subprocess
subprocess.run(['bash', 'dot_dat_cleanup.sh'])
subprocess.run(['make', 'clean'])
subprocess.run(['make'])
subprocess.run(['python', 'output.py'])
print("what should be the limit of dat files to generate? (default 10)")
limit = input()
if not limit.isdigit():
	limit = 10
subprocess.run(['./simulation' + ' ' + limit], shell=True)
subprocess.run(['python', 'interactive_heatmap.py'])
