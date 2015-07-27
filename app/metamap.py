import json, subprocess

textfile = 'madeupname.txt'

p = subprocess.Popen(['java', '-jar', 'MetaMapCaller.jar', textfile], 
	cwd='./java/lib/', stdout=subprocess.PIPE)
# JSON obj!
print(p.stdout.readline().decode('UTF-8'))