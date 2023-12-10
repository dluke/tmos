
import os
join = os.path.join
import sys

import json
import subprocess


root = os.path.normpath(os.path.join(os.path.dirname(__file__)))

notedir = join(root, 'notebook/')

sys.path.append(os.path.join(root, 'src/pyunit'))

vpath = os.path.join(root,'version.json')

def get_version():
	meta = get_meta()
	return meta['version']

def get_meta():
	if not os.path.exists(vpath):
		print('Warning: Could not identify pypili version number in {}'.format(
			os.path.abspath(vpath)))
		return None
	with open(vpath) as fj:
		metadata = json.load(fj)
	return metadata

def update_meta():
	revparse = 'git rev-parse --short HEAD'.split()
	commit_hash = subprocess.check_output(revparse, cwd=root).decode(sys.stdout.encoding).strip()
	meta = get_meta()
	print('updating version.json')
	meta['commit'] = commit_hash
	print(json.dumps(meta, indent=1))
	with open(vpath, 'w') as fj:
		json.dump(meta, fj, indent=1)

def write_version(chdir):
	# write version.json to current dir
	with open(vpath) as f:
		meta = json.load(f)
	with open(os.path.join(chdir,'version.json'), 'w') as f:
		json.dump(meta, f)

__version__ = get_version()

os.environ['PILI_ROOT'] = os.path.dirname(os.path.realpath(__file__))
