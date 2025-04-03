import argparse
import matplotlib.pyplot as plt
import numpy as np
import subprocess
from typing import Tuple

plt.rcParams['font.size'] = 14.0
plt.rcParams['legend.fontsize'] = 12.0

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--order', type=int, required=True, action='append')
parser.add_argument('-cs', '--use-compensated-sums', action='store_true')
args = parser.parse_args()

if args.use_compensated_sums:
  compensated_sums = '--use-compensated-sums'
else:
  compensated_sums = ''

dt0 = np.pi/2.0
N = 6
tfinal = 6.283185

def get_error(output: str) -> float:
  lines = output.split('\n')
  for line in lines:
    if not line.startswith('t'):
      continue
    time_field, solution_field, energy_field, error_field = line.split(',')
    time = float(time_field.replace('t = ', ''))
    if time == tfinal:
      error = float(error_field.replace('sol. err = ',''))
      return error
  raise Exception()

def get_evals(output: str) -> Tuple[int,int]:
  lines = output.split('\n')
  evals = [0,0]
  for i in range(2):
    _, eval_field = lines[-3+i].split(' = ')
    evals[i] = int(eval_field)

  return evals

f, ax = plt.subplots(1,4,sharey=True,figsize=(16,8))
color_dict = {}
for order in args.order:
  dt_list = []
  error_list = []
  evals_list = []
  for n in range(N):
    dt = dt0 / 2**n
    exe = f'./ark_harmonic_symplectic --nout 1 --dt {dt} --order {order} {compensated_sums}'
    print(f'Executing {exe}...')
    process = subprocess.run(exe, shell=True, check=True, capture_output=True,
      text=True)
    error = get_error(process.stdout)
    evals = get_evals(process.stdout)
    dt_list.append(dt)
    error_list.append(error)
    evals_list.append(evals)

  dt_list = np.array(dt_list)
  error_list = np.array(error_list)
  evals_list = np.array(evals_list)
  if order > 0:
    linespec = '-'
    label='SPRK'
  else:
    linespec = '--o'
    label='Split'
  if abs(order) in color_dict.keys():
    color = color_dict[abs(order)]
    ax[0].loglog(dt_list, error_list, linespec, color=color)
  else:
    line, = ax[0].loglog(dt_list, error_list, linespec)
    color = color_dict[abs(order)] = line.get_color()
  ax[1].loglog(np.sum(evals_list, axis=1), error_list, linespec, color=color, label=f'{label}{abs(order)}')
  ax[2].loglog(evals_list[:,0], error_list, linespec, color=color)
  ax[3].loglog(evals_list[:,1], error_list, linespec, color=color)
  if order > 5:
    ax[0].loglog(dt_list[:3], error_list[0]*(dt_list[:3]/dt_list[0])**abs(order), '--k')
  elif order > 0:
    ax[0].loglog(dt_list[-3:], error_list[-1]*(dt_list[-3:]/dt_list[-1])**abs(order), '--k')

ax[0].set_xlabel('dt')
ax[1].set_xlabel('# RHS evals')
ax[2].set_xlabel('# f1 evals')
ax[3].set_xlabel('# f2 evals')
ax[0].set_ylabel('RMSE')
ax[1].legend(loc='best')
f.tight_layout()
plt.show()



