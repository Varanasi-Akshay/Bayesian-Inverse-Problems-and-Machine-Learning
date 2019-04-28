import numpy as np
import sys
from os.path import dirname, realpath
parent = dirname(dirname(realpath(__file__)))
parent+='/sampleGenerator'
sys.path.append(parent)
from sampleGenerator import sampleGenerator
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import torch
from potentialNet import potential_net, cal_stiff, cal_stress, scaler_E, scaler_S
import numpy as np
import matplotlib.pyplot as plt
import IPython

print("start loading trained net")

potential_net.load_state_dict(torch.load('trained_net.pt'))

# n = 21
# E, S = sampleGenerator(n = n)

# E = scaler_E.transform(E)
# S = scaler_S.transform(S)

# E, S = torch.from_numpy(E), torch.from_numpy(S)
# E.requires_grad = True

# num_of_samples = len(E)





# plot biaxial traning test
n = 41
E, S = sampleGenerator(
    n = (n, n, 1)
)

# # plot biaxial validation test
# n = 101
# E, S = sampleGenerator(
#     Erange = (0.23, 0.29, 0.0),
#     n = (n, n, 1)
# )

I = [i+n*i for i in range(n)]
E, S = E[I], S[I]

E_np_transformed = scaler_E.transform(E)

E_torch_transformed = torch.from_numpy(E_np_transformed)
E_torch_transformed.requires_grad = True


Sp = []
Cp = []
detp = []

sample_size = len(E)

for i in range(sample_size):
    x = E_torch_transformed[i]
    fp = potential_net(x)
    sp = cal_stress(fp, x)
    cp = cal_stiff(sp, x) # need inverse transform
    dp = np.linalg.det(cp.data.numpy())

    sp = sp.data.numpy()
    cp = cp.data.numpy()
    sp = scaler_S.inverse_transform(sp)

    #TODO figure out a way to transform cp.
    Sp += [sp]
    Cp += [cp]
    detp += [dp]

Sp = np.array(Sp)
Cp = np.array(Cp)

eigp = []
for i in range(sample_size):
    _, s, _ = np.linalg.svd(Cp[i])
    eigp += [s]
eigp = np.array(eigp)

# # plot with data only
# plt.plot(E[:,0], S[:,0], '-', E[:,1], S[:,1], '-')
# plt.ylabel('S')
# plt.xlabel('E')
# plt.title('Biaxial Test: training data')
# plt.show()

# plot with data and curve
plt.plot(E[:,0], S[:,0], 'o', E[:,0], Sp[:,0], '-', E[:,1], S[:,1], 'o', E[:,1], Sp[:,1], '-')
plt.ylabel('S')
plt.xlabel('E')
plt.title('Biaxial Test: training data')
plt.show()

# # plot with validation data only
# J=E[:,0]>0.21
# plt.plot(E[J,0], S[J,0], 'o', E[J,0], Sp[J,0], '-', E[J,1], S[J,1], 'ro', E[J,1], Sp[J,1], '-')
# plt.ylabel('S')
# plt.xlabel('E')
# plt.title('Biaxial Test: validation data')
# plt.show()

IPython.embed()