import torch
import torch.nn as nn
import torch.nn.functional as functional
import torch.optim as optim
import sys
from os.path import dirname, realpath
parent = dirname(dirname(realpath(__file__)))
parent+='/sampleGenerator'
sys.path.append(parent)
from sampleGenerator import sampleGenerator
from sklearn.preprocessing import StandardScaler
from stressNet import stress_net, scaler_E, scaler_S, cal_stiff
import IPython

"""
input-strains-output-stress neural network

[1142] loss:     19.508
    rloss: 0.00528005
"""

# Customize the loss func
# https://spandan-madan.github.io/A-Collection-of-important-tasks-in-pytorch/


e, s = sampleGenerator()

E = scaler_E.transform(e)
S = scaler_S.transform(s)

E, S = torch.from_numpy(E), torch.from_numpy(S)

E.requires_grad = True
# S.requires_grad = True

criterion = nn.MSELoss()
optimizer = optim.Adam(stress_net.parameters(), amsgrad=True)

running_loss = 0.0

num_of_samples = len(E)

epoch = 0
rloss = 1.0
init_loss = 0.0
rloss_old = rloss

filename='trained_net.pt'
while rloss > 0.01 or running_loss < 1000:  # loop over the dataset multiple times
    running_loss = 0.0
    for i in range(num_of_samples):
        # zero the parameter gradients
        optimizer.zero_grad()

        # forward + backward + optimize
        x = E[i]
        outputs = stress_net(x)
        C = cal_stiff(outputs, x)

        #TODO: grad(outputs) <= outputs
        loss = criterion(outputs, S[i])
        loss += criterion(C, 1./2*(C + torch.t(C)))
        loss.backward()
        optimizer.step()

        # print statistics
        running_loss += loss.item()
        if epoch == 1:
            init_loss += criterion(
                S[i], 
                torch.zeros_like(S[i])
            ).item()

    if epoch >= 1:
        rloss = running_loss / init_loss
    print('[%d] loss: %10.3f' %(epoch + 1, running_loss))
    print('    rloss: %10.8f' %(rloss))
    print('  descent: '+str(rloss < rloss_old))
    rloss_old = rloss
    epoch += 1
    torch.save(stress_net.state_dict(),filename)


print('Finished Training')

IPython.embed()
